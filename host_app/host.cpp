#include <numeric>
#include <random>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <thread>
#include <future>

#include "include/OclApi.h"
#include "include/Types.h"
#include "include/DataSetup.h"
#include "include/Permutation.h"

int main(int argc, char **argv)
{
    unsigned long long int block_size_arg = 0;
    unsigned long long int num_snps = 0;
    unsigned int num_indiv = 0;
    unsigned int num_perms = 0;
    unsigned int min_perms_adp = 0;

    double phenotype_scale_factor = 1;
    perm_algo_t perm_algo = perm_regression_only;

    std::string genotype_file;
    std::string phenotype_file;
    std::string bim_file;

    std::string binaryFile;

    for (size_t i = 0; i < argc; i++)
    {
        std::string arg = argv[i];

        if ( arg == "-xclbin" )
        {
            std::cout << std::endl << "FPGA binary: " << argv[i + 1] << std::endl;
            binaryFile = argv[i + 1];
        }

        if ( arg == "-phenotype_scale_factor" )
        {
            phenotype_scale_factor = std::stod(argv[i + 1]);
            std::cout << "Phenotype scaling factor: " << phenotype_scale_factor << std::endl;
        }

        if ( arg == "-input_data" )
        {
            genotype_file = argv[i + 1];
            phenotype_file = argv[i + 2];
            bim_file = argv[i + 3];

            std::cout << "Genotype file: " << genotype_file << std::endl;
            std::cout << "Phenotype file: " << phenotype_file << std::endl;
            std::cout << "BIM file: " << bim_file << std::endl;
        }

        if ( arg == "-num_rows" )
        {
            num_snps = std::stoi( argv[i + 1] );
        }

        if ( arg == "-perm_algo" )
        {
            std::string algo_string = argv[i + 1];

            if (algo_string.compare("adp") == 0)
            {
                perm_algo = perm_adaptive;

                min_perms_adp = std::stoi( argv[i + 2] );
                num_perms = std::stoi( argv[i + 3] );

                std::cout << "Adaptive permutation: " << std::endl
                          << "Min perms per SNP = " << min_perms_adp << std::endl
                          << "Max perms per SNP = " << num_perms << std::endl;

                block_size_arg = 8388608;
            }
            else if (algo_string.compare("maxT") == 0)
            {
                perm_algo = perm_maxT;

                num_perms = std::stoi( argv[i + 2] );

                std::cout << "Performing " << num_perms << " maxT permutations" << std::endl;

                block_size_arg = 33554432;
            }
        }

        if ( arg == "-buf_size" )
        {
            block_size_arg = std::stoull( argv[i + 1] );
            std::cout << "Set buffer size to " << block_size_arg << std::endl;
        }
    }

    if ( binaryFile.empty() )
    {
        std::cout << "Usage: " << argv[0] << " -xclbin <XCLBIN File>" << std::endl;
        return EXIT_FAILURE;
    }

    if ( perm_algo == perm_regression_only )
    {
        std::cout << "Permutation testing option not selected, running GWAS analysis" << std::endl;
    }

    auto time = std::chrono::system_clock::now();
    auto time_out = std::chrono::system_clock::to_time_t(time);

    std::cout << std::endl << "Start time: " << std::ctime(&time_out) << std::endl;

    std::vector<double> pheno_file_data;
    std::vector<char> bed_file_data;
    size_t num_indiv_from_pheno_file = 0;
    size_t num_snps_from_bed_file = 0;

    auto start_import_time = std::chrono::high_resolution_clock::now();

    std::cout << "Reading phenotype file " << phenotype_file << std::endl;
    std::ifstream pheno_file_stream(phenotype_file);
    if (pheno_file_stream.is_open())
    {
        std::string line;
        while (std::getline(pheno_file_stream, line))
        {
            double num = std::stod(line);
            pheno_file_data.push_back(num);
        }
        pheno_file_stream.close();

        num_indiv_from_pheno_file = pheno_file_data.size();
    }
    std::cout << "Read phenotypes of " << num_indiv_from_pheno_file << " individuals" << std::endl;

    std::cout << "Reading .bim file " << bim_file << std::endl;
    std::ifstream inFile(bim_file);
    size_t num_snps_from_bim_file = std::count(std::istreambuf_iterator<char>(inFile),
                                                std::istreambuf_iterator<char>(), '\n');
    std::cout << "Read " << num_snps_from_bim_file << " SNPs from " << bim_file << std::endl;

    size_t expected_file_size = (num_snps_from_bim_file * ceil(num_indiv_from_pheno_file/4.0)) + 3;

    std::cout << "Reading .bed file " << genotype_file << std::endl;
    std::ifstream bed_file(genotype_file, std::ios::in | std::ios::binary);
    bed_file.seekg(0, std::ios::end);
    size_t genotype_filesize = bed_file.tellg();

    if ( genotype_filesize != expected_file_size )
    {
        std::cout << "Error: Expected .bed file size " << expected_file_size << std::endl;
        return 1;
    }

    //Skip magic number at the beginning of the .bed file
    bed_file.seekg(3, std::ios::beg);

    // Allocate vector to store genotype data
    bed_file_data.resize(genotype_filesize-3);
    bed_file.read((char *)bed_file_data.data(), genotype_filesize);

    num_snps_from_bed_file = (genotype_filesize-3) / ceil(num_indiv_from_pheno_file/4.0);
    std::cout << "Read " << num_snps_from_bed_file << " SNPs for " << num_indiv_from_pheno_file << " individuals" << std::endl;

    num_indiv = num_indiv_from_pheno_file;
    if (num_snps == 0)
    {
        num_snps = num_snps_from_bed_file;
    }

    std::cout << std::endl << "Analysing " << num_snps << " SNPs" << std::endl;

    // Buffer setup based on the matrix block size and the number of kernels
    unsigned int num_cu = 4;
    std::cout << "Using " << num_cu << " kernels" << std::endl;

    size_t matrix_block_size_bytes = ( block_size_arg > 0 ) ? block_size_arg : (4194304 * 2);

    size_t bytes_per_snp = round_to_multiple((size_t)ceil(num_indiv / 4.0), (INPUT_MAT_PAR_ENTRIES / 4));
    size_t snps_per_block = matrix_block_size_bytes / bytes_per_snp;
    size_t num_matrix_blocks = ceil( (double)num_snps / snps_per_block );
    size_t last_block_snps = num_snps - ( ( num_matrix_blocks - 1 ) * snps_per_block );

    size_t input_matrix_size_bytes = num_matrix_blocks * matrix_block_size_bytes;

    size_t pheno_num_elems = bytes_per_snp * 4;
    size_t pheno_size_bytes = INPUT_VEC_DATATYPE_SIZE_BYTES * pheno_num_elems;

    size_t in_mean_bytes_per_block = round_to_multiple( ( INPUT_MEAN_DATATYPE_SIZE_BYTES * snps_per_block ), 4096 );
    size_t in_mean_elems_per_block = in_mean_bytes_per_block / INPUT_MEAN_DATATYPE_SIZE_BYTES;
    size_t in_mean_elems = num_matrix_blocks * in_mean_elems_per_block;

    size_t out_vec_bytes_per_block = round_to_multiple( ( OUTPUT_DATATYPE_SIZE_BYTES * snps_per_block ), 4096 );
    size_t out_vec_elems_per_block = out_vec_bytes_per_block / OUTPUT_DATATYPE_SIZE_BYTES;
    size_t out_vec_elems = num_matrix_blocks * out_vec_elems_per_block;

    std::cout << "Bytes per SNP = " << bytes_per_snp << " bytes" << std::endl;
    std::cout << "Matrix block size = " << matrix_block_size_bytes << " bytes" << std::endl;
    std::cout << "Num matrix blocks = " << num_matrix_blocks << std::endl;
    std::cout << "SNPs per block = " << snps_per_block << std::endl;
    std::cout << "Last block SNPs = " << last_block_snps << std::endl;

    std::cout << "Matrix Size: \t\t" << input_matrix_size_bytes << " bytes" << std::endl;
    std::cout << "Phenotype Vector Size: \t" << pheno_size_bytes << " bytes" << std::endl;
    std::cout << "Input Mean Size: \t" << (num_matrix_blocks * in_mean_bytes_per_block) << " bytes" << std::endl;
    std::cout << "Output Vector Size: \t" << (num_matrix_blocks * out_vec_bytes_per_block) << " bytes" << std::endl;

    std::cout << "Input mean buffer size per block: \t" << in_mean_elems_per_block << "\t\t"
              << in_mean_bytes_per_block << " bytes" << std::endl;

    std::cout << "Input vector buffer size per block: \t" << pheno_num_elems << "\t\t"
              << pheno_size_bytes << " bytes" << std::endl;

    std::cout << "Output vector buffer size per block: \t" << out_vec_elems_per_block << "\t\t"
              << out_vec_bytes_per_block << " bytes" << std::endl;

    input_mean_vector_t     input_mean_buffer(in_mean_elems, 0);
    input_matrix_vector_t   input_matrix_buffer(input_matrix_size_bytes, 0);
    input_pheno_vector_t    input_pheno_buffer(pheno_num_elems, 0);
    output_vector_t         fpga_output_buffer(out_vec_elems, 0);

    vector_stats_t vector_stats(num_snps);

    vector_stats = generate_fpga_data(num_snps, num_indiv,
                                      matrix_block_size_bytes,
                                      bed_file_data, pheno_file_data,
                                      phenotype_scale_factor,
                                      input_mean_buffer, input_matrix_buffer, input_pheno_buffer);

    input_mean_vector_t *input_mean_ptr = &input_mean_buffer;
    input_matrix_vector_t *input_matrix_ptr = &input_matrix_buffer;

    std::cout << "Data import complete" << std::endl << std::endl;

    auto end_import_time = std::chrono::high_resolution_clock::now();

    // Observed results vector
    std::vector<double> gwas_result(num_snps, 0);

    // For maxT permutation
    std::vector<double> perm_maxF_res_vector(num_perms + 1);

    // For adaptive permutation
    int perm_adp_min_perms_per_snp = min_perms_adp;

    size_t perm_adp_drop_interval = perm_adp_min_perms_per_snp;
    size_t perm_adp_drop_counter = 0;
    size_t perm_adp_one_perm_per_cu_start = 0;

    perm_adp_results_t perm_adp_results(num_snps);

    // Buffers to store regenerated data for adaptive permutation
    input_mean_vector_t input_mean_adp_buffer;
    input_matrix_vector_t input_matrix_adp_buffer;
    if (perm_algo == perm_adaptive)
    {
        input_mean_adp_buffer.resize(in_mean_elems, 0);
        input_matrix_adp_buffer.resize(input_matrix_size_bytes, 0);
    }

    auto start_processing_time = std::chrono::high_resolution_clock::now();

    // Seed RNG for std::random_shuffle
    std::srand(std::time(NULL));

    // Open CL setup
    OclApi api(binaryFile);
    cl_int ocl_error_code;

    std::vector<ocl_task_buffers_t> kernel_buffers(num_cu);
    std::vector<OclTask> ocl_tasks(num_matrix_blocks);

    std::vector<cl::Kernel> kernels(num_cu);

    int nkernel = 0;
    OCL_CHECK(ocl_error_code, kernels[nkernel++] = cl::Kernel(api.program, "krnl_mvmul_0", &ocl_error_code));
    OCL_CHECK(ocl_error_code, kernels[nkernel++] = cl::Kernel(api.program, "krnl_mvmul_1", &ocl_error_code));
    OCL_CHECK(ocl_error_code, kernels[nkernel++] = cl::Kernel(api.program, "krnl_mvmul_2", &ocl_error_code));
    OCL_CHECK(ocl_error_code, kernels[nkernel++] = cl::Kernel(api.program, "krnl_mvmul_3", &ocl_error_code));

    for (size_t perm = 0; perm < (num_perms + 1); perm++)
    {
        if (perm > 0)
        {
            std::random_shuffle(input_pheno_buffer.begin(), input_pheno_buffer.begin() + num_indiv);
        }

        for (unsigned int i = 0; i < num_matrix_blocks; i++)
        {
            int kernel_num = i % num_cu;

            OCL_CHECK(ocl_error_code,
                      kernel_buffers[kernel_num].mat_buf =
                          cl::Buffer(api.context,
                                     CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                                     matrix_block_size_bytes,
                                     input_matrix_ptr->data() + (i * matrix_block_size_bytes),
                                     &ocl_error_code));

            OCL_CHECK(ocl_error_code,
                      kernel_buffers[kernel_num].mean_buf =
                          cl::Buffer(api.context,
                                     CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                                     in_mean_bytes_per_block,
                                     input_mean_ptr->data() + (i * in_mean_elems_per_block),
                                     &ocl_error_code));

            OCL_CHECK(ocl_error_code,
                      kernel_buffers[kernel_num].pheno_buf =
                          cl::Buffer(api.context,
                                     CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                                     pheno_size_bytes,
                                     input_pheno_buffer.data(),
                                     &ocl_error_code));

            OCL_CHECK(ocl_error_code,
                      kernel_buffers[kernel_num].out_buf =
                          cl::Buffer(api.context,
                                     CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY,
                                     out_vec_bytes_per_block,
                                     fpga_output_buffer.data() + (i * out_vec_elems_per_block),
                                     &ocl_error_code));

            unsigned int n_rows = (i == num_matrix_blocks - 1) ? last_block_snps : snps_per_block;

            int narg = 0;
            OCL_CHECK(ocl_error_code, ocl_error_code = kernels[kernel_num].setArg(narg++, (const unsigned int)n_rows));
            OCL_CHECK(ocl_error_code, ocl_error_code = kernels[kernel_num].setArg(narg++, (const unsigned int)pheno_num_elems));
            OCL_CHECK(ocl_error_code, ocl_error_code = kernels[kernel_num].setArg(narg++, kernel_buffers[kernel_num].mean_buf));
            OCL_CHECK(ocl_error_code, ocl_error_code = kernels[kernel_num].setArg(narg++, kernel_buffers[kernel_num].mat_buf));
            OCL_CHECK(ocl_error_code, ocl_error_code = kernels[kernel_num].setArg(narg++, kernel_buffers[kernel_num].pheno_buf));
            OCL_CHECK(ocl_error_code, ocl_error_code = kernels[kernel_num].setArg(narg++, kernel_buffers[kernel_num].out_buf));

            if (i < num_cu)
            {
                ocl_tasks[i].run(api, kernels[kernel_num], kernel_buffers[kernel_num]);
            }
            else
            {
                ocl_tasks[i].run(api, kernels[kernel_num], kernel_buffers[kernel_num], ocl_tasks[i - num_cu].getDoneEv());
            }
        }

        OCL_CHECK(ocl_error_code, ocl_error_code = api.queue.finish());

        if (perm_algo == perm_regression_only)
        {
            break;
        }

        // Update the permutation results
        if (perm == 0)
        {
            if (perm_algo == perm_adaptive)
            {
                perm_adp_update_results(perm,
                                        num_snps,
                                        snps_per_block,
                                        out_vec_elems_per_block,
                                        std::cref(fpga_output_buffer),
                                        std::cref(vector_stats),
                                        std::ref(gwas_result),
                                        std::ref(perm_adp_results));
            }
            else if (perm_algo == perm_maxT)
            {
                perm_maxT_update_results(perm,
                                         num_snps,
                                         snps_per_block,
                                         out_vec_elems_per_block,
                                         std::cref(fpga_output_buffer),
                                         std::cref(vector_stats),
                                         std::ref(gwas_result),
                                         std::ref(perm_maxF_res_vector));
            }
        }
        else
        {
            if (perm_algo == perm_adaptive)
            {
                std::thread result_thread(perm_adp_update_results,
                            perm,
                            num_snps,
                            snps_per_block,
                            out_vec_elems_per_block,
                            fpga_output_buffer,
                            std::cref(vector_stats),
                            std::ref(gwas_result),
                            std::ref(perm_adp_results));

                if ( ( perm % perm_adp_drop_interval == 0 ) )
                {
                    result_thread.join();
                }
                else
                {
                    result_thread.detach();
                }

                if ( perm % perm_adp_drop_interval == 0 )
                {
                    size_t num_dropped_snps = perm_adp_drop_snps(perm,
                                                                 num_snps,
                                                                 perm_adp_min_perms_per_snp,
                                                                 perm_adp_results);

                    size_t new_num_snps = num_snps - num_dropped_snps;
                    size_t new_num_blocks = ceil((double)new_num_snps / snps_per_block);

                    size_t new_input_matrix_size = new_num_blocks * matrix_block_size_bytes;
                    size_t new_last_block_rows = new_num_snps - ((new_num_blocks - 1) * snps_per_block);

                    input_matrix_adp_buffer.resize(new_input_matrix_size);

                    perm_adp_regen_input_data(bytes_per_snp,
                                              snps_per_block,
                                              matrix_block_size_bytes,
                                              snps_per_block,
                                              matrix_block_size_bytes,
                                              input_mean_buffer,
                                              input_matrix_buffer,
                                              input_mean_adp_buffer,
                                              input_matrix_adp_buffer,
                                              perm_adp_results);

                    input_matrix_ptr = &input_matrix_adp_buffer;
                    input_mean_ptr = &input_mean_adp_buffer;

                    num_matrix_blocks = new_num_blocks;
                    last_block_snps = new_last_block_rows;

                    std::cout << "Perm = " << perm
                              << " Dropped SNPs = " << num_dropped_snps
                              << " Remaining SNPs = " << (num_snps - num_dropped_snps)
                              << " Num blocks = " << new_num_blocks << std::endl;

                    // If the number of blocks is less than the number of kernels transition to
                    // one permutation per kernel
                    if (new_num_blocks <= num_cu)
                    {
                        perm_adp_one_perm_per_cu_start = perm;
                        perm_adp_drop_counter = 0;

                        last_block_snps = (num_snps - num_dropped_snps);

                        std::cout << "Starting 1 perm per kernel with " << last_block_snps << " SNPs" << std::endl;
                        break;
                    }

                    if (perm == perm_adp_min_perms_per_snp)
                    {
                        perm_adp_drop_interval = perm_adp_min_perms_per_snp / 2;
                    }

                    // Every 5 drops increase the drop interval by 20%
                    // These parameters can probably be optimised
                    if (perm_adp_drop_counter == 5)
                    {
                        perm_adp_drop_interval *= 1.2;
                        perm_adp_drop_counter = 0;

                        std::cout << "New drop interval = " << perm_adp_drop_interval << std::endl;
                    }

                    perm_adp_drop_counter++;
                }
            }
            else if (perm_algo == perm_maxT)
            {
                std::thread result_thread(perm_maxT_update_results,
                                          perm,
                                          num_snps,
                                          snps_per_block,
                                          out_vec_elems_per_block,
                                          fpga_output_buffer,
                                          std::cref(vector_stats),
                                          std::ref(gwas_result),
                                          std::ref(perm_maxF_res_vector) );

                if ( perm == num_perms )
                {
                    result_thread.join();
                }
                else
                {
                    result_thread.detach();
                }

                if (perm % 10 == 0)
                {
                    std::cout << "\rProgress: " << ((perm * 100) / num_perms) << " %" << std::flush;
                }
            }
        }
    }

    if ( ( perm_algo == perm_adaptive ) && ( perm_adp_one_perm_per_cu_start > 0 ) )
    {
        std::vector<size_t> perm_adp_significant_snp_indexes(last_block_snps, 0);

        int count = 0;
        for (size_t i = 0; i < num_snps; i++)
        {
            if (perm_adp_results.is_snp_dropped[i] == 0)
            {
                perm_adp_significant_snp_indexes[count++] = i;
            }
        }

        size_t cu_inc_factor = 64;
        perm_adp_drop_interval *= 20;
        num_matrix_blocks = num_cu * cu_inc_factor;

        ocl_tasks.resize(num_matrix_blocks);

        size_t perm_adp_rows_per_block = last_block_snps;
        size_t perm_adp_block_size_bytes = last_block_snps * bytes_per_snp;

        in_mean_elems_per_block = round_to_multiple(last_block_snps, INPUT_MEAN_PAR_ENTRIES);
        in_mean_bytes_per_block = in_mean_elems_per_block * INPUT_MEAN_DATATYPE_SIZE_BYTES;

        out_vec_elems_per_block = round_to_multiple(last_block_snps , OUTPUT_PAR_ENTRIES);
        out_vec_bytes_per_block = out_vec_elems_per_block * OUTPUT_DATATYPE_SIZE_BYTES;

        // std::cout << "Remaining SNPs = " << last_block_snps << std::endl;
        std::cout << "New matrix block size = " << perm_adp_block_size_bytes << " bytes" << std::endl;
        std::cout << "New SNPs per block = " << perm_adp_rows_per_block << std::endl;
        std::cout << "New mean elems per block = " << in_mean_elems_per_block << std::endl;
        std::cout << "New output elems per block = " << out_vec_elems_per_block << std::endl;
        std::cout << "New drop interval = " << perm_adp_drop_interval << std::endl;

        input_matrix_vector_t input_matrix_temp(perm_adp_block_size_bytes, 0);
        input_mean_vector_t input_mean_temp(in_mean_elems_per_block, 0);

        perm_adp_regen_input_data(bytes_per_snp,
                                  snps_per_block,
                                  matrix_block_size_bytes,
                                  perm_adp_rows_per_block,
                                  perm_adp_block_size_bytes,
                                  input_mean_buffer,
                                  input_matrix_buffer,
                                  input_mean_temp,
                                  input_matrix_temp,
                                  perm_adp_results);

        // std::cout << "Starting long term permutation with " << last_block_snps << " SNPs" << std::endl;
        // std::cout << "Starting from the " << perm_adp_one_perm_per_cu_start + 1 << "th permutation" << std::endl;
        // std::cout << "Matrix Size: \t\t" << perm_adp_block_size_bytes << " bytes" << std::endl;
        // std::cout << "Input Mean Size: \t" << in_mean_elems_per_block * INPUT_MEAN_DATATYPE_SIZE_BYTES << " bytes" << std::endl;
        // std::cout << "Output Vector Size: \t" << out_vec_elems_per_block * OUTPUT_DATATYPE_SIZE_BYTES << " bytes" << std::endl;

        std::vector<input_mean_vector_t> input_mean_replicated(num_matrix_blocks, input_mean_vector_t( in_mean_elems_per_block, 0 ) );
        std::vector<input_matrix_vector_t> input_matrix_replicated(num_matrix_blocks, input_matrix_vector_t( perm_adp_block_size_bytes, 0 ) );
        std::vector<input_pheno_vector_t> input_pheno_replicated(num_matrix_blocks, input_pheno_vector_t( pheno_num_elems, 0 ) );
        std::vector<output_vector_t> output_vector_replicated(num_matrix_blocks, output_vector_t( out_vec_elems_per_block, 0 ) );
        std::vector<ocl_task_buffers_t> kernel_bufs(num_matrix_blocks);

        std::cout << "Replicated vectors"
                  << " --- Mean: " << input_mean_replicated.size() << " " << input_mean_replicated[0].size()
                  << " --- Matrix: " << input_matrix_replicated.size() << " " << input_matrix_replicated[0].size()
                  << " --- Pheno: " << input_pheno_replicated.size() << " " << input_pheno_replicated[0].size()
                  << " --- Output: " << output_vector_replicated.size() << " " << output_vector_replicated[0].size() << std::endl;

#pragma omp parallel for
        for (size_t i = 0; i < num_matrix_blocks; i++)
        {
            std::copy(input_mean_temp.begin(),
                      input_mean_temp.end(),
                      input_mean_replicated[i].begin());

            std::copy(input_matrix_temp.begin(),
                      input_matrix_temp.end(),
                      input_matrix_replicated[i].begin());

            std::copy(input_pheno_buffer.begin(),
                      input_pheno_buffer.end(),
                      input_pheno_replicated[i].begin());

            OCL_CHECK(ocl_error_code,
                      kernel_bufs[i].mat_buf =
                          cl::Buffer(api.context,
                                     CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                                     perm_adp_block_size_bytes,
                                     input_matrix_replicated[i].data(),
                                     &ocl_error_code));

            OCL_CHECK(ocl_error_code,
                      kernel_bufs[i].mean_buf =
                          cl::Buffer(api.context,
                                     CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                                     in_mean_bytes_per_block,
                                     input_mean_replicated[i].data(),
                                     &ocl_error_code));

            OCL_CHECK(ocl_error_code,
                      kernel_bufs[i].pheno_buf =
                          cl::Buffer(api.context,
                                     CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                                     pheno_size_bytes,
                                     input_pheno_replicated[i].data(),
                                     &ocl_error_code));

            OCL_CHECK(ocl_error_code,
                      kernel_bufs[i].out_buf =
                          cl::Buffer(api.context,
                                     CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY,
                                     out_vec_bytes_per_block,
                                     output_vector_replicated[i].data(),
                                     &ocl_error_code));
        }

        for (size_t perm = perm_adp_one_perm_per_cu_start + 1; perm < num_perms; perm+=num_matrix_blocks)
        {
            for (unsigned int i = 0; i < num_matrix_blocks; i++)
            {
                std::random_shuffle(input_pheno_replicated[i].begin(), input_pheno_replicated[i].begin() + num_indiv);

                int kernel_num = i % num_cu;

                int narg = 0;
                OCL_CHECK(ocl_error_code, ocl_error_code = kernels[kernel_num].setArg(narg++, (const unsigned int)last_block_snps));
                OCL_CHECK(ocl_error_code, ocl_error_code = kernels[kernel_num].setArg(narg++, (const unsigned int)pheno_num_elems));
                OCL_CHECK(ocl_error_code, ocl_error_code = kernels[kernel_num].setArg(narg++, kernel_bufs[i].mean_buf));
                OCL_CHECK(ocl_error_code, ocl_error_code = kernels[kernel_num].setArg(narg++, kernel_bufs[i].mat_buf));
                OCL_CHECK(ocl_error_code, ocl_error_code = kernels[kernel_num].setArg(narg++, kernel_bufs[i].pheno_buf));
                OCL_CHECK(ocl_error_code, ocl_error_code = kernels[kernel_num].setArg(narg++, kernel_bufs[i].out_buf));

                if (i < num_cu)
                {
                    ocl_tasks[i].run(api, kernels[kernel_num], kernel_bufs[i]);
                }
                else
                {
                    ocl_tasks[i].run(api, kernels[kernel_num], kernel_bufs[i], ocl_tasks[i - num_cu].getDoneEv());
                }
            }

            OCL_CHECK(ocl_error_code, ocl_error_code = api.queue.finish());

            std::thread result_thread(perm_adp_update_results_multi_vec,
                                      output_vector_replicated,
                                      std::cref(vector_stats),
                                      std::cref(gwas_result),
                                      std::ref(perm_adp_results),
                                      std::cref(perm_adp_significant_snp_indexes) );

            if ( ( (perm % perm_adp_drop_interval) < num_matrix_blocks ) )
            {
                result_thread.join();
            }
            else
            {
                result_thread.detach();
            }

            if ( ( (perm % perm_adp_drop_interval) < num_matrix_blocks ) )
            {
                auto prev_snps = last_block_snps;

                size_t num_dropped_snps = perm_adp_drop_snps(perm,
                                                             num_snps,
                                                             perm_adp_min_perms_per_snp,
                                                             perm_adp_results);

                last_block_snps = num_snps - num_dropped_snps;

                size_t min_buf_size = bytes_per_snp * 10;
                if (last_block_snps < prev_snps && perm_adp_block_size_bytes > min_buf_size )
                {
                    perm_adp_block_size_bytes = last_block_snps * bytes_per_snp;

                    in_mean_elems_per_block = round_to_multiple(last_block_snps, INPUT_MEAN_PAR_ENTRIES);
                    in_mean_bytes_per_block = in_mean_elems_per_block * INPUT_MEAN_DATATYPE_SIZE_BYTES;

                    out_vec_elems_per_block = round_to_multiple(last_block_snps , OUTPUT_PAR_ENTRIES);
                    out_vec_bytes_per_block = out_vec_elems_per_block * OUTPUT_DATATYPE_SIZE_BYTES;

#pragma omp parallel for
                    for (size_t i = 0; i < input_matrix_replicated.size(); i++)
                    {
                        input_matrix_replicated[i].resize(perm_adp_block_size_bytes);
                        input_mean_replicated[i].resize(in_mean_elems_per_block);
                        output_vector_replicated[i].resize(out_vec_elems_per_block);

                        OCL_CHECK(ocl_error_code,
                                  kernel_bufs[i].mat_buf =
                                      cl::Buffer(api.context,
                                                 CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                                                 perm_adp_block_size_bytes,
                                                 input_matrix_replicated[i].data(),
                                                 &ocl_error_code));

                        OCL_CHECK(ocl_error_code,
                                  kernel_bufs[i].mean_buf =
                                      cl::Buffer(api.context,
                                                 CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                                                 in_mean_bytes_per_block,
                                                 input_mean_replicated[i].data(),
                                                 &ocl_error_code));

                        OCL_CHECK(ocl_error_code,
                                  kernel_bufs[i].out_buf =
                                      cl::Buffer(api.context,
                                                 CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY,
                                                 out_vec_bytes_per_block,
                                                 output_vector_replicated[i].data(),
                                                 &ocl_error_code));
                    }

                    std::cout << "New buffer size = " << perm_adp_block_size_bytes
                        << " New SNPs per buffer = " << ( perm_adp_block_size_bytes / bytes_per_snp ) << std::endl;
                }

                perm_adp_regen_and_replicate_input_data(bytes_per_snp,
                                                        snps_per_block,
                                                        matrix_block_size_bytes,
                                                        in_mean_elems_per_block,
                                                        perm_adp_block_size_bytes,
                                                        input_mean_buffer,
                                                        input_matrix_buffer,
                                                        perm_adp_results,
                                                        perm_adp_significant_snp_indexes,
                                                        input_mean_replicated,
                                                        input_matrix_replicated);

                std::cout << " Perm = " << perm
                          << " Dropped SNPs = " << num_dropped_snps
                          << " Remaining SNPs = " << (num_snps - num_dropped_snps) << std::endl;

                if (perm_adp_drop_counter == 10)
                {
                    perm_adp_drop_interval *= 2;

                    if (perm_adp_drop_interval >= 5000000)
                    {
                        perm_adp_drop_interval = 5000000;
                    }

                    std::cout << "Drop interval = " << perm_adp_drop_interval << std::endl;
                    perm_adp_drop_counter = 0;

                    // Write intermediate results to file
                    if (perm > 50000000)
                    {
                        std::ostringstream iss;
                        iss << "perm_adaptive." << perm << ".txt";
                        std::string filename = iss.str();

                        std::ofstream res_file(filename);
                        std::ifstream bim_file_stream(bim_file);

                        time = std::chrono::system_clock::now();
                        time_out = std::chrono::system_clock::to_time_t(time);

                        res_file << "Write time: " << std::ctime(&time_out) << std::endl;
                        size_t idx = 0;
                        for (size_t i = 0; i < num_snps; i++)
                        {
                            std::string line;
                            std::getline(bim_file_stream, line);

                            if (i != perm_adp_significant_snp_indexes[idx])
                                continue;

                            size_t num_perms;
                            double p_val;
                            if (perm_adp_results.perm_dropped[i] == 0)
                            {
                                num_perms = perm;
                            }
                            else
                            {
                                num_perms = perm_adp_results.perm_dropped[i];
                            }

                            p_val = (double)(perm_adp_results.perm_count[i] + 1) / (num_perms + 1);

                            std::string buf;
                            std::string snp_name;
                            std::stringstream stream(line);
                            stream >> buf;
                            stream >> snp_name;
                            res_file << std::left << std::setw(30) << snp_name.c_str()
                                     << std::left << std::setw(25) << gwas_result[i]
                                     << std::left << std::setw(25) << p_val
                                     << std::left << std::setw(25) << num_perms << std::endl;

                            idx++;
                        }
                    }
                }

                perm_adp_drop_counter++;
            }
        }

        size_t num_dropped_snps = perm_adp_drop_snps(num_perms,
                                                     num_snps,
                                                     perm_adp_min_perms_per_snp,
                                                     perm_adp_results);

        std::cout << num_perms << " adaptive permutations complete." << std::endl
                  << "Dropped SNPs = " << num_dropped_snps << std::endl
                  << "Remaining SNPs = " << (num_snps - num_dropped_snps) << std::endl;
    }

    auto end_processing_time = std::chrono::high_resolution_clock::now();

    auto start_results_export_time = std::chrono::high_resolution_clock::now();

    if ( perm_algo == perm_adaptive )
    {
        std::cout << "Writing adaptive permutation results to perm_adaptive.res.txt" << std::endl;

        std::vector<double> perm_result(num_snps);
        std::vector<size_t> num_perms_vec(num_snps);

        std::ofstream res_file("perm_adaptive.res.txt");

#pragma omp parallel for
        for (size_t i = 0; i < num_snps; i++)
        {
            double p_val;

            if (perm_adp_results.perm_dropped[i] == 0)
            {
                p_val = (double)(perm_adp_results.perm_count[i] + 1) / (num_perms + 1);
                num_perms_vec[i] = num_perms;
            }
            else
            {
                p_val = (double)(perm_adp_results.perm_count[i] + 1) / (perm_adp_results.perm_dropped[i] + 1);
                num_perms_vec[i] = perm_adp_results.perm_dropped[i];
            }

            perm_result[i] = p_val;
        }

        std::ifstream bim_file_stream(bim_file);
        // size_t idx_count = 0;
        for (size_t i = 0; i < num_snps; i++)
        {
            std::string line;
            std::getline(bim_file_stream, line);

            // if (i != perm_adp_significant_snp_indexes[idx_count]) continue;

            std::string buf;
            std::string snp_name;
            std::stringstream stream(line);
            stream >> buf;
            stream >> snp_name;
            res_file << std::left << std::setw(30) << snp_name.c_str()
                << std::left << std::setw(25) << gwas_result[i]
                << std::left << std::setw(25) << perm_result[i]
                << std::left << std::setw(25) << num_perms_vec[i] << std::endl;

            // idx_count++;
        }
    }
    else if ( perm_algo == perm_maxT )
    {
        // std::ofstream significant_95("perm_maxT.95th.txt");
        // std::ofstream significant_99("perm_maxT.99th.txt");
        // int percentile_95_pos = ((double)perm_maxF_res_vector.size() * 0.05);
        // int percentile_99_pos = ((double)perm_maxF_res_vector.size() * 0.01);

        // double percentile_95 = perm_maxF_res_vector[percentile_95_pos];
        // double percentile_99 = perm_maxF_res_vector[percentile_99_pos];

        // std::ofstream best_file("perm_maxT.best.txt");
        // for (size_t i = 0; i < perm_maxF_res_vector.size(); i++)
        // {
        //     best_file << std::left << std::setw(10) << i << "Max F = " << perm_maxF_res_vector[i] << std::endl;
        // }

        std::cout << "Writing maxT permutation results to perm_maxT.res.txt" << std::endl;

        std::ofstream res_file("perm_maxT.res.txt");

        std::vector<double> perm_result(num_snps);

        perm_maxF_res_vector.erase(perm_maxF_res_vector.begin());
        std::sort(perm_maxF_res_vector.begin(), perm_maxF_res_vector.end(), std::greater<double>());

#pragma omp parallel for
        for (size_t i = 0; i < num_snps; i++)
        {
            double test_val = std::abs(gwas_result[i]);
            int rank = num_perms + 1;

            for (size_t j = 0; j < perm_maxF_res_vector.size(); j++)
            {
                if (test_val >= perm_maxF_res_vector[j])
                {
                    rank = j + 1;
                    break;
                }
            }

            double p_val = (double)rank / (num_perms + 1);
            perm_result[i] = p_val;
        }

        // significant_95 << "95th percentile = " << percentile_95 << std::endl << std::endl;
        // significant_99 << "99th percentile = " << percentile_99 << std::endl << std::endl;

        std::ifstream bim_file_stream(bim_file);
        for (size_t i = 0; i < num_snps; i++)
        {
            std::string line;
            std::getline(bim_file_stream, line);

            std::string buf;
            std::string snp_name;
            std::stringstream stream(line);
            stream >> buf;
            stream >> snp_name;
            res_file << std::left << std::setw(30) << snp_name.c_str()
                << std::left << std::setw(25) << gwas_result[i]
                << std::left << std::setw(25) << perm_result[i] << std::endl;

            // if ( std::abs(gwas_result[i]) > percentile_95 )
            // {
            //     significant_95 << std::left << std::setw(30) << snp_name.c_str()
            //         << std::left << std::setw(25) << gwas_result[i]
            //         << std::left << std::setw(25) << perm_result[i] << std::endl;
            // }

            // if ( std::abs(gwas_result[i]) > percentile_99 )
            // {
            //     significant_99 << std::left << std::setw(30) << snp_name.c_str()
            //         << std::left << std::setw(25) << gwas_result[i]
            //         << std::left << std::setw(25) << perm_result[i] << std::endl;
            // }
        }
    }
    else if ( perm_algo == perm_regression_only )
    {
        std::cout << "Writing GWAS results to gwas_results.txt" << std::endl;

        std::ofstream res_file("gwas_results.txt");

        std::ifstream bim_file_stream(bim_file);
        for (size_t i = 0; i < num_snps; i++)
        {
            std::string line;
            std::getline(bim_file_stream, line);

            std::string buf;
            std::string snp_name;
            std::stringstream stream(line);
            stream >> buf;
            stream >> snp_name;

            unsigned int idx_start = ( i / snps_per_block ) * out_vec_elems_per_block;
            unsigned int output_idx = idx_start + (i % snps_per_block);

            double dot_prod = fpga_output_buffer[output_idx].to_double();
            double beta = ( dot_prod / vector_stats.X[i].sq_sum ) / phenotype_scale_factor;
            double denom = vector_stats.X[i].std_dev * vector_stats.Y.std_dev;
            double corr_sq = ( dot_prod * dot_prod ) / ( denom * denom );
            double F = std::sqrt((corr_sq / (1 - corr_sq)) * (vector_stats.X[i].count - 2));

            res_file << std::left << std::setw(20) << snp_name.c_str()
                     << std::left << std::setw(8) << vector_stats.X[i].count
                     // << std::left << std::setw(10) << vector_stats.X[i].sq_sum
                     // << std::left << std::setw(10) << vector_stats.X[i].std_dev
                     // << std::left << std::setw(10) << vector_stats.Y.std_dev
                     // << std::left << "Dot prod = " << std::setw(15) << dot_prod
                     // << std::left << "Denom = " << std::setw(10) << denom
                     << std::left << std::setw(15) << beta
                     << std::left << F << std::endl;
        }
    }

    auto end_results_export_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> import_time = end_import_time - start_import_time;
    std::chrono::duration<double> processing_time = end_processing_time - start_processing_time;
    std::chrono::duration<double> results_export_time = end_results_export_time - start_results_export_time;

    std::cout << std::endl << "Data import time: " << import_time.count() << " s" << std::endl;
    std::cout << "FPGA time: " << processing_time.count() << " s" << std::endl;
    std::cout << "Results export time: " << results_export_time.count() << " s" << std::endl;
    std::cout << "Total runtime: " << ( processing_time.count() + import_time.count() + results_export_time.count() ) << " s" << std::endl;

    time = std::chrono::system_clock::now();
    time_out = std::chrono::system_clock::to_time_t(time);

    std::cout << "End time: " << std::ctime(&time_out) << std::endl;

    return 0;
}
