#ifndef DATA_SETUP_H
#define DATA_SETUP_H

#include <vector>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <random>
#include <bitset>
#include <fstream>
#include <iomanip>

#include "ap_int.h"
#include "ap_fixed.h"

#include "Types.h"

template <typename T>
void read_file_to_vector(char *file_name, std::vector<T> &v)
{
    std::ifstream file_stream(file_name, std::ios::in | std::ios::binary);
    file_stream.seekg(0, std::ios::end);
    size_t filesize = file_stream.tellg();
    file_stream.seekg(0, std::ios::beg);

    file_stream.read((char *)v.data(), filesize);
}

inline int round_to_multiple(int input, int multiple)
{
    return ((input + multiple - 1) / multiple) * multiple;
}

inline void set_bits(std::bitset<8> &in_byte, int pos, int num)
{
    switch (num)
    {
    case 3:
        in_byte[pos] = 1;
        in_byte[pos + 1] = 1;
        break;

    case 2:
        in_byte[pos] = 0;
        in_byte[pos + 1] = 1;
        break;

    case 1:
        in_byte[pos] = 1;
        in_byte[pos + 1] = 0;
        break;

    case 0:
        in_byte[pos] = 0;
        in_byte[pos + 1] = 0;
        break;

    default:
        break;
    }
}

inline int convert_bed_val(int bed_val)
{
    switch (bed_val)
    {
    case 3: //Homozygous for second (major) allele in .bim file
        return 0;
        break;

    case 2: //Heterozygous
        return 1;
        break;

    case 1: //Missing genotype i.e. skip in all calcs
        return 3;
        break;

    case 0: //Homozygous for first (minor) allele in .bim file
        return 2;
        break;

    default:
        return 0;
        break;
    }
}

inline void convert_bed_byte_to_fpga_byte(
    char bed_byte,
    std::bitset<8> &fpga_byte,
    int &byte_sum,
    int &byte_count,
    std::vector<char> &converted_snp_data,
    size_t compressed_col_idx,
    size_t num_snp_vals = 4)
{
    // int num_snp_vals = (compressed_col_idx == num_compressed_cols - 1) ? (converted_snp_data.size() % 4) : 4;

    for (size_t i = 0; i < num_snp_vals; i++)
    {
        int bed_val = (bed_byte >> (i * 2)) & 0b00000011;
        int snp_val = convert_bed_val(bed_val);

        converted_snp_data[(compressed_col_idx * 4) + i] = snp_val;
        set_bits(fpga_byte, (i * 2), snp_val);

        if (snp_val != 3)
        {
            byte_sum += snp_val;
            byte_count += 1;
        }
    }
}

struct square
{
    square(double mean) : mean(mean) {}

    inline double operator()(const double &Left, const char &Right) const
    {
        if (Right != 3)
        {
            double temp = Right - mean;
            return (Left + temp * temp);
        }
        else
        {
            return Left;
        }
    }

private:
    double mean;
};

vector_stats_t generate_fpga_data(
    unsigned int num_rows,
    unsigned int num_cols,
    unsigned long int block_size_bytes,
    std::vector<char> bed_file_buffer,
    std::vector<double> phenotype_buffer,
    double phenotype_scale_factor,
    input_mean_vector_t &input_mean_buffer,
    input_matrix_vector_t &input_matrix_buffer,
    input_pheno_vector_t &input_pheno_buffer)
{
    vector_stats_t vector_stats(num_rows);

    vector_stats.Y.count = num_cols;
    vector_stats.Y.sq_sum = 0;

    for (size_t i = 0; i < num_cols; i++)
    {
        double phenotype = phenotype_buffer[i] * phenotype_scale_factor;
        input_pheno_buffer[i] = phenotype;

        vector_stats.Y.sq_sum += phenotype * phenotype;
    }

    vector_stats.Y.std_dev = std::sqrt(vector_stats.Y.sq_sum);

    unsigned int input_bytes_per_snp = ceil(num_cols / 4.0);
    unsigned int output_bytes_per_snp = round_to_multiple(input_bytes_per_snp, (INPUT_MAT_PAR_ELEMS / 4));

    unsigned int snps_per_block = block_size_bytes / output_bytes_per_snp;
    unsigned int row_incr = round_to_multiple(snps_per_block, (4096 / INPUT_MEAN_DATATYPE_SIZE_BYTES));

    unsigned int last_byte_snps = (num_cols % 4 == 0) ? 4 : (num_cols % 4);

#pragma omp parallel for
    for (size_t row = 0; row < num_rows; row++)
    {
        auto snp_start = bed_file_buffer.begin() + (row * input_bytes_per_snp);
        auto snp_end = bed_file_buffer.begin() + ((row + 1) * input_bytes_per_snp);

        std::vector<char> bed_file_snp(snp_start, snp_end);
        std::vector<char> converted_snp(num_cols);

        unsigned int row_sum = 0;
        unsigned int row_count = 0;

        // #pragma omp parallel for reduction(+: row_sum, row_count)
        for (size_t byte_idx = 0; byte_idx < input_bytes_per_snp; byte_idx++)
        {
            unsigned char bed_file_byte = bed_file_snp[byte_idx];
            std::bitset<8> converted_byte;

            int byte_sum = 0;
            int byte_count = 0;

            if (byte_idx == input_bytes_per_snp - 1)
            {
                convert_bed_byte_to_fpga_byte(bed_file_byte, converted_byte, byte_sum, byte_count, converted_snp, byte_idx, last_byte_snps);
            }
            else
            {
                convert_bed_byte_to_fpga_byte(bed_file_byte, converted_byte, byte_sum, byte_count, converted_snp, byte_idx);
            }

            row_sum += byte_sum;
            row_count += byte_count;

            size_t block_start_index = (row / snps_per_block) * block_size_bytes;
            size_t block_offset = ((row % snps_per_block) * output_bytes_per_snp);
            size_t matrix_index = block_start_index + block_offset + byte_idx;
            unsigned char fpga_byte = converted_byte.to_ulong() & 0xFF;
            input_matrix_buffer[matrix_index] = fpga_byte;
        }

        mean_stddev_t mean_stddev;
        mean_stddev.count = row_count;
        mean_stddev.mean = (double)row_sum / row_count;
        mean_stddev.sq_sum = std::accumulate(converted_snp.begin(), converted_snp.end(), 0.0, square(mean_stddev.mean));

        mean_stddev.std_dev = std::sqrt(mean_stddev.sq_sum);

        vector_stats.X[row] = mean_stddev;

        unsigned int mean_start = (row / snps_per_block) * row_incr;
        unsigned int mean_index = mean_start + (row % snps_per_block);
        if (std::isnan(vector_stats.X[row].mean))
        {
            input_mean_buffer[mean_index] = 0;
        }
        else
        {
            input_mean_buffer[mean_index] = vector_stats.X[row].mean;
        }
    }

    if (0)
    {
        std::cout << std::endl
                  << "Writing " << 1 << " * " << input_matrix_buffer.size() << " bytes to input_matrix_buffer.bin" << std::endl;
        std::ofstream fout_mat("input_matrix_buffer.bin", std::ios::out | std::ios::binary);
        fout_mat.write((char *)&input_matrix_buffer[0], input_matrix_buffer.size());
        fout_mat.close();

        std::cout << std::endl
                  << "Writing " << sizeof(input_pheno_data_type_t) << " * " << input_pheno_buffer.size() << " bytes to input_pheno_buffer.bin" << std::endl;
        std::ofstream fout_vec("input_pheno_buffer.bin", std::ios::out | std::ios::binary);
        fout_vec.write((char *)&input_pheno_buffer[0], input_pheno_buffer.size() * sizeof(input_pheno_data_type_t));
        fout_vec.close();

        std::cout << std::endl
                  << "Writing " << sizeof(input_mean_data_type_t) << " * " << input_mean_buffer.size() << " bytes to input_mean_buffer.bin" << std::endl;
        std::ofstream fout_mean("input_mean_buffer.bin", std::ios::out | std::ios::binary);
        fout_mean.write((char *)&input_mean_buffer[0], input_mean_buffer.size() * sizeof(input_mean_data_type_t));
        fout_mean.close();
    }

    return vector_stats;
}

#endif