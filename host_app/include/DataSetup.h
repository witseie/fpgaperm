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

template<typename T>
void read_file_to_vector(char * file_name, std::vector<T> &v)
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
    case 3:  //Homozygous for second (major) allele in .bim file
        return 0;
        break;

    case 2:  //Heterozygous
        return 1;
        break;

    case 1:  //Missing genotype i.e. skip in all calcs
        return 3;
        break;

    case 0:  //Homozygous for first (minor) allele in .bim file
        return 2;
        break;

    default:
        return 0;
        break;
    }
}

void convert_bed_byte_to_fpga_byte(
    char bed_byte,
    std::bitset<8>& fpga_byte,
    int& byte_sum,
    int& byte_count,
    std::vector<char> &converted_snp_data,
    size_t compressed_col_idx,
    size_t num_compressed_cols
    )
{
    int num_snp_vals = (compressed_col_idx == num_compressed_cols - 1) ? (converted_snp_data.size() % 4) : 4;

    for (size_t i = 0; i < num_snp_vals; i++)
    {
        int bed_val = ( bed_byte >> ( i * 2 ) ) & 0b00000011;
        int snp_val = convert_bed_val(bed_val);

        converted_snp_data[(compressed_col_idx * 4) + i] = snp_val;
        set_bits(fpga_byte, ( i * 2 ), snp_val);

        if ( snp_val != 3 )
        {
            byte_sum += snp_val;
            byte_count += 1;      
        }
    }
}

struct square
{
    square(double mean) : mean(mean) {}

    double operator()(const double &Left, const char &Right) const
    {
        if (Right < 3)
        {
            double temp = Right - mean;
            return ( Left + temp * temp );
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
    input_mean_vector_t &input_mean,
    input_matrix_vector_t &input_matrix,
    input_pheno_vector_t &input_pheno)
{
    vector_stats_t vector_stats(num_rows);

    vector_stats.Y.count = num_cols;
    vector_stats.Y.sq_sum = 0;
    for (size_t i = 0; i < num_cols; i++)
    {
        input_pheno[i] = phenotype_buffer[i];

        vector_stats.Y.sq_sum += phenotype_buffer[i] * phenotype_buffer[i];
    }

    vector_stats.Y.std_dev = std::sqrt(vector_stats.Y.sq_sum);

    unsigned int bytes_per_row = ceil(num_cols / 4.0);
    unsigned int cols_per_block = round_to_multiple(bytes_per_row, INPUT_MAT_PAR_ENTRIES);

    unsigned int rows_per_block = block_size_bytes / cols_per_block;
    unsigned int row_incr = round_to_multiple(rows_per_block, ( 4096 / INPUT_MEAN_DATATYPE_SIZE_BYTES ) );

#pragma omp parallel for
    for (size_t row = 0; row < num_rows; row++)
    {
        auto start_row = bed_file_buffer.begin() + ( row * bytes_per_row );
        auto end_row = bed_file_buffer.begin() + ( ( row + 1 )  * bytes_per_row );

        std::vector<char> bed_row(start_row, end_row);
        std::vector<char> converted_row(num_cols);

        unsigned int row_sum = 0;
        unsigned int row_count = 0;

#pragma omp parallel for reduction(+: row_sum, row_count)
        for (size_t byte_idx = 0; byte_idx < bytes_per_row; byte_idx++)
        {
            unsigned char bed_file_byte = bed_row[byte_idx];
            std::bitset<8> converted_byte;

            int byte_sum = 0;
            int byte_count = 0;

            convert_bed_byte_to_fpga_byte(bed_file_byte, converted_byte, byte_sum, byte_count, converted_row, byte_idx, bytes_per_row);

            row_sum += byte_sum;
            row_count += byte_count;

            size_t block_start_index = ( row / rows_per_block ) * block_size_bytes;
            size_t block_offset = ( (row % rows_per_block) * cols_per_block );
            size_t matrix_index = block_start_index + block_offset + byte_idx;
            unsigned char fpga_byte = converted_byte.to_ulong() & 0xFF;
            input_matrix[matrix_index] = fpga_byte;
        }

        mean_stddev_t mean_stddev;
        mean_stddev.count = row_count;
        mean_stddev.mean = (double)row_sum / row_count;
        mean_stddev.sq_sum = std::accumulate( converted_row.begin(), converted_row.end(), 0.0, square(mean_stddev.mean) );

        mean_stddev.std_dev = std::sqrt(mean_stddev.sq_sum);

        vector_stats.X[row] = mean_stddev;

        unsigned int mean_start = ( row / rows_per_block ) * row_incr;
        unsigned int mean_index = mean_start + (row % rows_per_block);
        if (std::isnan(vector_stats.X[row].mean))
        {
            input_mean[mean_index] = 0;
        }
        else
        {
            input_mean[mean_index] = vector_stats.X[row].mean;
        }
    }

    if (0)
    {
        std::cout << std::endl
                  << "Writing " << 1 << " * " << input_matrix.size() << " bytes to input_matrix.bin" << std::endl;
        std::ofstream fout_mat("input_matrix_test.bin", std::ios::out | std::ios::binary);
        fout_mat.write((char *)&input_matrix[0], input_matrix.size());
        fout_mat.close();

        std::cout << std::endl
                  << "Writing " << sizeof(input_pheno_data_type_t) << " * " << input_pheno.size() << " bytes to input_pheno.bin" << std::endl;
        std::ofstream fout_vec("input_pheno.bin", std::ios::out | std::ios::binary);
        fout_vec.write((char *)&input_pheno[0], input_pheno.size() * sizeof(input_pheno_data_type_t));
        fout_vec.close();

        std::cout << std::endl
                  << "Writing " << sizeof(input_mean_data_type_t) << " * " << input_mean.size() << " bytes to input_mean.bin" << std::endl;
        std::ofstream fout_mean("input_mean.bin", std::ios::out | std::ios::binary);
        fout_mean.write((char *)&input_mean[0], input_mean.size() * sizeof(input_mean_data_type_t));
        fout_mean.close();
    }

    return vector_stats;
}

#endif