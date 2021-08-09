#ifndef PERMUTATION_H
#define PERMUTATION_H

#include <vector>
#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "Types.h"
#include "DataSetup.h"

void perm_adp_regen_and_replicate_input_data(
    size_t bytes_per_row,
    size_t rows_per_block,
    size_t block_size_bytes,
    size_t num_mean_elems,
    size_t num_matrix_elems,
    const input_mean_vector_t &input_mean,
    const input_matrix_vector_t &input_matrix,
    const perm_adp_results_t &perm_adp_results,
    const std::vector<size_t> &perm_adp_significant_snp_indexes,
    std::vector<input_mean_vector_t> &input_mean_replicated,
    std::vector<input_matrix_vector_t> &input_matrix_replicated)
{
    input_mean_vector_t mean_regen(num_mean_elems, 0);
    input_matrix_vector_t matrix_regen(num_matrix_elems, 0);

    size_t input_mean_bytes_per_block = round_to_multiple(rows_per_block, (4096 / INPUT_MEAN_DATATYPE_SIZE_BYTES));

    size_t output_idx = 0;
    for (size_t i = 0; i < perm_adp_significant_snp_indexes.size(); i++)
    {
        size_t snp_index = perm_adp_significant_snp_indexes[i];

        if (perm_adp_results.is_snp_dropped[snp_index])
        {
            continue;
        }

        // size_t input_mat_idx = i * bytes_per_row;
        size_t input_block = (snp_index / rows_per_block);

        size_t input_mat_start = input_block * block_size_bytes;
        size_t input_mat_idx = input_mat_start + ((snp_index % rows_per_block) * bytes_per_row);

        size_t output_mat_idx = output_idx * bytes_per_row;

        std::copy(input_matrix.data() + input_mat_idx,
                  input_matrix.data() + input_mat_idx + bytes_per_row,
                  matrix_regen.data() + output_mat_idx);

        size_t input_mean_start = input_block * input_mean_bytes_per_block;
        size_t input_mean_idx = input_mean_start + (snp_index % rows_per_block);

        mean_regen[output_idx] = input_mean[input_mean_idx];

        output_idx++;
    }

    int num_buffers = input_mean_replicated.size();

    for (size_t i = 0; i < num_buffers; i++)
    {
        std::copy(mean_regen.begin(),
                  mean_regen.end(),
                  input_mean_replicated[i].begin());

        std::copy(matrix_regen.begin(),
                  matrix_regen.end(),
                  input_matrix_replicated[i].begin());
    }
}

void perm_adp_regen_input_data(
    size_t bytes_per_row,
    size_t rows_per_block,
    size_t block_size_bytes,
    size_t output_rows_per_block,
    size_t output_block_size_bytes,
    const input_mean_vector_t &input_mean,
    const input_matrix_vector_t &input_matrix,
    input_mean_vector_t &input_mean_regen,
    input_matrix_vector_t &input_matrix_regen,
    const perm_adp_results_t &perm_adp_results)
{
    size_t input_mean_bytes_per_block = round_to_multiple(rows_per_block, (4096 / INPUT_MEAN_DATATYPE_SIZE_BYTES));
    size_t output_mean_bytes_per_block = round_to_multiple(output_rows_per_block, (4096 / INPUT_MEAN_DATATYPE_SIZE_BYTES));

    size_t output_idx = 0;
    for (size_t i = 0; i < perm_adp_results.is_snp_dropped.size(); i++)
    {
        if (perm_adp_results.is_snp_dropped[i])
        {
            continue;
        }

        size_t input_block = (i / rows_per_block);
        size_t input_mat_start = input_block * block_size_bytes;
        size_t input_mat_idx = input_mat_start + ((i % rows_per_block) * bytes_per_row);

        size_t output_block = (output_idx / output_rows_per_block);
        size_t output_mat_start = output_block * output_block_size_bytes;
        size_t output_mat_idx = output_mat_start + ((output_idx % output_rows_per_block) * bytes_per_row);

        std::copy(input_matrix.data() + input_mat_idx,
                  input_matrix.data() + input_mat_idx + bytes_per_row,
                  input_matrix_regen.data() + output_mat_idx);

        size_t input_mean_start = input_block * input_mean_bytes_per_block;
        size_t input_mean_idx = input_mean_start + (i % rows_per_block);

        size_t output_mean_start = output_block * output_mean_bytes_per_block;
        size_t output_mean_idx = output_mean_start + (output_idx % output_rows_per_block);

        input_mean_regen[output_mean_idx] = input_mean[input_mean_idx];

        output_idx++;
    }
}

int perm_adp_drop_snps(
    size_t perm,
    size_t num_rows,
    size_t min_perms_per_snp,
    perm_adp_results_t &perm_adp_results)
{
    size_t num_dropped_snps = 0;

    for (size_t i = 0; i < num_rows; i++)
    {
        if (perm_adp_results.is_snp_dropped[i] == 1)
        {
            num_dropped_snps++;
        }
        else
        {
            if (perm_adp_results.perm_count[i] >= min_perms_per_snp)
            {
                perm_adp_results.perm_dropped[i] = perm;
                perm_adp_results.is_snp_dropped[i] = 1;
                num_dropped_snps++;
            }
        }
    }

    return num_dropped_snps;
}

void perm_adp_update_results_multi_vec(
    const std::vector<output_vector_t> output_vector_fpga,
    const vector_stats_t &vector_stats,
    const std::vector<double> &gwas_result,
    perm_adp_results_t &perm_adp_results,
    const std::vector<size_t> &perm_adp_significant_snp_indexes)
{
    for (size_t j = 0; j < output_vector_fpga.size(); j++)
    {
        auto vec = output_vector_fpga[j];
        size_t output_idx = 0;

        for (size_t i = 0; i < perm_adp_significant_snp_indexes.size(); i++)
        {
            size_t snp_index = perm_adp_significant_snp_indexes[i];

            if (perm_adp_results.is_snp_dropped[snp_index])
                continue;

            double dot_prod = vec[output_idx].to_double();
            double denom = vector_stats.X[snp_index].std_dev * vector_stats.Y.std_dev;
            double corr_sq = (dot_prod * dot_prod) / (denom * denom);

            double F = std::sqrt((corr_sq / (1 - corr_sq)) * (vector_stats.X[snp_index].count - 2));

            output_idx++;

            if (F > std::abs(gwas_result[snp_index]))
            {
                perm_adp_results.perm_count[snp_index]++;
            }
        }
    }
}

void perm_adp_update_results(
    size_t perm,
    size_t num_rows,
    size_t rows_per_block,
    size_t elems_per_block,
    const output_vector_t output_vector_fpga,
    const vector_stats_t &vector_stats,
    std::vector<double> &gwas_result,
    perm_adp_results_t &perm_adp_results)
{
    size_t out_vec_index = 0;

    for (size_t i = 0; i < num_rows; i++)
    {
        if (perm_adp_results.is_snp_dropped[i])
        {
            continue;
        }

        size_t idx_start = (out_vec_index / rows_per_block) * elems_per_block;
        size_t output_idx = idx_start + (out_vec_index % rows_per_block);
        out_vec_index++;

        double dot_prod = output_vector_fpga[output_idx].to_double();
        // double beta = ( dot_prod / vector_stats.X[i].sq_sum );
        double denom = vector_stats.X[i].std_dev * vector_stats.Y.std_dev;
        double corr_sq = (dot_prod * dot_prod) / (denom * denom);

        double F = std::sqrt((corr_sq / (1 - corr_sq)) * (vector_stats.X[i].count - 2));

        if (perm == 0)
        {
            if (std::isnan(F))
            {
                gwas_result[i] = F;
                perm_adp_results.is_snp_dropped[i] = 1;
            }
            else
            {
                double neg_sign = (std::signbit(dot_prod)) ? -1.0 : 1.0;
                gwas_result[i] = neg_sign * F;
            }
        }
        else
        {
            if (F > std::abs(gwas_result[i]))
            {
                perm_adp_results.perm_count[i]++;
            }
        }
    }
}

void perm_maxT_update_results(
    size_t perm,
    size_t num_rows,
    size_t rows_per_block,
    size_t elems_per_block,
    const output_vector_t output_vector_fpga,
    const vector_stats_t &vector_stats,
    std::vector<double> &gwas_result,
    std::vector<double> &perm_maxF_res_vector)
{
    double max_F = 0.0;

    for (size_t i = 0; i < num_rows; i++)
    {
        unsigned int idx_start = (i / rows_per_block) * elems_per_block;
        unsigned int output_idx = idx_start + (i % rows_per_block);

        double dot_prod = output_vector_fpga[output_idx].to_double();
        // double beta = ( dot_prod / vector_stats.X[i].sq_sum );
        double denom = vector_stats.X[i].std_dev * vector_stats.Y.std_dev;
        double corr_sq = (dot_prod * dot_prod) / (denom * denom);

        double F = std::sqrt((corr_sq / (1 - corr_sq)) * (vector_stats.X[i].count - 2));

        if (F > max_F)
        {
            max_F = F;
        }

        if (perm == 0)
        {
            double neg_sign = (std::signbit(dot_prod)) ? -1.0 : 1.0;
            gwas_result[i] = neg_sign * F;
        }

        // std::cout << i << std::setw (8)
        // 	<< vector_stats.X[i].count << std::setw (15)
        // 	<< vector_stats.X[i].sq_sum << std::setw (15)
        // 	<< vector_stats.X[i].std_dev << std::setw (15)
        // 	<< vector_stats.Y.std_dev << std::setw (15)
        //  << "Dot prod = " << dot_prod << std::setw (15)
        //  << "Denom = " << denom << std::setw (15)
        //  << "Beta = " << beta << std::setw (15)
        //  << "F = " << gwas_result[i] << std::endl;
    }

    perm_maxF_res_vector[perm] = max_F;

    // std::cout << perm << " " << max_F << " " << perm_maxF_res_vector[perm] << std::endl;
}

#endif