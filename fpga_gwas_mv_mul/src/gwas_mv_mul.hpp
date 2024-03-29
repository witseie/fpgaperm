#ifndef GWAS_MV_MUL_H
#define GWAS_MV_MUL_H

#include "fpga_types.hpp"

void mv_mul(const unsigned int num_rows,
            const unsigned int num_cols,
            hls::stream<mean_wide_type_t> &mean_stream,
            hls::stream<matrix_wide_type_t> &matrix_stream,
            ap_int<DATA_WIDTH> *vector_buffer,
            hls::stream<output_wide_type_t> &output_stream)
{
    const unsigned int blocks_per_row = num_cols / INPUT_MAT_PAR_ELEMS;

    output_wide_type_t output_wt;
#pragma HLS ARRAY_PARTITION variable = output_wt complete dim = 1

mv_mul:
    for (unsigned int row = 0; row < num_rows; ++row)
    {
        acc_data_type_t dot_prod = 0;

    dot_prod:
        for (unsigned int i = 0; i < blocks_per_row; ++i)
        {
#pragma HLS PIPELINE

            mean_wide_type_t mean_wt;
#pragma HLS ARRAY_PARTITION variable = mean_wt complete dim = 1

            if ((i == 0) && (row % INPUT_MEAN_PAR_ELEMS == 0))
            {
                mean_wt = mean_stream.read();
            }

            matrix_wide_type_t mat_wt = matrix_stream.read();
#pragma HLS ARRAY_PARTITION variable = mat_wt complete dim = 1

        dot_prod_block:
            for (unsigned int j = 0; j < INPUT_MAT_PAR_ELEMS; ++j)
            {
#pragma HLS UNROLL
                unsigned int vec_wt_index = ((i * INPUT_MAT_PAR_ELEMS) + j) / INPUT_VEC_PAR_ELEMS;

                vector_wide_type_t vec_wt = vector_buffer[vec_wt_index];
#pragma HLS ARRAY_PARTITION variable = vec_wt complete dim = 1

                input_matrix_data_type_t temp_mat = mat_wt[j];

                mul_data_type_t temp_product;
                if (temp_mat == 3)
                {
                    temp_product = 0;
                }
                else
                {
                    input_mean_data_type_t temp_mean = mean_wt[row % INPUT_MEAN_PAR_ELEMS];
                    int_vec_data_type_t temp_vec = vec_wt[j % INPUT_VEC_PAR_ELEMS];
                    int_snp_data_type_t temp_mat_sub = temp_mat - temp_mean;

                    temp_product = temp_mat_sub * temp_vec;
                }

                dot_prod += temp_product;
            }

            if (i == (blocks_per_row - 1))
            {
                unsigned int wt_out_index = row % OUTPUT_PAR_ELEMS;

                output_wt[wt_out_index] = dot_prod;

                if (((wt_out_index + 1) == OUTPUT_PAR_ELEMS) || (row == num_rows - 1))
                {
                    output_stream.write(output_wt);
                }
            }
        }
    }
}

#endif