#ifndef DATA_TRANSFER_H
#define DATA_TRANSFER_H

#include "hls_stream.h"

#include "fpga_types.hpp"

void readMeanToStream(
    const unsigned int num_rows,
    const ap_int<DATA_WIDTH> *in_ptr,
    hls::stream<mean_wide_type_t> &out_stream)
{
    double rows = num_rows;
    const unsigned int num_blocks = hls::ceil(rows / INPUT_MEAN_PAR_ENTRIES);

read_mean:
    for (unsigned int i = 0; i < num_blocks; ++i)
    {
#pragma HLS PIPELINE
        mean_wide_type_t l_val = in_ptr[i];
        out_stream.write(l_val);
    }
}

void readMatToStream(
    const unsigned int num_rows,
    const unsigned int num_cols,
    const ap_int<DATA_WIDTH> *in_ptr,
    hls::stream<matrix_wide_type_t> &out_stream)
{
    const unsigned int num_blocks = num_rows * num_cols / INPUT_MAT_PAR_ENTRIES;

read_mat:
    for (unsigned int i = 0; i < num_blocks; ++i)
    {
#pragma HLS PIPELINE
        ap_int<DATA_WIDTH> in_wt = in_ptr[i];
        matrix_wide_type_t mat_wt;

        for (unsigned int j = 0; j < INPUT_MAT_PAR_ENTRIES; ++j)
        {
            mat_wt[j] = in_wt.range(INPUT_MAT_DATATYPE_SIZE * (j + 1) - 1, INPUT_MAT_DATATYPE_SIZE * j);
        }
        out_stream.write(mat_wt);
    }
}

void readVecToBuffer(
    const unsigned int num_cols,
    const ap_int<DATA_WIDTH> *in_ptr,
    ap_int<DATA_WIDTH> *buffer_ptr)
{
    const unsigned int num_blocks = num_cols / INPUT_VEC_PAR_ENTRIES;

read_vec:
    for (unsigned int i = 0; i < num_blocks; ++i)
    {
#pragma HLS PIPELINE
        buffer_ptr[i] = in_ptr[i];
    }
}

void writeStreamToVec(
    const unsigned int num_rows,
    hls::stream<output_wide_type_t> &in_stream,
    ap_int<DATA_WIDTH> *out_ptr)
{
    double rows = num_rows;
    const unsigned int num_blocks = hls::ceil(rows / OUTPUT_PAR_ENTRIES);

write_vec:
    for (unsigned int i = 0; i < num_blocks; ++i)
    {
#pragma HLS PIPELINE
        out_ptr[i] = in_stream.read();
    }
}

#endif