#include "fpga_types.hpp"
#include "data_transfer.hpp"
#include "gwas_mv_mul.hpp"

extern "C"
{
    void krnl_mvmul_3(
        const unsigned int num_rows,
        const unsigned int num_cols,
        const ap_int<DATA_WIDTH> *mean,
        const ap_int<DATA_WIDTH> *mat,
        const ap_int<DATA_WIDTH> *vec,
        ap_int<DATA_WIDTH> *result)
    {
#pragma HLS INTERFACE m_axi port = mat offset = slave bundle = gmem0
#pragma HLS INTERFACE m_axi port = vec offset = slave bundle = gmem1
#pragma HLS INTERFACE m_axi port = mean offset = slave bundle = gmem2
#pragma HLS INTERFACE m_axi port = result offset = slave bundle = gmem3

#pragma HLS INTERFACE s_axilite port = num_rows bundle = control
#pragma HLS INTERFACE s_axilite port = num_cols bundle = control
#pragma HLS INTERFACE s_axilite port = mean bundle = control
#pragma HLS INTERFACE s_axilite port = mat bundle = control
#pragma HLS INTERFACE s_axilite port = vec bundle = control
#pragma HLS INTERFACE s_axilite port = result bundle = control
#pragma HLS INTERFACE s_axilite port = return bundle = control

#pragma HLS DATAFLOW

        hls::stream<mean_wide_type_t> meanStream;
#pragma HLS data_pack variable = meanStream

        hls::stream<matrix_wide_type_t> matStream;
#pragma HLS data_pack variable = matStream

        ap_int<DATA_WIDTH> vecBuffer[NUM_COLS / INPUT_VEC_PAR_ENTRIES];
#pragma HLS ARRAY_PARTITION variable = vecBuffer cyclic factor = 8 dim = 1

        hls::stream<output_wide_type_t> outVecStream;
#pragma HLS data_pack variable = outVecStream

#pragma HLS DATAFLOW
        readVecToBuffer(num_cols, vec, vecBuffer);
        readMeanToStream(num_rows, mean, meanStream);
        readMatToStream(num_rows, num_cols, mat, matStream);

        mv_mul(num_rows, num_cols, meanStream, matStream, vecBuffer, outVecStream);

        writeStreamToVec(num_rows, outVecStream, result);
    }
}