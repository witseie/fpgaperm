#ifndef FPGA_TYPES_H
#define FPGA_TYPES_H

#include "ap_int.h"
#include "ap_fixed.h"

#include "hls_stream.h"
#include "wide_type.hpp"

#define NUM_COLS   262144

#define DATA_WIDTH 512

#define INPUT_VEC_DATATYPE_SIZE  16
#define INPUT_MEAN_DATATYPE_SIZE 16
#define INPUT_MAT_DATATYPE_SIZE  2
#define OUTPUT_DATATYPE_SIZE     32

#define INPUT_VEC_PAR_ENTRIES    (DATA_WIDTH / INPUT_VEC_DATATYPE_SIZE)
#define INPUT_MEAN_PAR_ENTRIES   (DATA_WIDTH / INPUT_MEAN_DATATYPE_SIZE)
#define INPUT_MAT_PAR_ENTRIES    (DATA_WIDTH / INPUT_MAT_DATATYPE_SIZE)
#define OUTPUT_PAR_ENTRIES       (DATA_WIDTH / OUTPUT_DATATYPE_SIZE)

#define input_matrix_data_type_t ap_uint<2>
#define input_mean_data_type_t   ap_fixed<16, 2>
#define input_vector_data_type_t ap_fixed<16, 11>
#define output_data_type_t       ap_fixed<32, 20, AP_TRN, AP_SAT>

#define matrix_wide_type_t xf::blas::WideType<input_matrix_data_type_t, INPUT_MAT_PAR_ENTRIES, INPUT_MAT_DATATYPE_SIZE>
#define mean_wide_type_t   xf::blas::WideType<input_mean_data_type_t, INPUT_MEAN_PAR_ENTRIES>
#define vector_wide_type_t xf::blas::WideType<input_vector_data_type_t, INPUT_VEC_PAR_ENTRIES>
#define output_wide_type_t xf::blas::WideType<output_data_type_t, OUTPUT_PAR_ENTRIES>

#define int_snp_data_type_t ap_fixed<18, 2>
#define int_vec_data_type_t ap_fixed<18, 11>
#define acc_data_type_t     ap_fixed<48, 30>
#define mul_data_type_t     ap_fixed<48, 30>

#endif