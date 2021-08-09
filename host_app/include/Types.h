#ifndef TYPES_H
#define TYPES_H

#include <vector>

#include <ap_int.h>
#include <ap_fixed.h>

#include "CL/opencl.h"
#include <CL/cl2.hpp>

#define DATA_WIDTH 512

#define INPUT_MAT_DATATYPE_SIZE 2
#define INPUT_MAT_PAR_ELEMS     (DATA_WIDTH / INPUT_MAT_DATATYPE_SIZE) // (512/2 = 256)

#define INPUT_VEC_DATATYPE_SIZE       16
#define INPUT_VEC_DATATYPE_SIZE_BYTES (INPUT_VEC_DATATYPE_SIZE >> 3)
#define INPUT_VEC_PAR_ELEMS           (DATA_WIDTH / INPUT_VEC_DATATYPE_SIZE) // (512/16 = 32)

#define INPUT_MEAN_DATATYPE_SIZE       16
#define INPUT_MEAN_DATATYPE_SIZE_BYTES (INPUT_MEAN_DATATYPE_SIZE >> 3)
#define INPUT_MEAN_PAR_ELEMS           (DATA_WIDTH / INPUT_MEAN_DATATYPE_SIZE) // (512/16 = 32)

#define OUTPUT_DATATYPE_SIZE       32
#define OUTPUT_DATATYPE_SIZE_BYTES (OUTPUT_DATATYPE_SIZE >> 3)
#define OUTPUT_PAR_ELEMS           (DATA_WIDTH / OUTPUT_DATATYPE_SIZE) // (512/32 = 16)

#define input_matrix_data_type_t ap_int<2>
#define input_mean_data_type_t   ap_fixed<16, 2>   // Range: -4      to 3.99987793          Resolution = 0.00012207
#define input_pheno_data_type_t  ap_fixed<16, 11>  // Range: -1024   to 1023.96875          Resolution = 0.03125
#define output_data_type_t       ap_fixed<32, 20>  // Range: -524288 to 524287.999755859    Resolution = 0.000244141

#define input_matrix_vector_t std::vector<unsigned char, AlignedAllocator<unsigned char>>
#define input_mean_vector_t   std::vector<input_mean_data_type_t, AlignedAllocator<input_mean_data_type_t>>
#define input_pheno_vector_t  std::vector<input_pheno_data_type_t, AlignedAllocator<input_pheno_data_type_t>>
#define output_vector_t       std::vector<output_data_type_t, AlignedAllocator<output_data_type_t>>

// 4kB aligned memory allocator for efficient memory transfers
template <typename T>
struct AlignedAllocator
{
    using value_type = T;
    T *allocate(std::size_t num)
    {
        void *ptr = nullptr;
        if (posix_memalign(&ptr, 4096, num * sizeof(T)))
            throw std::bad_alloc();
        return reinterpret_cast<T *>(ptr);
    }
    void deallocate(T *p, std::size_t num)
    {
        free(p);
    }
};

enum perm_algo_t
{
    perm_adaptive = 0,
    perm_maxT = 1,
    perm_regression_only = 2
};

struct perm_adp_results_t
{
    std::vector<unsigned int> perm_count;
    std::vector<unsigned int> perm_dropped;
    std::vector<char> is_snp_dropped;

    perm_adp_results_t(int num_rows = 0)
    {
        perm_count = std::vector<unsigned int>(num_rows, 0);
        perm_dropped = std::vector<unsigned int>(num_rows, 0);
        is_snp_dropped = std::vector<char>(num_rows, 0);
    }
};

struct ocl_task_buffers_t
{
    cl::Buffer mean_buf;
    cl::Buffer mat_buf;
    cl::Buffer pheno_buf;
    cl::Buffer out_buf;
};

struct mean_stddev_t
{
    size_t count;
    double sq_sum;
    double mean;
    double std_dev;
};

struct vector_stats_t
{
    mean_stddev_t Y;
    std::vector<mean_stddev_t> X;

    vector_stats_t(int num_rows = 0)
    {
        X = std::vector<mean_stddev_t>(num_rows);
    }
};

#endif
