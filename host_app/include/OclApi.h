#ifndef __APIHANDLE_H__
#define __APIHANDLE_H__

#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY 1
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include "CL/opencl.h"
#include <CL/cl2.hpp>

#include "Types.h"

#define OCL_CHECK(error, call)                                       \
    call;                                                            \
    if (error != CL_SUCCESS)                                         \
    {                                                                \
        printf("%s:%d Error calling " #call ", error code is: %d\n", \
               __FILE__, __LINE__, error);                           \
        exit(EXIT_FAILURE);                                          \
    }

class OclApi
{
    char *read_binary_file(const std::string &xclbin_file_name, unsigned &nb)
    {
        std::cout << "Reading " << xclbin_file_name << std::endl;

        if (access(xclbin_file_name.c_str(), R_OK) != 0)
        {
            printf("ERROR: %s xclbin not found\n", xclbin_file_name.c_str());
            exit(EXIT_FAILURE);
        }

        // Load XCL Bin into char buffer
        std::ifstream bin_file(xclbin_file_name.c_str(), std::ifstream::binary);
        bin_file.seekg(0, bin_file.end);
        nb = bin_file.tellg();
        bin_file.seekg(0, bin_file.beg);
        char *buf = new char[nb];
        bin_file.read(buf, nb);
        return buf;
    }

public:
    cl::Context context;
    cl::CommandQueue queue;
    cl::Program program;
    std::vector<cl::Device> devices;
    cl::Device device;

    OclApi(std::string xclbin_file_name)
    {
        cl_int err;
        std::vector<cl::Platform> platforms;

        // Get OCL platforms
        OCL_CHECK(err, err = cl::Platform::get(&platforms));

        size_t i;
        cl::Platform platform;
        for (i = 0; i < platforms.size(); i++)
        {
            platform = platforms[i];

            OCL_CHECK(err, std::string platformName = platform.getInfo<CL_PLATFORM_NAME>(&err));

            if (platformName == "Xilinx")
            {
                // std::cout << "Found Platform" << std::endl;
                std::cout << "Found platform: " << platformName.c_str() << std::endl;
                break;
            }
        }
        if (i == platforms.size())
        {
            std::cout << "Platform Not Found" << std::endl;
            exit(err);
        }

        //Get FPGA devices and select 1st device
        OCL_CHECK(err, err = platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devices));

        device = devices[0];
        OCL_CHECK(err, std::string cl_device_name = device.getInfo<CL_DEVICE_NAME>(&err));

        std::cout << "Using device: " << cl_device_name.c_str() << std::endl;

        char *krnl_bin;
        unsigned krnl_bin_size;

        //Create context and command queue
        OCL_CHECK(err, context = cl::Context(device, NULL, NULL, NULL, &err));

        // Load bitstream from XCLbin
        krnl_bin = read_binary_file(xclbin_file_name, krnl_bin_size);

        cl::Program::Binaries bins{{krnl_bin, krnl_bin_size}};

        devices.resize(1);
        OCL_CHECK(err, program = cl::Program(context, {device}, bins, NULL, &err));

        // Create Out of Order command queue
        OCL_CHECK(err, queue = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err));

        delete[] krnl_bin;

        std::cout << "FPGA Setup Complete" << std::endl
                  << std::endl;
    }
};

class OclTask
{
    std::vector<cl::Event> write_data_event = std::vector<cl::Event>(1);
    std::vector<cl::Event> enqueue_event = std::vector<cl::Event>(1);
    std::vector<cl::Event> read_data_event = std::vector<cl::Event>(1);

    cl_int err;

public:
    std::vector<cl::Event> *getDoneEv() { return &read_data_event; }

    void run(
        OclApi &api,
        cl::Kernel kernel,
        ocl_task_buffers_t task_buffers,
        std::vector<cl::Event> *prevEvent = nullptr)
    {
        // Copy data to the FPGA
        OCL_CHECK(err, err = api.queue.enqueueMigrateMemObjects({task_buffers.mat_buf, task_buffers.pheno_buf, task_buffers.mean_buf},
                                                                0 /* 0 means from host*/,
                                                                prevEvent, &write_data_event[0]));

        // Enqueue kernel execution
        OCL_CHECK(err, err = api.queue.enqueueTask(kernel, &write_data_event, &enqueue_event[0]));

        // Copy result from FPGA
        OCL_CHECK(err, err = api.queue.enqueueMigrateMemObjects({task_buffers.out_buf},
                                                                CL_MIGRATE_MEM_OBJECT_HOST,
                                                                &enqueue_event, &read_data_event[0]));
    }
};

#endif
