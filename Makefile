.PHONY: help

# Run Target:
#   hw  - Compile for hardware
#   sw_emu/hw_emu - Compile for software/hardware emulation

help::
	@echo "Makefile Usage:"
	@echo "  make all TARGET=<sw_emu/hw_emu/hw>"
	@echo "      Build xclbin for target FPGA (default hw) and compile host app."
	@echo ""
	@echo "  make xclbin TARGET=<sw_emu/hw_emu/hw>"
	@echo "      Build xclbin for target FPGA (default hw)."
	@echo ""
	@echo "  make clean"
	@echo "      Remove files."
	@echo ""
	@echo "  make app"
	@echo "      Compile the host app."
	@echo ""
	@echo "  make app-debug"
	@echo "      Compile the host app with debug prints."
	@echo ""

CURRENT_DIR = $(shell pwd)

#Checks for XILINX_VITIS
ifndef XILINX_VITIS
$(error XILINX_VITIS variable is not set, please set correctly and rerun)
endif

TARGET := hw

DSA := xilinx_aws-vu9p-f1_shell-v04261818_201920_2

CXX := g++
VPP := v++

BUILD_DIR := fpga_gwas_mv_mul/mv_mul.$(TARGET).xilinx_aws_vu9p_f1

HOST_SRC = host_app/host.cpp

HOST_EXE = host

KERNEL_SRC_0 := fpga_gwas_mv_mul/src/krnl_mv_mul_0.cpp
KERNEL_SRC_1 := fpga_gwas_mv_mul/src/krnl_mv_mul_1.cpp 
KERNEL_SRC_2 := fpga_gwas_mv_mul/src/krnl_mv_mul_2.cpp 
KERNEL_SRC_3 := fpga_gwas_mv_mul/src/krnl_mv_mul_3.cpp 

KERNEL_XO_0 := $(BUILD_DIR)/mv_mul_0.$(TARGET).xilinx_aws_vu9p.xo
KERNEL_XO_1 := $(BUILD_DIR)/mv_mul_1.$(TARGET).xilinx_aws_vu9p.xo
KERNEL_XO_2 := $(BUILD_DIR)/mv_mul_2.$(TARGET).xilinx_aws_vu9p.xo
KERNEL_XO_3 := $(BUILD_DIR)/mv_mul_3.$(TARGET).xilinx_aws_vu9p.xo

XCLBIN := mv_mul.$(TARGET).xilinx_aws_vu9p.xclbin

# g++ compiler arguments
CXXFLAGS = -I$(XILINX_XRT)/include -I$(XILINX_VIVADO)/include/
CXXFLAGS += -Ihost_app/ -Ihost_app/include/ -O3 -std=c++11 -fopenmp
LDFLAGS = -lOpenCL -lpthread -lrt -lstdc++ -lgomp -L$(XILINX_XRT)/lib/

# v++ compiler arguments
VPP_COMPILE_OPTS := -t $(TARGET) -I'fpga_gwas_mv_mul/src/'

.PHONY: xclbin app all

xclbin: $(XCLBIN)

app: $(HOST_EXE)

app-debug: CXXFLAGS += -DDEBUG
app-debug: $(HOST_EXE)

all: xclbin app

clean:
	rm -rf $(HOST_EXE) $(XCLBIN)* $(BUILD_DIR) v++_* x*.log

XCLBIN_XO =
ifeq ($(TARGET), hw)
	VPP_LINK_OPTS := -t $(TARGET) --config fpga_gwas_mv_mul/link.cfg
	XCLBIN_XO = $(KERNEL_XO_0) $(KERNEL_XO_1) $(KERNEL_XO_2) $(KERNEL_XO_3)
else
	VPP_LINK_OPTS := -t $(TARGET) --config fpga_gwas_mv_mul/link_emu.cfg
	XCLBIN_XO = $(KERNEL_XO_0)
endif

$(XCLBIN): $(XCLBIN_XO)
	$(VPP) $(VPP_LINK_OPTS) -l -o $@ $+ --temp_dir $(BUILD_DIR)

$(KERNEL_XO_0): $(KERNEL_SRC_0)
	$(VPP) $(VPP_COMPILE_OPTS) --config fpga_gwas_mv_mul/kernel0.cfg -c -k krnl_mvmul_0 -o $@ $(KERNEL_SRC_0) --temp_dir $(BUILD_DIR)

$(KERNEL_XO_1): $(KERNEL_SRC_1)
	$(VPP) $(VPP_COMPILE_OPTS) --config fpga_gwas_mv_mul/kernel1.cfg -c -k krnl_mvmul_1 -o $@ $(KERNEL_SRC_1) --temp_dir $(BUILD_DIR)

$(KERNEL_XO_2): $(KERNEL_SRC_2)
	$(VPP) $(VPP_COMPILE_OPTS) --config fpga_gwas_mv_mul/kernel2.cfg -c -k krnl_mvmul_2 -o $@ $(KERNEL_SRC_2) --temp_dir $(BUILD_DIR)

$(KERNEL_XO_3): $(KERNEL_SRC_3)
	$(VPP) $(VPP_COMPILE_OPTS) --config fpga_gwas_mv_mul/kernel3.cfg -c -k krnl_mvmul_3 -o $@ $(KERNEL_SRC_3) --temp_dir $(BUILD_DIR)

$(HOST_EXE): $(HOST_SRC)
	$(CXX) $(HOST_SRC) $(HOST_INC) $(LDFLAGS) $(CXXFLAGS) -o '$(HOST_EXE)' 
