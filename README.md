# fpga-perm

An FPGA-accelerated implementation of GWAS permutation testing designed for AWS EC2 F1 (`f1.2xlarge`) instances.

Two different permutation testing algorithms are supported:
- maxT permutation testing
- Adaptive permutation testing

The accelerator consists of three components:
- The FPGA hardware component which accelerates the compute-heavy part of the algorithm
- A C++ host application which has a numer of functions including:
    - Data preprocessing
    - Managing the FPGA execution
    - Running the permutation testing algorithms
- A Python script (`runPermTest.py`) which calculates phenotype residuals using the `statsmodels` package and runs the host application

**NB: The accelerator only supports PLINK-formatted input files (.bed, .bim, .fam)**

**NB: The accelerator only supports GWAS permutation testing using the additive genotypic model (i.e. it expects variant genotypes encoded as 0, 1, 2 or missing using PLINK's .bed file format)**

**NB: The host application does not filter the input data so all data in the supplied files is used for permutation testing**

## Instance Setup
A Linux CentOS (username=`centos`) AMI is provided for convenience. The working directory of FPGA_perm is `/home/centos/FPGA_perm` and the setup script `/home/centos/FPGA_perm/setup.sh` should be run to initialise the FPGA_perm environment.

## Usage
The `runPermTest.py` script should be used to run FPGA_perm on the selected dataset.

The script expects the host application executable (`host`) and the FPGA binary (`mv_mul_4CU.hw.xilinx_aws_vu9p.awsxclbin`) to be in the FPGA_perm working directory.

The Python script uses `numpy`, `pandas` and `statsmodels`

### Input Data
`-i [input_pattern]` or `--input_data [input_pattern]` is used to reference the input files i.e. `input_pattern`.bed, `input_pattern`.bim, `input_pattern`.fam

As the accelerator performs association tests using simple linear regression, only continuous phenotypes are supported. If case/control (0 = control, 1 = case) phenotypes are detected in the .fam file, they are linearised by `runPermTest.py` using a general linear model.

Missing phenotypes are not supported.

Multiple phenotypes are not supported i.e. only the first column of phenotypes in the .fam file is used

### Covariate Data
`-c [filename]` or `--covar [filename]` is an optional argument used to reference covariate data stored in `filename`

The order of the covariates should correspond to the order of the phenotypes in the .fam file.

All of the covariate data in `filename` is used.

In order to include covariates in the permutation procedure, the residuals of the phenotype data regressed on the covariate data (with an appended column of ones representing the y-intercept) are used instead of the phenotype data itself.

Missing covariates are not supported.

### Permutation Testing
`-p` or `--perm_algo` is used to select the permutation testing algorithm

#### Adaptive Permutation Testing
To select adaptive permutation testing, use the following:
```
-p adp [min_perms] [max_perms]
```
- `min_perms` is the minimum permutations per SNP
- `max_perms` is the maximum permutations per SNP

The host application will generate a text file (`perm_adaptive.res.txt`) with the following fields:

- SNP ID
- Observed test statistic
- Permutation p-value
- Number of permutations

#### maxT Permutation Testing
To select maxT permutation testing, use the following:
```
-p maxT [num_perms]
```
- `num_perms` is the number of permutations

The host application will generate a text file (`perm_maxT.res.txt`) with the following fields:

- SNP ID
- Observed test statistic
- Permutation p-value

#### Association Testing
If the `--perm_algo` argument is not supplied to `runPermTest.py`, basic association testing is performed (although this has not been extensively tested)

The host application will generate a text file (`gwas_results.txt`) with the following fields:

- SNP ID
- Number of non-missing genotypes
- Linear regression coefficient (beta)
- Test statistic

## Make
A Makefile is provided to compile the host application (using GCC) and the FPGA design (using Xilinx Vitis/Vivado 2020.2).

In order to compile the host application, the OpenCL, OpenMP and Xilinx XRT libraries must be installed

The Makefile allows an FPGA design to be compiled for software emulation (sw_emu), hardware emulation (hw_emu) or hardware (hw). The FPGA design is specifically compiled for the Xilinx FPGA provided by the `f1.2xlarge` instance (the Xilinx VU9P).

To compile the host application use

    make app

To compile the FPGA design use

    make xclbin