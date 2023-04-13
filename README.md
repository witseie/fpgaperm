# fpga-perm

An FPGA-accelerated implementation of GWAS permutation testing designed for AWS EC2 F1 (`f1.2xlarge`) instances.

Two different permutation testing algorithms are supported:
- maxT permutation testing
- Adaptive permutation testing

The accelerator consists of three components:
- The FPGA hardware component which accelerates the compute-heavy part of the algorithm
- A C++ host application which has a number of functions including:
    - Data preprocessing
    - Managing the FPGA execution
    - Running the permutation testing algorithms
- A Python 3 script (`runPermTest.py`) which calculates phenotype residuals and runs the host application

**NB: The accelerator only supports PLINK-formatted input files (.bed, .bim, .fam)**

**NB: The accelerator only supports GWAS permutation testing using the additive genotypic model (i.e. it expects variant genotypes encoded as 0, 1, 2 or missing using PLINK's .bed file format)**

**NB: The host application does not filter the input data so all data in the supplied files is used for permutation testing**

## Installation


The accelerator is designed to run on AWS F1 FPGA instances launched with the [FPGA Developer AMI](https://aws.amazon.com/marketplace/pp/prodview-iehshpgi7hcjg?sr=0-2&ref_=beagle&applicationId=AWSMPContessa). Version 1.0.5 of this AMI does not require subscription (e.g., ami-056e4346b21bf5cb1 in us-east-1, ami-0d5a2e5dd1e802d8e in af-south-1, ami-0bda4348cc5912e6f in eu-west-1). The code should run on later versions of the AMI but this has not been fully.

- Launch an instance of the appropriate AMI
- To initialise the FPGA_perm environment, source the setup script `setup.sh`.

## Usage
Once the FPGA_perm environment has been initialised, the `runPermTest.py` Python3 script can be used to run FPGA_perm on the selected dataset e.g.
```
python3 ./runPermTest.py -i [input_data] -p maxT [num_perms]
```

### Input Data
`-i [input_pattern]` or `--input_data [input_pattern]` is used to reference the input files i.e. `input_pattern`.bed, `input_pattern`.bim, `input_pattern`.fam

As the accelerator performs association tests using simple linear regression, only continuous phenotypes are supported. If case/control (0 = control, 1 = case) phenotypes are detected in the .fam file, they are linearised by `runPermTest.py` using a general linear model.

Missing phenotypes are not supported.

Multiple phenotypes are not supported i.e. only the first column of phenotypes in the .fam file is used

The Python script reads the phenotype data from the 5th column of the .fam file.

### Covariate Data
`-c [filename]` or `--covar [filename]` is an optional argument used to reference covariate data stored in `filename`.

In order to include covariates in the permutation procedure, the residuals of the phenotype data regressed on the covariate data (with an appended column of ones representing the y-intercept) are used instead of the phenotype data itself.

The covariate data should be stored in space-separated columns with the first column representing the first covariate.

The order of the covariates should correspond to the order of the phenotypes in the .fam file.

All of the covariate data in `filename` is used. 

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

## Requirements
The `runPermTest.py` script expects the host application executable (`host`) and the FPGA binary (`mv_mul.hw.xilinx_aws_vu9p.awsxclbin`) to be in the FPGA_perm working directory.

The `runPermTest.py` script uses `numpy`, `pandas` and `statsmodels`.

## Make
A Makefile is included to compile the host application (using GCC) and the FPGA design (using Xilinx Vitis/Vivado 2020.2).

In order to compile the host application, the OpenCL, OpenMP and Xilinx XRT libraries must be installed.

The Makefile allows an FPGA design to be compiled for software emulation (sw_emu), hardware emulation (hw_emu) or hardware (hw). The FPGA design is specifically compiled for the Xilinx FPGA provided by the `f1.2xlarge` instance (the Xilinx Virtex UltraScale+ VU9P).

To compile the host application use

    make app

To compile the FPGA design use

    make xclbin
