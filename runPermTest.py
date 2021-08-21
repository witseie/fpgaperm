import sys, getopt
import argparse
import subprocess
import os.path
import statsmodels.api as sm
import numpy as np
import pandas as pd

def has_header(file):
    with open(file) as f:
        header = f.readline().split()
        has_header = pd.Series(pd.to_numeric(header, errors='coerce')).isnull().all()

        return has_header

def get_pheno_residuals(fam_file, covar_file = None):
    if has_header(fam_file):
        pheno_df = pd.read_csv(fam_file, delimiter = " ")
    else:
        pheno_df = pd.read_csv(fam_file, delimiter = " ", header=None)

    pheno_data = pheno_df[5]

    if (pheno_data.isnull().values.any()):
        print('Error: Missing phenotypes detected')
        sys.exit(2)

    num_rows = pheno_data.shape[0]

    covar_data = pd.DataFrame(index=range(num_rows))

    if covar_file is not None:
        if has_header(covar_file):
            covar_df = pd.read_csv(covar_file, delimiter = " ")
        else:
            covar_df = pd.read_csv(covar_file, delimiter = " ", header=None)

        if covar_df.shape[0] != num_rows:
            print('Error: Invalid covariate file')
            sys.exit(2)

        covar_data = covar_df

        if (covar_data.isnull().values.any()):
            print('Error: Missing covariates detected')
            sys.exit(2)

        covar_data.insert(loc=0, value=1, column='y_int')
    else:
        covar_data['y_int'] = 1

    # Check if the data is categorical
    is_categorical = np.in1d([0,1], pheno_data).all()

    if is_categorical:
        glm_res = sm.GLM(pheno_data, covar_data, family=sm.families.Binomial()).fit()
        residuals = glm_res.resid_working
    else:
        # Check if the phenotype data is centred
        if np.isclose(np.mean(pheno_data), 0, atol=1e-3) and covar_file is None:
            residuals = pheno_data
        else:
            ols_res = sm.OLS(pheno_data, covar_data).fit()
            residuals = ols_res.resid

    max_val = np.max(np.abs(residuals))
    min_val = np.min(np.abs(residuals))

    # Get the phenotype scaling factor to pass to the host application
    scaling_factor = np.round(100 / max_val, 2)

    location = fam_file.rsplit('/', 1)[0]

    resids_file = location + "/pheno_residuals"
    np.savetxt(resids_file, residuals, fmt='%.10f')

    return resids_file, scaling_factor

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_data', required=True)
    parser.add_argument('-p', '--perm_algo', nargs='+')
    parser.add_argument('-c', '--covar', nargs='?')

    args = parser.parse_args()

    host = '/home/centos/FPGA_perm/host'
    xclbin = '/home/centos/FPGA_perm/mv_mul.hw.xilinx_vu9p.awsxclbin'

    if not os.path.isfile(host):
        print('Error: Cannot find the host application executable ' + host)
        sys.exit(2)

    if not os.path.isfile(xclbin):
        print('Error: Cannot find the FPGA config file ' + xclbin)
        sys.exit(2)

    input_pattern = args.input_data

    bed_file = '{0}.bed'.format(input_pattern)
    fam_file = '{0}.fam'.format(input_pattern)
    bim_file = '{0}.bim'.format(input_pattern)

    if not os.path.isfile(bed_file):
        print('Error: invalid .bed file')
        sys.exit(2)
    else:
        print('Genotype file: ' + bed_file)

    if not os.path.isfile(fam_file):
        print('Error: invalid .fam file')
        sys.exit(2)
    else:
        print('Phenotype file: ' + fam_file)

    if not os.path.isfile(bim_file):
        print('Error: invalid .bim file')
        sys.exit(2)
    else:
        print('.bim file: ' + bim_file)

    covar_file = args.covar
    if covar_file is not None:
        if not os.path.isfile(covar_file):
            print('Error: invalid covariate file')
            sys.exit(2)
        else:
            print('Covariate file: ' + covar_file)

    resids_file, phenotype_scaling_factor = get_pheno_residuals(fam_file, covar_file)

    if args.perm_algo:
        if args.perm_algo[0] == 'maxT':
            if not args.perm_algo[1].isdigit():
                print('Error: ', args.perm_algo[1], ' is not a number. Usage: -p maxT [num_perms]')
                sys.exit(2)

            proc = subprocess.run([host,
                                '-xclbin', xclbin,
                                '-perm_algo', args.perm_algo[0], args.perm_algo[1],
                                '-phenotype_scale_factor', str(phenotype_scaling_factor),
                                '-input_data', bed_file, resids_file, bim_file],
                                stdout=sys.stdout)
        elif args.perm_algo[0] == 'adp':
            if not args.perm_algo[1].isdigit():
                print('Error: ', args.perm_algo[1], ' is not a number. Usage: -p adp [min_perms] [max_perms]')
                sys.exit(2)

            if not args.perm_algo[2].isdigit():
                print('Error: ', args.perm_algo[2], ' is not a number. Usage: -p adp [min_perms] [max_perms]')
                sys.exit(2)

            proc = subprocess.run([host,
                                '-xclbin', xclbin,
                                '-perm_algo', args.perm_algo[0], args.perm_algo[1], args.perm_algo[2],
                                '-phenotype_scale_factor', str(phenotype_scaling_factor),
                                '-input_data', bed_file, resids_file, bim_file],
                                stdout=sys.stdout)
    else:
        proc = subprocess.run([host,
                               '-xclbin', xclbin,
                               '-phenotype_scale_factor', str(phenotype_scaling_factor),
                               '-input_data', bed_file, resids_file, bim_file],
                              stdout=sys.stdout)

if __name__ == '__main__':
    main()
