import numpy as np
import matplotlib.pyplot as plt
import sys, getopt
import argparse
import subprocess

def get_pheno_residuals(phenoFile, covarFile = None):
    with open(phenoFile) as f:
        lines = f.read().splitlines()

    lines = [line.split() for line in lines]

    pheno = np.array(lines)
    pheno = (pheno[:, -1:]).astype('float32')

    num_rows = np.shape(pheno)[0]

    if covarFile is not None:
        with open(covarFile) as f:
            lines = f.read().splitlines()

        lines.pop(0)
        lines = [line.split() for line in lines]

        covars = np.array(lines)
        covars = (covars[:, 2:]).astype('float32')

        ones_col = np.ones((num_rows , 1)).astype('float32')

        covars = np.concatenate((ones_col, covars), axis=1)

        covars_pinv = np.linalg.pinv(covars)

        fit = np.dot(covars,(np.dot(covars_pinv, pheno)))

        pheno_resids = pheno - fit
    else:
        pheno_mean = np.mean(pheno)
        pheno_resids = pheno - pheno_mean

    location = phenoFile.rsplit('/',1)[0]

    resids_file = location + "/pheno_resids"
    np.savetxt(resids_file, pheno_resids, fmt='%.10f')

    return resids_file

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_data', nargs=3)
    parser.add_argument('-r', '--num_rows')
    parser.add_argument('-p', '--perm_algo', nargs=2)
    parser.add_argument('-c', '--covar')
    parser.add_argument('-x', '--xclbin')
    parser.add_argument('-b', '--buf_size')

    parser.add_argument('-a', '--host_app')

    args = parser.parse_args()

    host = './host'
    geno_file = [s for s in args.input_data if '.bed' in s][0]
    pheno_file = [s for s in args.input_data if '.fam' in s][0]
    snp_file = [s for s in args.input_data if '.bim' in s][0]

    input_files = [geno_file, pheno_file, snp_file]

    if any(input_files) is None:
        print('Missing input file')
        sys.exit(2)

    covar_file = args.covar

    resids_file = get_pheno_residuals(pheno_file, covar_file)

    proc = subprocess.run([host,
                            '-xclbin', args.xclbin,
                            #'-num_rows', args.num_rows,
                            '-perm_algo', args.perm_algo[0], args.perm_algo[1],
                            '-buf_size', args.buf_size,
                            '-input_data', geno_file, resids_file, snp_file], stdout=sys.stdout)

    if args.perm_algo[0] == 'maxT':
        host_location = host.rsplit('/',1)[0]

        out_file = host_location + '/perm.best.txt'
        with open(out_file) as f:
            lines = f.read().splitlines()
        
        lines = [line.split() for line in lines]
        results = np.array(lines)
        results = (results[1:, -1:]).astype('float32')
        perc_95 = np.percentile(results, 95)
        perc_99 = np.percentile(results, 99)

        plt.hist(results, bins=30)
        plt.axvline(x=perc_95,color='k', linestyle='--')
        plt.axvline(x=perc_99,color='k', linestyle='--')

        # plt.show()
        plt.savefig('hist.png', bbox_inches='tight')

    # out_file = '/home/yaniv/Desktop/mv_mul_6_unconstrained_6CU_32KB_buf/Results/Adaptive/1000000 Perms/2MB_min_buf_41_56.txt'
    # x = []
    # y = []
    # with open(out_file,"r") as f:
    #     for ln in f:
    #         if ln.startswith("Perm = "):
    #             test = ln.split()
    #             x.append(test[2])
    #             y.append(test[6])
    #             asd = 1

    # x = np.asarray(x, dtype=np.int32)
    # y = np.asarray(y, dtype=np.int32)
    # plt.plot(x, y, 'b.')
    # plt.show()





