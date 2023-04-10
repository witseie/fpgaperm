export PATH=$PATH:/usr/local/sbin
export PATH=$PATH:/usr/sbin
export FPGA_PERM_PATH=/home/centos/FPGA_perm
export LC_ALL="C"

pip3 install --user statsmodels

git clone https://github.com/aws/aws-fpga.git $AWS_FPGA_REPO_DIR

source $AWS_FPGA_REPO_DIR/vitis_setup.sh
source $AWS_FPGA_REPO_DIR/vitis_runtime_setup.sh

export PLATFORM_REPO_PATHS=$AWS_FPGA_REPO_DIR/Vitis/aws_platform/xilinx_aws-vu9p-f1_shell-v04261818_201920_2

# Clean up environment variables
PATH="$(perl -e 'print join(":", grep { not $seen{$_}++ } split(/:/, $ENV{PATH}))')"
PYTHONPATH="$(perl -e 'print join(":", grep { not $seen{$_}++ } split(/:/, $ENV{PYTHONPATH}))')"
LD_LIBRARY_PATH="$(perl -e 'print join(":", grep { not $seen{$_}++ } split(/:/, $ENV{LD_LIBRARY_PATH}))')"

sudo systemctl start mpd
