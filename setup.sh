export PATH=$PATH:/usr/local/sbin
export PATH=$PATH:/usr/sbin
export FPGA_PERM_PATH=/home/centos/FPGA_perm

pip3 install numpy --user
pip3 install pandas --user
pip3 install statsmodels --user

git clone https://github.com/aws/aws-fpga.git $AWS_FPGA_REPO_DIR
git clone https://github.com/witseie/fpgaperm.git $FPGA_PERM_PATH # Check this

source $AWS_FPGA_REPO_DIR/vitis_setup.sh
source $AWS_FPGA_REPO_DIR/vitis_runtime_setup.sh

export PLATFORM_REPO_PATHS=$AWS_FPGA_REPO_DIR/Vitis/aws_platform/xilinx_aws-vu9p-f1_shell-v04261818_201920_2

sudo systemctl start mpd
