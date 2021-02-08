export PATH=$PATH:/usr/local/sbin
export PATH=$PATH:/usr/sbin

if ! [ -x "$(command -v python3)" ]; then
  sudo yum install -y python3
else
  echo "python3 available"
fi

pip3 install numpy --user
pip3 install matplotlib --user

git clone https://github.com/aws/aws-fpga.git $AWS_FPGA_REPO_DIR

source $AWS_FPGA_REPO_DIR/vitis_setup.sh
source $AWS_FPGA_REPO_DIR/vitis_runtime_setup.sh

export PLATFORM_REPO_PATHS=$AWS_FPGA_REPO_DIR/Vitis/aws_platform/xilinx_aws-vu9p-f1_shell-v04261818_201920_2

sudo systemctl start mpd
