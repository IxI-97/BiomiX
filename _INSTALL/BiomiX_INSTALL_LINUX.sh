#!/bin/bash

conda env list | grep -q '^BiomiX-env' && conda env remove -n BiomiX-env --yes

source ~/miniconda3/etc/profile.d/conda.sh

conda create -n BiomiX-env python=3.9 -y
conda activate BiomiX-env
conda config --add channels conda-forge
conda install mamba -y
mamba install -c conda-forge r-base=4.4.1 -y
mamba install -c conda-forge r-systemfonts=1.1.0 -y
mamba install -c conda-forge r-rcpp=1.0.13 -y
python3 BiomiX_INSTALL_LINUX.py
