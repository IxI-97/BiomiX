#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda create -n BiomiX-env python=3.9 -y
conda activate BiomiX-env
conda config --add channels conda-forge
conda install -c conda-forge r-base=4.4.1 -y
python3 BiomiX_INSTALL_LINUX.py
