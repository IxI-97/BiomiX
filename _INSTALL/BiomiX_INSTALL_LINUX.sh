#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda create -n BiomiX-env python=3.9 r-base=4.2.0 r-essentials
conda activate BiomiX-env
python3 BiomiX_INSTALL_LINUX.py
