#!/bin/bash

#source ~/miniconda3/etc/profile.d/conda.sh
SCRIPT_DIR=$(dirname "$0")
PYTHON_SCRIPT="BiomiX_INSTALL_OS.py"
conda create -n BiomiX-env python=3.9 -y
conda activate BiomiX-env
#/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
conda config --add channels conda-forge
#conda install --solver=classic conda-forge::conda-libmamba-solver conda-forge::libmamba conda-#forge::libmambapy conda-forge::libarchive -y

conda install -c conda-forge r-base=4.4.1 -y
conda install conda-forge::glib=2.82.1 -y
brew install libgit2
#brew install freetype
#brew install harfbuzz
#brew install fribidi
#brew install libpng
 

PYTHON_SCRIPT_PATH="$SCRIPT_DIR/$PYTHON_SCRIPT"
echo $PYTHON_SCRIPT_PATH
python3 $PYTHON_SCRIPT_PATH
