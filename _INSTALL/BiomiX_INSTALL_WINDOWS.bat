
call %HOMEDRIVE%%HOMEPATH%\miniconda3\Scripts\activate.bat %HOMEDRIVE%%HOMEPATH%\miniconda3
call conda create -n BiomiX-env python=3.9 -y
call conda activate BiomiX-env
call conda install mamba -y
call mamba install -c conda-forge r-base=4.4.1 -y
call mamba install -c conda-forge r-igraph -y
call python BiomiX_INSTALL_WINDOWS.py

pause
