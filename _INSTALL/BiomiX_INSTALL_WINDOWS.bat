
call %HOMEDRIVE%%HOMEPATH%\miniconda3\Scripts\activate.bat %HOMEDRIVE%%HOMEPATH%\miniconda3
call conda create -n BiomiX-env python=3.9
call conda activate BiomiX-env
call conda install -c conda-forge r-base=4.1.3
call conda install -c conda-forge r-igraph
call python BiomiX_INSTALL_WINDOWS.py

pause

