import os
import subprocess

# List of commands to execute
commands = [
    'mamba install -c conda-forge pkg-config=0.29.2 -y',
    'mamba install -c anaconda libcurl=8.11.1 -y',
    'mamba install -c conda-forge libnetcdf=4.9.2 -y',
    'mamba install -c anaconda libxcb=1.17.0 -y',
#   'conda install -c conda-forge pyqt -y',
    'pip install PyQT5==5.12.3',
    'mamba install -c anaconda pandas=2.2.3 -y',
    'pip install scikit-learn==1.6.1',
    'pip install xlrd==2.0.1',
    'pip install openpyxl==3.1.5',
    'pip install mofapy2==0.7.1',
    'mamba install -c conda-forge r-ncdf4==1.23 -y'
]

# Execute commands
for command in commands:
    subprocess.run(command, shell=True)

# Get the path of the current Python script
python_script_path = os.path.realpath(__file__)

# Construct the path to the R script in the same directory as the Python script
r_script_path = os.path.join(os.path.dirname(python_script_path), "INSTALL_BiomiX_LINUX.r")

print("INSTALL R PACKAGE")
# Run the R script using the subprocess module
subprocess.run(["Rscript", r_script_path])

#print("VERIFY R PACKAGE MISSED")
# Run the R script using the subprocess module
#subprocess.run(["Rscript", r_script_path])
