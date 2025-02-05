import os
import subprocess

# List of commands to execute
commands = [
    'conda install -c conda-forge pkg-config=0.29.2 -y',
    'conda install -c anaconda libcurl -y',
    'conda install -c conda-forge libnetcdf -y',
    'conda install -c anaconda libxcb -y',
#   'conda install -c conda-forge pyqt -y',
    'pip install PyQT5==5.12.3',
    'conda install -c anaconda pandas -y',
    'pip install scikit-learn',
    'pip install xlrd',
    'pip install openpyxl',
    'pip install mofapy2==0.7.1',
    'conda install -c conda-forge r-ncdf4 -y'
]

# Execute commands
for command in commands:
    subprocess.run(command, shell=True)

# Get the path of the current Python script
python_script_path = os.path.realpath(__file__)

# Construct the path to the R script in the same directory as the Python script
r_script_path = os.path.join(os.path.dirname(python_script_path), "INSTALL_BiomiX_OS.r")

print("INSTALL R PACKAGE")
# Run the R script using the subprocess module
subprocess.run(["Rscript", r_script_path])

#print("VERIFY R PACKAGE MISSED")
# Run the R script using the subprocess module
#subprocess.run(["Rscript", r_script_path])
