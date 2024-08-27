### Required Dependencies
Magma   
Cado NFS  
Sage   

Sage is used as a Python library.
The following Python packages are not included in Sage and must be install: pytest, tqdm.
This can be done with pip using ```python -m pip install PACKAGE_NAME```

### Setup
- Install the required dependencies
- Open paths_temp.py and constants_temp.py and follow the instructions in them

### Parameters and defining polynomial files 
Like CADO-NFS, two kinds of configuration files are used.
Parameters for solving a specific instance are stored in file with suffix "_params".
Defining polynomial coefficients are stored in a file with suffix ".poly"
