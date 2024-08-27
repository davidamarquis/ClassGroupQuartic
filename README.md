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

## Configuration
### Top level params
smoothness bound,
sieve radius,
number of primes,

Some discussion of some of the tuning params can be found in 
/Users/davidmarquis/git/thesis/src/data_workflow/sane_parameter_defaults.py

### Parameters for solving a specific instance
Parameters for solving a specific instance are stored in file with names ending in "_params"
The naming convention for the prefix after "task" is that if it's a constant used for a single stage in [init, sieve, filter, linalg] then use that
and constants that apply to more than one task proceed directly to the name.
Examples,
```
tasks.init.max_num_primes = 
tasks.sieve.mfb0 = 
```
### Init task
init task includes 2 steps: selecting a least_integer and initializing the sieve

### Constants that apply at the server level
Constants that apply at the server level are stored in ```constants.py```

## Data Workflow
Data is generated and processed by the scripts in ```src.data_workflow```

We start by constructing polys.
```src.data_workflow.def_polys``` currently there are files for making a specific urank OR using a construction.
These are stored in src.polredabs_polys.

There are 3 stages of data: Raw, derived and tidy. 
Raw data==a log dump from the algorithm
Derived data==a csv of data scraped from the raw data.
Tidy data==a csv of a subset of data from Derived that I can use in tables.

Unfortunately, because there is a server involved I need 2 steps for generating the raw data.
So we have 4 steps:
S1 Puts files on the server so that a job can be started. Local. Modules: magma_init, SI_init. 
S2 Runs job and outputs raw data. Remote. Text. modules: magma_raw, SI_raw. 
S3 Downloads data and outputs derived data. Local. Csv
S4 Run locally. 

### Usage 
magma raw: Runs on remote and takes a long time. Think carefully before firing off this module. 
S3 can be run often to save partial work done on the server.

magma init runs locally and runs fairly fast so editing file and rerunning is fine.

## Data layout
Data for my thesis is stored in a separate repo.
The code for the data workflow assumes a certain directory structure so here's a description:
```
data/ for raw data downloaded from server
data_dervied/ for derived
data_tidy/ for tidy
```

path constants need to be set 
basepath_save points to data/
basepath_save_remote points to data/ on the server
basepath_artifacts A server path that says where the intermediate stuff goes. This could be under /tmp

Experimental data is dumped in /Users/davidmarquis/git/thesis/data_for_exprs 

### data/
Results about an individual nf, like a Magma log, should be in nf_label/...
Previously, I had been collecting data in separate dirs for Magma and SI so there's old junk dirs from that.

### Data Structures 
Relations are stored as rows of the matrix
Its more natural to look at matrices with more rows than columns than the other way around: 
-on a printed page there's more vertical space than horizontal space
-on a screen scrolling vertically is easiest
