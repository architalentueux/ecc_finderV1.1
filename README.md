# ecc_finderV1.1
Tool for eccdna analysis using ONT only 
![Optional Text](../main/data/workflow_eccfinder1.png)
## I- Install Miniconda3 First


  - mkdir -p ~/miniconda3
  - wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
  - bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
  - rm ~/miniconda3/miniconda.sh

### After instalation
   source ~/miniconda3/bin/activate

### Config chanel Conda for Bioconda

   - conda config --add channels defaults
   - conda config --add channels bioconda
   - conda config --add channels conda-forge
   - conda config --set channel_priority strict

## II- Process To launch Eccfinder on ont map-ont

### 0- download the tool

   - git clone https://github.com/architalentueux/ecc_finderV1.1.git

### 1- create environmment
 
   - conda create -n eccdna_master

### 2- Activate the Environment

   - conda activate eccdna_master

### 3- upload the packages

   - conda env update -f packages.yml –prune

### 4- config file config.yaml

   - configuration with parameters you want in config.yaml

### 5- execute the script 

   - python exe.py 
