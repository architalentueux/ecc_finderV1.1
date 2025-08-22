#I- This First Step it's if you doesn't have Miniconda3 you can skip that


mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh

#After instalation
source ~/miniconda3/bin/activate

#Config chanel Conda for Bioconda

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# II- Process To launch Eccfinder on ont map-ont
0- download the software from the repository

   git clone https://github.com/architalentueux/ecc_finderV1.1.git

1- create environmment
 
   conda create -n eccdna_master

2- Activate de Environment

   conda activate eccdna_master

3- upload the packages

   conda env update -f packages.yml –prune

4- config file config.yaml for configuration of parameters you wan
     

5- execute the script with the data in config.yaml (there is some default parameters in config.yaml)

   python exe.py 
