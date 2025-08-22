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

1- create environmment
 
   conda create -n eccdna_master

2- Activate de Environment

   conda activate eccdna_master

3- upload the packages

   conda env update -f packages.yml –prune

4- config file config.yaml for configuration of parameters you wan
     

<<<<<<< HEAD
5- execute the script with the data in config.yaml (there is some default parameters in config.yaml)
=======
python exe_ecc.py map-ont /data/PAU29426_S13_sup_head30000.fastq -r /data/S.meliloti2011complete.fasta
>>>>>>> c1925922080ddff0f88064c0998162bf83a1c392

   python exe.py 
