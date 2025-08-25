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

   cd ecc_finderV1.1

1- create environmment
 
   conda create -n eccdna_master

2- Activate de Environment

   conda activate eccdna_master

3- upload the packages

   conda env update -f packages.yml –prune

4- config file config.yaml for configuration of parameters you wan
     

5- execute the script with the data in config.yaml (there is some default parameters in config.yaml)

   python exe.py 

6- Parameters in config.yaml that you can change before launch eccfinder

Function|	Parameters_eccfinder|	value|	description|
-----------------------------------------------------------------
Options for peak-calling:GenRich	-yv		
	-l|	200|	Minpeak Minimum length of a peak|
	-g|	100|	Maxdist Maximum distance between signif. sites (def. 100)|
	-p|	0.05|	Maxpvalue Maximum p-value (def. 0.01)|    
Not Using in eccfinder|	 -a| 	 200.0|	Minimum AUC for a peak (def. 200.0)|
TideHunter Tandem Criteria|	 -n| 	2|	minimum copy number of tandem repeat in a long read (>=2) [2]|
	-e|	0.25|	maximum allowed divergence rate between two consecutive repeats| 
	 -s	30	minimum period size of tandem repeat (>=2) [30]
Not Using in eccfinder|	-M| 		 match score [2]|
Not Using in eccfinder| 	-X|		mismatch penalty [4]|
Not Using in eccfinder|	-O| 		gap opening penalty (O1,O2) [4,24]|
FilteringBedFunction|	--min-read|	3|	filter locus by unique mapped read number|
	--min-bound|	0.8|	filter locus at regions by boundary coverage|
	--min-cov|	10|	minimum coverage of detected eccDNA loci|
