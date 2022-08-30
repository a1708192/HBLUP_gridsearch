# HBLUP_gridsearch
The performance of HBLUP using various hyper-parameters such as blending, tuning and scale factor in simulated data.

The H-matrix best linear unbiased prediction (HBLUP) method has been widely used in livestock breeding programs. The existing HBLUP method (e.g., that implemented in BLUPf90 software) requires hyper-parameters that should be adequately optimised as otherwise the genomic prediction accuracy may decrease. In this study, we assess the performance of HBLUP using various hyper-parameters such as blending, tuning and scale factor in simulated data.

The code has been developed on a LINUX server and a bash script is provided for running the code on a server (test_sim_gadi.sh)

#!/bin/bash
 
#PBS -l ncpus=4   
#PBS -l mem=6GB  
#PBS -l jobfs=7GB  
#PBS -P eu82  
#PBS -l walltime=48:00:00  
#PBS -M name.family@uni.edu.au  
#PBS -l wd  
#PBS -m abe  
#PBS -q normal  
module load R  
R CMD BATCH --no-save run_all.R  
#----------------------------------------------------  
**Software requirments:**


The main R script code is entitled 'run_all.R'. In this file,  
