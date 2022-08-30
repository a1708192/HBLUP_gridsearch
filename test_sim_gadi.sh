#!/bin/bash
 
#PBS -l ncpus=4
#PBS -l mem=6GB
#PBS -l jobfs=7GB
#PBS -P eu82
#PBS -l walltime=48:00:00
#PBS -M Mehdi.Neshat@unisa.edu.au
#PBS -l wd
#PBS -m abe
#PBS -q normal
module load R
R CMD BATCH --no-save run_all.R



  