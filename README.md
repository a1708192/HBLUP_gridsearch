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
**Software requirements:**

**MTG2** : https://sites.google.com/site/honglee0707/mtg2  
**PLINK**: https://www.cog-genomics.org/plink/  
**QMSim**: https://animalbiosciences.uoguelph.ca/~msargol/qmsim/  

#---------------------------------------------------  

**The main R script:**

In this script entitled 'run_all.R', the main core of the grid search is developed to evaluate all feasible configurations of the HBLUP hyper-parameters. It is noted that if the resolution of the alpha parameter changes, the relevant rtmx_parameters script should be developed.   

#----------------------------------------------------  
**QMSim profile:**  

The technical details and control parameters of the simulations of the historical populations can be seen in '*.prm' script.    
/*******************************  
  **     Global parameters     **  
  *******************************/  
  title = "Simulation historical population";  
nrep  = 1;                    //Number of replicates  
h2    = 0.8;                  //Heritability  
qtlh2 = 0.0;                  //QTL heritability 
phvar = 1.0;                  //Phenotypic variance  
  
  /*******************************  
    **   Historical population   **  
    *******************************/  
    begin_hp;  
    hg_size = 100 [0]    
               100 [95]  
               1000 [100];      //Size of the historical generations (with effective population size Ne)  
   nmlhg   = 50;              //Number of males in the last generation  
   nmfhg   = 50;              // Number of males in the first generation   
  end_hp;  
  
  /*******************************  
    **        Populations        **  
    *******************************/  
    begin_pop = "p1";  
     begin_founder;  
       male [n = 50, pop = "hp", select = rnd];  
       female [n = 500, pop = "hp", select= rnd];  
     end_founder;  
     ng = 5;  
     ls = 2;          //Litter size  
     pmp = 0.5 /fix; //Proportion of male progeny  
     md = rnd;       //Mating design (control of inbreeding) ;rnd, rnd_ug, minf, maxf, p_assort, n_assort  
     sd = rnd;       //Selection design; rnd, phen, tbv and ebv  
     ebv_est = blup ; // Breeding value estimation method; The best linear unbiased prediction (BLUP) of breeding values  
     
      begin_popoutput;  
       data;  
       stat;   //Save brief statistics on simulated data  
       genotype /snp_code ;  
       allele_freq;   //Save allele frequencies 
       ld;          //Save linkage disequilibrium stat  
     end_popoutput;  
     
    end_pop;  
  
  /*******************************  
    **          Genome           **  
    *******************************/  
    begin_genome;  
    
      begin_chr = 30;  
         chrlen = 100;     //Chromosome length (cM)  
         nmloci = 300;     //Number of markers 9000/30  
         mpos   = rnd;     //Marker positions  
         nma    = all 2;    //Number of marker alleles  
         maf    = eql;     //Marker allele frequencies  
         nqloci = 90;         //Number of QTL 90 990 9990 60000  
         qpos   = rnd;       //QTL positions  
         nqa    = all 2;     //Number of QTL alleles  
         qaf    = eql;       //QTL allele frequencies  
         qae    = rndg 0.4;  //QTL allele effects % effects are sampled from gamma distribution  
     end_chr;
      mmutr  = 2.5e-8 /recurrent;      // Marker mutation rate  
      qmutr  = 2.5e-8;     //QTL mutation rate  
      r_mpos_g; //Randomize marker positions across genome  
      r_qpos_g; //Randomize QTL positions across genome  
  end_genome;  
  
  /*******************************  
    **       Output options      **  
    *******************************/  
    begin_output;  
     linkage_map; // Marker and QTL linkage map (GWAS)  
     allele_effect;  
    end_output;  
#--------------------------------------------------------------  

