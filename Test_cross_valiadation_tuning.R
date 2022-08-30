
Test_cross_valiadation_tuning=function(job_iteration,Tune, Blend, alpha_value) {
  
# HBLUP

print (paste("============= iter=",job_iteration))

#load the function
source("Make_Map_file_PLINK.R")
source("Make_Ped_file_PLINK.R")
source("sel_snp.r")
source("combine.r")
source("Estimate_H_Matrix.R")
source("prediction_accuracy_2.R")
  
Num_Chr=30
nmloci=300
nmSnps=9000
QMSim_filename="Hist_pop_QMSim"

 if (job_iteration==1) {

  #call the QMSim software
  
  system("./QMSim Hist_pop_QMSim.prm")
  
  # Make map plink file
  
  Make_Map_file_PLINK(Num_Chr,nmloci,nmSnps) 
  
  
  
  # make ped plink file
  
  Make_Ped_file_PLINK(Num_Chr,nmloci,nmSnps,QMSim_filename)
  
  
  print("The QMSim folder is removing...............")
  unlink("r_Hist_pop_QMSim", recursive = TRUE) # will delete directory of QMSim outputs
  
  
  print("QC process")
  system("./plink --noweb --file Snps1 --make-bed --out snp_binary --chr-set 30")
  system("./plink --bfile snp_binary --maf 0.01 --make-bed --out snps_maf_01 --chr-set 30 --nonfounders")
  system("./plink --bfile snps_maf_01 --hwe 0.001 --make-bed --out snps_hwe_001 --chr-set 30")
  
  
  
  print ("Extract the target file from the whole population based on the last two generations")
    
  # make the ID file for the subpopulation (FID, IID)
  fam_last2gen_temp=read.table("snps_hwe_001.fam",skip=3550)
  fam_last2gen=data.frame(V1 = fam_last2gen_temp[, 'V1'], V2 = fam_last2gen_temp[, 'V2'])
  write.table(fam_last2gen, file="fam_last2gen.txt", append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE, quote=FALSE)

   # Extract the plink file from the original population file after QC
  system("./plink --bfile snps_hwe_001 --keep fam_last2gen.txt --make-bed --out last2gen_binary --chr-set 30")

  # Make fam file of target and QC process
  print("Make fam file of the target population")
  #system("./plink --noweb --file last2gen --make-bed --out last2gen_binary --chr-set 30")
  system("./plink --bfile last2gen_binary --maf 0.01 --make-bed --out last2gen_maf_01 --chr-set 30 --nonfounders")
  system("./plink --bfile last2gen_maf_01 --hwe 0.001 --make-bed --out last2gen_hwe_001 --chr-set 30")
  
  
  print(" make population_information.pop for different scale factor")
  temp=read.table("snps_hwe_001.fam")
  pop = matrix(1,nrow(temp),1)
  pop_file = data.frame(temp$V1,temp$V2,pop)
  write.table(pop_file, file="population_information.pop", append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE, quote=FALSE)

  } # end if for runing QMSim
 #-----------------------------------------------------------
  
# simulate phenotypes based on the given genotyped data

 if (job_iteration==1) {
  
  system("./mtg2 -plink snps_hwe_001 -frq 1")
  
  # Obtain snps list from allele frequency
  sel_snp()
  
  #Generate .bv file (genetic value with mean 0 and variance 1)
  system("./mtg2 -plink snps_hwe_001 -simreal snp.lst")
  
  #simulated phenotype will be generated using the residual values
  combine()
  }

#----------------------------------------
     
  # Estimate H matrix 
  
 

 
  if (alpha_value == -1.5) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_1 -out GRM_last2gen")
}
  
  if (alpha_value == -1.4) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_2 -out GRM_last2gen")
}

if (alpha_value == -1.3) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_3 -out GRM_last2gen")
}

if (alpha_value == -1.2) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_4 -out GRM_last2gen")
}

if (alpha_value == -1.1) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_5 -out GRM_last2gen")
}


if (alpha_value == -1) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_6 -out GRM_last2gen")
}

if (alpha_value == -0.9) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_7 -out GRM_last2gen")
}

if (alpha_value == -0.8) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_8 -out GRM_last2gen")
}

if (alpha_value == -0.7) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_9 -out GRM_last2gen")
}

if (alpha_value == -0.6) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_10 -out GRM_last2gen")
}

if (alpha_value == -0.5) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_11 -out GRM_last2gen")
}

if (alpha_value == -0.4) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_12 -out GRM_last2gen")
}

if (alpha_value == -0.3) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_13 -out GRM_last2gen")
}

if (alpha_value == -0.2) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_14 -out GRM_last2gen")
}

if (alpha_value == -0.1) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_15 -out GRM_last2gen")
}

if (alpha_value == 0) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_16 -out GRM_last2gen")
}

if (alpha_value == 0.1) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_17 -out GRM_last2gen")
}

if (alpha_value == 0.2) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_18 -out GRM_last2gen")
}

if (alpha_value == 0.3) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_19 -out GRM_last2gen")
}

if (alpha_value == 0.4) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_20 -out GRM_last2gen")
}

if (alpha_value == 0.5) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_21 -out GRM_last2gen")
}

if (alpha_value == 0.6) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_22 -out GRM_last2gen")
}

if (alpha_value == 0.7) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_23 -out GRM_last2gen")
}

if (alpha_value == 0.8) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_24 -out GRM_last2gen")
}

if (alpha_value == 0.9) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_25 -out GRM_last2gen")
}

if (alpha_value == 1) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_26 -out GRM_last2gen")
}

if (alpha_value == 1.1) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_27 -out GRM_last2gen")
}

if (alpha_value == 1.2) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_28 -out GRM_last2gen")
}

if (alpha_value == 1.3) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_29 -out GRM_last2gen")
}

if (alpha_value == 1.4) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_30 -out GRM_last2gen")
}

if (alpha_value == 1.5) {
  print(paste(" Estimating Genomic relationship matrix for different scale factors:: alpha value is ", alpha_value))
  system("./mtg2 -plink last2gen_hwe_001 -rtmx2 rtmx_parameters_31 -out GRM_last2gen")
}


#--------------------------------------
  print(paste("Tune=",Tune,"::  Blend=",Blend,":: alpha=",alpha_value,"::job_iteration=",job_iteration))

  Estimate_H_Matrix(job_iteration,Tune, Blend,alpha_value)

}