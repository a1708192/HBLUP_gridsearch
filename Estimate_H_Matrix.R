Estimate_H_Matrix=function(op_iter,Tune,Blend,alpha_value)
{
  
  print(paste("Tune=",Tune,"::  Blend=",Blend,":: alpha=",alpha_value,"::job_iteration=",op_iter))
  #make_command = paste("./mtg2 -p last2gen_hwe_001.fam -zg GRM_last2gen.grm.gz -hrtmx snps_hwe_001.fam -tune", Tune, "-blend", Blend," -thread 10  -out last2gen_H.hrtmx -grmout Modified_GRM")
  make_command = paste("./mtg2 -p last2gen_hwe_001.fam -g GRM_last2gen -hrtmx snps_hwe_001.fam -tune", Tune, "-blend", Blend," -thread 20  -out last2gen_H.hrtmx")


  # Estimate H-matrix using the tune and blend
  system(make_command)
  
  
  

  #  Validation
  Index_sampling=read.table('Index_sampling.txt') 
  head(Index_sampling)
  N_phen = nrow(Index_sampling)

  fold=round(N_phen/5 )-1

  print(paste0('iteration= ',op_iter,' ================================================'))
  print("five-fold cross validation")


  pheno=read.table("snps_hwe_001.phen")
 #------------ five-fold cross validation --------
for (r in 1:5){
  print (paste('iteration of 5-fold cross validation=',r,' :: evaluation number = ',op_iter ))
  s_index = (r-1)*fold+1
  e_index = s_index-1+fold
  print(paste('row number=',nrow(Index_sampling),':: column number=',ncol(Index_sampling),':: start=',s_index,':: end=',e_index))
  target_index  =  Index_sampling[s_index:e_index,1]  #sample(3550:5550, 1000, replace=F) # 1000 random samples from the last two generations
  Whole_pop_index=1:5550
  reference_index= Whole_pop_index[-target_index]
  
  
  for (j in 1:length( target_index)){
    pheno$V3[target_index[j]]='NA'
  }
  
  i=1
  file_name= paste0("snps_hwe_00",i,"_ref.phen")
  write.table(pheno, file=file_name, append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE, quote=FALSE)
  
  
  print("HGREML: Estimating the heritability based on the phenotypic correlation and sample correlation")
  REF_NAME = paste0('snps_hwe_00',i,'_ref.phen')
  BV_NAME = paste0('snps_hwe_00',i,'.bv')
  MTG_command = paste('./mtg2 -p snps_hwe_001.fam  -d ',REF_NAME,' -g last2gen_H.hrtmx thread 20 -out snp_GREML -bv ',BV_NAME)
  # system("./mtg2 -p snps_hwe_001.fam  -d snps_hwe_001_ref.phen  -zg GRM_snps.grm.gz thread 20 -out snp_GREML -bv snps_hwe_001.bv")
  system (MTG_command)
  
  

   prediction_accuracy_2 (r,op_iter,Tune,Blend,alpha_value,target_index)  

} # end 5-fold cross

} # end function
