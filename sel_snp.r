sel_snp=function(){
# clear workspace
#rm(list=ls())

#set working directory
#setwd("C:/Users/neshatm/R_code_plink/plink_win64_20210416/")

  bim=read.table("snps_hwe_001.freq")
  v1=seq(1:length(bim$V2))
  v2=sample(v1,1000)

  v12=rnorm(length(v2),0,1)
  v12=v12*sqrt(1/(bim$V3[v2]*(1-bim$V3[v2]/2))^0)

  v3=cbind(v2,as.character(bim$V2[v2]),v12,bim$V3[v2])
  sink("snp.lst")
  write.table(v3,quote=F,col.name=F,row.name=F)
  sink()

}