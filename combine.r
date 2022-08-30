
combine=function(){
#install.packages("lava")
# clear workspace
#rm(list=ls())

#set working directory
 # setwd("C:/Users/neshatm/R_code_plink/QMSim_Linux/QMSim_process_test/Estimate_Phen_MTG2/")
  
 
fam=read.table("snps_hwe_001.fam")
yv=read.table("snps_hwe_001.bv")

# residual variance
rv = 0.5
  v2=rnorm(length(yv$V1),0,rv)  # (if rv=1, h2=0.5), (if rv=2, h2=0.2), (if rv=0.5, h2=0.8) 

  out_y=yv$V1+v2
  v1=cbind(as.character(fam$V1),as.character(fam$V2),out_y)

  sink("snps_hwe_001.phen")
  write.table (v1,row.names=F,col.names=F,quote=F)
  sink()
  
  jpeg(file="hist_genetic_value.jpeg")
  hist(yv[,1],breaks=10,main="Histogram of genetic values mu=0 and sigma=1")
  dev.off()

  jpeg(file="hist_residual_value.png")
  hist(v2,breaks=10,main="Histogram of residual values mu=0 and sigma=1")
  dev.off()
 
  jpeg(file="hist_phenotype_value.png")
  hist(out_y,breaks=10,main="Histogram of phenotype ")
  dev.off()

  cat("phenotypic variance           :",var(out_y),"\n")
}