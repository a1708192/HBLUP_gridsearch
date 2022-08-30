

Make_Map_file_PLINK=function(Num_Chr,nmloci,nmSnps)
{
# Num_Chr=30, nmloci=500, nmSnps=15000, 
#set working directory
#setwd("C:/Users/neshatm/R_code_plink/QMSim_Linux/QMSim_process_test/")
  print("------------------ Making Map file is running --------------")
  # print the genetic variables
  print(paste(Num_Chr,"------Number of Chromosome----"))
  print(paste(nmloci,"------Number of loci----"))
  print(paste(nmSnps,"------Number of SNPs----"))
  
Chr= rep(c(1:Num_Chr), each = nmloci)

SN = 1:nmSnps

GD= rep(0,times=nmSnps)

Bp= seq(100000, nmloci*100000, by=100000) 

BPP= rep(Bp,times=Num_Chr)

snps_map <- data.frame("Chr" = Chr, "SN" = SN, "GD" = GD, "BPP" = BPP,stringsAsFactors = FALSE)

filename = paste0("Snps",1,".map")
#write.csv(snps_map, "snps_map.csv", row.names=FALSE, quote=FALSE) 
write.table(snps_map, file=filename, append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE, quote=FALSE)

}
#"snps_map_T.txt"