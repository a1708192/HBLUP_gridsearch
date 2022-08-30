
Make_Ped_file_PLINK=function(Num_Chr,nmloci,nmSnps,QMSim_filename)
{
#  R CMD BATCH --no-save convert_QMSim_PLINK.r
#install.packages(c( "foreach", "doParallel") ) 
# clear workspace
#rm(list=ls())
id=1
    print("------------------ Making Ped file is running --------------")

#set working directory
#setwd("C:/Users/neshatm/R_code_plink/QMSim_Linux/QMSim_process_test/")
  
pathfile= getwd()
pathfile

data_filename= paste0(pathfile,"/r_",QMSim_filename,"/p1_data_00",id,".txt")
data_filename

# read necessary QMSim files
P_data <- read.table(data_filename, header = FALSE, sep = "", dec = ".")

#data_filename= paste0("p1_mrk_00",id,".txt")
data_filename= paste0(pathfile,"/r_",QMSim_filename,"/p1_mrk_00",id,".txt")
data_filename

P_mrk <- read.delim(data_filename, header = TRUE, sep = "\t", dec = ".")

Size_P_mrk = nrow(P_mrk)
print(paste(Size_P_mrk, " is number of P_mrk row ")) 


Dim = dim(P_data)
print(paste(Dim[1], " is number of P_data row:: ", Dim[2], " is number of P_data column:: " ))

# Make the .ped file dataframe
FID <- c(1:(Dim[1]-1))
IID <- c(1:(Dim[1]-1))

Sex_temp <- P_data[[4]]
Sex= 0
for(i in 2:length(Sex_temp)){
   if (Sex_temp[i]=="M") {  Sex[i-1]= 1}
     else {Sex[i-1]=2}
}

rm(Sex_temp)

PID_TEMP <- P_data[[2]]
PID_TEMP=PID_TEMP[-c(1)]
PID = as.numeric(as.character(PID_TEMP))
#PID = rep(0, Size_P_mrk) 
rm(PID_TEMP)

MID_TEMP <- P_data[[3]]
MID_TEMP=MID_TEMP[-c(1)]
MID = as.numeric(as.character(MID_TEMP))
#MID = rep(0, Size_P_mrk) 
rm(MID_TEMP)

P <- rep(-9,times=(Dim[1]-1))
Len_snps= nmSnps*2 #nchar(row1)
mat_snps = matrix(0,Size_P_mrk ,Len_snps )

for (ii in 1:Size_P_mrk) { 
 row1=P_mrk[ii,1]
 te=strsplit(as.character(row1), "\\s+")[[1]]
 rm(row1)
 row1=te[2]
 k=1
 
  for (j in 1:nmSnps){
   Ch=substr(row1, j, j)
   if (Ch!=" "){
     if (Ch==0){
        mat_snps[ii,k]=1
        k=k+1
        mat_snps[ii,k]=1
        k=k+1}
     else if (Ch==2){
        mat_snps[ii,k]=2
        k=k+1
        mat_snps[ii,k]=2
        k=k+1}
     else if (Ch==3){
        mat_snps[ii,k]=1
        k=k+1
        mat_snps[ii,k]=2
        k=k+1}
     else if (Ch==4){
        mat_snps[ii,k]=2
        k=k+1
        mat_snps[ii,k]=1
        k=k+1}
     else if (Ch==5){
        mat_snps[ii,k]='NA' 
        k=k+1
        mat_snps[ii,k]='NA' 
        k=k+1}
     
   }
     
} # end for
} # end for
Phen=P_data$V10
Phen=Phen[-1]
Phen=as.numeric(as.character(Phen))
filename=paste0("Snps",id,".phen")
 snps_phen <- data.frame("FID" = FID, "IID" = IID,"PHEN"=Phen ,stringsAsFactors = FALSE)
 write.table(snps_phen, file=filename, append = FALSE, sep = " ", dec = ".",
             row.names = FALSE, col.names = FALSE)
 
 filename=paste0("Snps",id,".ped")
 snps_ped <- data.frame("FID" = FID, "IID" = IID, "PID" = PID, "MID" = MID, "Sex" = Sex, "P" = P, "SNPs" =mat_snps, stringsAsFactors = FALSE)
  #a good reference: https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/write.table
 write.table(snps_ped, file=filename, append = FALSE, sep = " ", dec = ".",
             row.names = FALSE, col.names = FALSE)
#temp1=snps_ped[3551:5550,]
# filename=paste0("Snps_last2",id,".ped")
#write.table(temp1, file=filename, append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)

}