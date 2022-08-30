prediction_accuracy_2=function(id,op_iter,tune,blend,alpha_value,target_index)
{
  
 
  
print(paste('Tune=', tune, ':: blend=', blend, ':: alpha value=', alpha_value))
test_data_nan=read.table("snps_hwe_001_ref.phen")
predict_phen=read.table("snps_hwe_001.bv", skip=1)
test_data =read.table("snps_hwe_001.phen")
estimate_h2=read.table("snp_GREML", nrow=3)

index_nan=target_index #which(is.na(test_data_nan$V3))
true_phen= test_data$V3[index_nan]
index_m= match(test_data$V1[index_nan],test_data_nan$V1)

predictions = predict_phen$V2[index_m]
R_value = cor(predictions, true_phen)
R2_value = R_value^2
RMSE_ALL = sqrt(mean((true_phen - predictions)^2))
MAE_ALL  = mean(abs(true_phen - predictions))

file_name=paste0("accuracy_prediction_tune_",tune,"_blend_",blend,"alpha_value_",alpha_value,"_sim_grm.txt")

if (id==1) {
  accuracy_pre= data.frame( R = R_value, R2= R2_value,
                            RMSE_value = RMSE_ALL,
                            MAE_value = MAE_ALL ,
                            h2 = estimate_h2$V2[3],
                            h2_se = estimate_h2$V3[3],
                            Va = estimate_h2$V2[2],
                            Va_se = estimate_h2$V3[2],
                            Ve = estimate_h2$V2[1],
                            Ve_se = estimate_h2$V3[1]) }

else{
  accuracy_pre =read.table(file_name, header=F)
  temp=c(R_value,R2_value,
         RMSE_ALL,
         MAE_ALL ,
         h2 = estimate_h2$V2[3],
         h2_se = estimate_h2$V3[3],
         Va = estimate_h2$V2[2],
         Va_se = estimate_h2$V3[2],
         Ve = estimate_h2$V2[1],
         Ve_se = estimate_h2$V3[1])
  
  accuracy_pre = rbind(accuracy_pre,temp)
}


if (id %% 1 == 0) {
write.table(accuracy_pre, file=file_name, append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE, quote=FALSE)

}

if (id==5){
  test_data =read.table(file_name)
  R = mean(test_data$V1)
  R2= mean(test_data$V2)
  RMSE_value = mean(test_data$V3)
  MAE_value = mean(test_data$V4) 
  h2 = mean(test_data$V5)
  h2_se = mean(test_data$V6)
  Va = mean(test_data$V7)
  Va_se = mean(test_data$V8)
  Ve = mean(test_data$V9)
  Ve_se = mean(test_data$V10)
  
  
  
  if(op_iter==1){
    op_accuracy=data.frame(Tune=tune, Blend=blend, Alpha_value=alpha_value,R_all=R,R2_all=R2,
                    RMSE_ALL=RMSE_value,
                    MAE_ALL=MAE_value,
                    H2_ALL=h2,
                    H2_SE_ALL=h2_se,
                    VA_ALL=Va,
                    VA_SE_ALL=Va_se,
                    VE_ALL=Ve,
                    VE_SE_ALL=Ve_se )
    
    write.table(op_accuracy, file="random_search_op_hyper.txt", append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE, quote=FALSE)

  print("5-fild cross valiadation file  is removing...............")
 # unlink(file_name, recursive = TRUE)

  } else{
    op_accuracy =read.table("random_search_op_hyper.txt")
    temp=c(tune,blend,alpha_value,R,R2,RMSE_value,MAE_value,h2,h2_se,Va,Va_se,Ve,Ve_se )
    op_accuracy = rbind(op_accuracy,temp)
    write.table(op_accuracy, file="random_search_op_hyper.txt", append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE, quote=FALSE)
    
    print("5-fild cross valiadation file  is removing...............")
    unlink(file_name, recursive = TRUE)

  }
  
}

#system("awk 'n!=1 && $1!='LKH' {printf $2' '}; $1=='LKH' {n=1;printf $2'\n'}' snp_GREML >> GREML_h2.out")
#system("awk 'n!=1 && $1!='LKH'' {printf $3' ''}; $1=='LKH' {n=1;printf $2'\n'}' snp_GREML >> GREML_SE.out")


}
