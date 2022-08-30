print("***************************************")
print(" Grid search for tuning the hyper-parameters of HBLUP for simulated data") 
# Developer : Dr Mehdi Neshat (neshat.mehdi@gmail.com)
# Supervisor : Prof Hong Lee
# University of South Australia
# Date : 30/ 07/ 2022
#-----------------------------------------

#load the function
source("Test_cross_valiadation_tuning.R")
s=1
alpha       = seq(-1.5,1.5,0.1)
blend       = seq(0,1,0.1)
tune        = seq(0,2,1)
 

  for (i in 1:3) { # tuning
 
  for(j in 1:length(blend)){
    
    for(k in 1:length(alpha)){

 Tune_t  = round(tune[i],digits=0)
 Blend_t = round(blend[j],digits=2)
 Alpha_t = round(alpha[k],digits=1)
 op_iter = s
 print(paste('Tune =', Tune_t,'::: Blend=', Blend_t, ':::Alpha_t=',Alpha_t))
  start_time <- Sys.time()
 
  Test_cross_valiadation_tuning(op_iter,Tune_t, Blend_t, Alpha_t)
  end_time <- Sys.time()
  print(paste("runtime of this replication is around=",end_time - start_time))
  
s=s+1
  }
     } 
         }


