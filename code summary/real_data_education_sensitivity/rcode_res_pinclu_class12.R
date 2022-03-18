load("data_class_as_clu_list.RData")
load("nmig12_class_posterior.RData") # class as cluster
load("nmig12_class_posterior_st2.RData")
load("nmig12_class_posterior_st3.RData")
load("nmig12_class_posterior_st4.RData")
load("nmig12_class_posterior_st5.RData")
load("nmig12_class_posterior_st6.RData")
load("nmig12_class_posterior_st7.RData")
load("nmig12_class_posterior_st8.RData")
load("nmig12_class_posterior_st9.RData")
load("nmig12_class_posterior_st10.RData")



library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

iter=15000
burn=5000
thin=10

###### here dataid=subject
get_pinclu_res = function(dataid,data_list,iter,burn){
  set.seed(dataid)
  data_obs = data_list[[dataid]]$data_obs
  N = dim(data_obs)[1]
  m = length(unique(data_obs$clu))
  n = as.vector(table(data_obs$clu))
####### recall result from first MCMC 
  ####### recall result from first MCMC 
  pstra11_all = rbind(nmig12_class_posterior[[dataid]]$prob_inc11[seq((burn+1),iter,thin),],
                  nmig12_class_posterior_st2[[dataid]]$prob_inc11[seq((burn+1),iter,thin),],
                  nmig12_class_posterior_st3[[dataid]]$prob_inc11[seq((burn+1),iter,thin),],
                  nmig12_class_posterior_st4[[dataid]]$prob_inc11[seq((burn+1),iter,thin),],
                  nmig12_class_posterior_st5[[dataid]]$prob_inc11[seq((burn+1),iter,thin),],
                  nmig12_class_posterior_st6[[dataid]]$prob_inc11[seq((burn+1),iter,thin),],
                  nmig12_class_posterior_st7[[dataid]]$prob_inc11[seq((burn+1),iter,thin),],
                  nmig12_class_posterior_st8[[dataid]]$prob_inc11[seq((burn+1),iter,thin),],
                  nmig12_class_posterior_st9[[dataid]]$prob_inc11[seq((burn+1),iter,thin),],
                  nmig12_class_posterior_st10[[dataid]]$prob_inc11[seq((burn+1),iter,thin),]) 
  pstra11 = apply(pstra11_all,2,mean)
  
  pstra10_all = rbind(nmig12_class_posterior[[dataid]]$prob_inc10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st2[[dataid]]$prob_inc10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st3[[dataid]]$prob_inc10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st4[[dataid]]$prob_inc10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st5[[dataid]]$prob_inc10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st6[[dataid]]$prob_inc10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st7[[dataid]]$prob_inc10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st8[[dataid]]$prob_inc10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st9[[dataid]]$prob_inc10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st10[[dataid]]$prob_inc10[seq((burn+1),iter,thin),]) 
  pstra10 = apply(pstra10_all,2,mean)

  pstra00_all = rbind(nmig12_class_posterior[[dataid]]$prob_inc00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st2[[dataid]]$prob_inc00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st3[[dataid]]$prob_inc00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st4[[dataid]]$prob_inc00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st5[[dataid]]$prob_inc00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st6[[dataid]]$prob_inc00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st7[[dataid]]$prob_inc00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st8[[dataid]]$prob_inc00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st9[[dataid]]$prob_inc00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st10[[dataid]]$prob_inc00[seq((burn+1),iter,thin),]) 
  pstra00 = apply(pstra00_all,2,mean) 

  pstra_y11_all = rbind(nmig12_class_posterior[[dataid]]$prob_inc_y11[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st2[[dataid]]$prob_inc_y11[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st3[[dataid]]$prob_inc_y11[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st4[[dataid]]$prob_inc_y11[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st5[[dataid]]$prob_inc_y11[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st6[[dataid]]$prob_inc_y11[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st7[[dataid]]$prob_inc_y11[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st8[[dataid]]$prob_inc_y11[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st9[[dataid]]$prob_inc_y11[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st10[[dataid]]$prob_inc_y11[seq((burn+1),iter,thin),]) 
  pstra_y11 = apply(pstra_y11_all,2,mean)
  
  pstra_y10_all = rbind(nmig12_class_posterior[[dataid]]$prob_inc_y10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st2[[dataid]]$prob_inc_y10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st3[[dataid]]$prob_inc_y10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st4[[dataid]]$prob_inc_y10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st5[[dataid]]$prob_inc_y10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st6[[dataid]]$prob_inc_y10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st7[[dataid]]$prob_inc_y10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st8[[dataid]]$prob_inc_y10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st9[[dataid]]$prob_inc_y10[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st10[[dataid]]$prob_inc_y10[seq((burn+1),iter,thin),]) 
  pstra_y10 = apply(pstra_y10_all,2,mean)
  
  pstra_y00_all = rbind(nmig12_class_posterior[[dataid]]$prob_inc_y00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st2[[dataid]]$prob_inc_y00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st3[[dataid]]$prob_inc_y00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st4[[dataid]]$prob_inc_y00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st5[[dataid]]$prob_inc_y00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st6[[dataid]]$prob_inc_y00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st7[[dataid]]$prob_inc_y00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st8[[dataid]]$prob_inc_y00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st9[[dataid]]$prob_inc_y00[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st10[[dataid]]$prob_inc_y00[seq((burn+1),iter,thin),]) 
  pstra_y00 = apply(pstra_y00_all,2,mean)
  
  pstra_y01_all = rbind(nmig12_class_posterior[[dataid]]$prob_inc_y01[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st2[[dataid]]$prob_inc_y01[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st3[[dataid]]$prob_inc_y01[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st4[[dataid]]$prob_inc_y01[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st5[[dataid]]$prob_inc_y01[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st6[[dataid]]$prob_inc_y01[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st7[[dataid]]$prob_inc_y01[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st8[[dataid]]$prob_inc_y01[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st9[[dataid]]$prob_inc_y01[seq((burn+1),iter,thin),],
                nmig12_class_posterior_st10[[dataid]]$prob_inc_y01[seq((burn+1),iter,thin),]) 
  pstra_y01 = apply(pstra_y01_all,2,mean)
  
  ### selected each if prob>0.5 to slab, otherwise to spike 
  p_inc11 = ifelse(pstra11>= 0.5,1,0)
  p_inc10 = ifelse(pstra10>= 0.5,1,0)
  p_inc00 = ifelse(pstra00>= 0.5,1,0)
  p_inc_y11 = ifelse(pstra_y11>= 0.5,1,0)
  p_inc_y10 = ifelse(pstra_y10>= 0.5,1,0)
  p_inc_y00 = ifelse(pstra_y00>= 0.5,1,0)
  p_inc_y01 = ifelse(pstra_y01>= 0.5,1,0)

  
  print(dataid)
  
  return(list(p_inc11 = p_inc11, p_inc10 = p_inc10, p_inc00 = p_inc00, p_inc_y11 = p_inc_y11, p_inc_y10 = p_inc_y10, p_inc_y00 = p_inc_y00, p_inc_y01 = p_inc_y01))
}
get_pinclu_12_class_comb_thin10 = mclapply(1,get_pinclu_res,data_class_as_clu_list,iter=iter,burn=burn,mc.cores=20)
save(get_pinclu_12_class_comb_thin10,file="get_pinclu_12_class_comb_thin10.RData")
