load("data_class_as_clu_list.RData")
load("nmig12_class_posterior.RData") # class as cluster
poster_nmig_first = nmig12_class_posterior


library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

iter=11000
burn=1000

###### here dataid=subject
get_pinclu_res = function(dataid,data_list,poster_nmig_first,iter,burn){
  set.seed(dataid)
  data_obs = data_list[[dataid]]$data_obs
  N = dim(data_obs)[1]
  m = length(unique(data_obs$clu))
  n = as.vector(table(data_obs$clu))
####### recall result from first MCMC 
  pstra11 = apply(poster_nmig_first[[dataid]]$prob_inc11[(burn+1):iter,],2,mean)
  pstra10 = apply(poster_nmig_first[[dataid]]$prob_inc10[(burn+1):iter,],2,mean)
  pstra00 = apply(poster_nmig_first[[dataid]]$prob_inc00[(burn+1):iter,],2,mean)
  pstra_y11 = apply(poster_nmig_first[[dataid]]$prob_inc_y11[(burn+1):iter,],2,mean)
  pstra_y10 = apply(poster_nmig_first[[dataid]]$prob_inc_y10[(burn+1):iter,],2,mean)
  pstra_y00 = apply(poster_nmig_first[[dataid]]$prob_inc_y00[(burn+1):iter,],2,mean)
  pstra_y01 = apply(poster_nmig_first[[dataid]]$prob_inc_y01[(burn+1):iter,],2,mean)
  
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
get_pinclu_12_class_single = mclapply(1:4,get_pinclu_res,data_class_as_clu_list,poster_nmig_first,iter=iter,burn=burn,mc.cores=20)
save(get_pinclu_12_class_single,file="get_pinclu_12_class_single.RData")
