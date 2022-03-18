#########################################################
######### combine the first mcmc 3 chains, thin=10 ######
#########################################################

load("simdata_m10.RData")
#load("poster_nmig_resm10_st3_probinclu.RData") # first MCMC result
load("poster_nmig_resm10.RData")
load("poster_nmig_resm10_st2.RData")
load("poster_nmig_resm10_st3.RData")

library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

iter=11000
burn=1000
thin=10

get_pinclu_res = function(dataid,simdata,iter,burn){
  set.seed(dataid)
  dataG = simdata[[dataid]]$data$G
  data_obs = simdata[[dataid]]$data_obs
  N = dim(data_obs)[1]
  m = simdata[[dataid]]$m
  n = simdata[[dataid]]$n
  ui = simdata[[dataid]]$ui
####### recall result from first MCMC 
  pstra11_all = rbind(poster_nmig_resm10[[dataid]]$prob_inc11[seq((burn+1),iter,thin),],
                  poster_nmig_resm10_st2[[dataid]]$prob_inc11[seq((burn+1),iter,thin),],
                  poster_nmig_resm10_st3[[dataid]]$prob_inc11[seq((burn+1),iter,thin),]) 
  pstra11 = apply(pstra11_all,2,mean)
  
  pstra10_all = rbind(poster_nmig_resm10[[dataid]]$prob_inc10[seq((burn+1),iter,thin),],
                poster_nmig_resm10_st2[[dataid]]$prob_inc10[seq((burn+1),iter,thin),],
                poster_nmig_resm10_st3[[dataid]]$prob_inc10[seq((burn+1),iter,thin),]) 
  pstra10 = apply(pstra10_all,2,mean)

  pstra00_all = rbind(poster_nmig_resm10[[dataid]]$prob_inc00[seq((burn+1),iter,thin),],
                poster_nmig_resm10_st2[[dataid]]$prob_inc00[seq((burn+1),iter,thin),],
                poster_nmig_resm10_st3[[dataid]]$prob_inc00[seq((burn+1),iter,thin),]) 
  pstra00 = apply(pstra00_all,2,mean) 

  pstra_y11_all = rbind(poster_nmig_resm10[[dataid]]$prob_inc_y11[seq((burn+1),iter,thin),],
                poster_nmig_resm10_st2[[dataid]]$prob_inc_y11[seq((burn+1),iter,thin),],
                poster_nmig_resm10_st3[[dataid]]$prob_inc_y11[seq((burn+1),iter,thin),]) 
  pstra_y11 = apply(pstra_y11_all,2,mean)
  
  pstra_y10_all = rbind(poster_nmig_resm10[[dataid]]$prob_inc_y10[seq((burn+1),iter,thin),],
                poster_nmig_resm10_st2[[dataid]]$prob_inc_y10[seq((burn+1),iter,thin),],
                poster_nmig_resm10_st3[[dataid]]$prob_inc_y10[seq((burn+1),iter,thin),]) 
  pstra_y10 = apply(pstra_y10_all,2,mean)
  
  pstra_y00_all = rbind(poster_nmig_resm10[[dataid]]$prob_inc_y00[seq((burn+1),iter,thin),],
                poster_nmig_resm10_st2[[dataid]]$prob_inc_y00[seq((burn+1),iter,thin),],
                poster_nmig_resm10_st3[[dataid]]$prob_inc_y00[seq((burn+1),iter,thin),]) 
  pstra_y00 = apply(pstra_y00_all,2,mean)
  
  pstra_y01_all = rbind(poster_nmig_resm10[[dataid]]$prob_inc_y01[seq((burn+1),iter,thin),],
                poster_nmig_resm10_st2[[dataid]]$prob_inc_y01[seq((burn+1),iter,thin),],
                poster_nmig_resm10_st3[[dataid]]$prob_inc_y01[seq((burn+1),iter,thin),]) 
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
get_pinclu_resm10_thincomb = mclapply(1:100,get_pinclu_res,simdata_m10,iter=iter,burn=burn,mc.cores=20)
save(get_pinclu_resm10_thincomb,file="get_pinclu_resm10_thincomb.RData")
