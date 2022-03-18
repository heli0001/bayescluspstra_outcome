load("simdata_m10.RData")
load("poster_naive_resm10.RData")
load("poster_naive_resm10_st2.RData")
load("poster_naive_resm10_st3.RData")
library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)
true_para = as.numeric(as.matrix(read.csv("true parameters v3 prins_out.csv",header=T,sep=",",fill=T))[,2])

## for principal stratification
beta00_t = true_para[1:3]
beta10_t = true_para[4:6]
beta11_t = true_para[7:9]

## for outcomes 
alpha00_t = true_para[10:13]
alpha10_t = true_para[14:17]
alpha01_t = true_para[18:21]
alpha11_t = true_para[22:25]
sigsq00_t = true_para[26]
sigsq01_t = true_para[27]
sigsq10_t = true_para[28]
sigsq11_t = true_para[29]


iter=11000
burn=1000
thin = 10
###########################################################################
######################### effect estimand #################################
###########################################################################

effect_naive_func = function(dataid,simdata){
  set.seed(dataid)
  data_obs = simdata[[dataid]]$data_obs
  N = dim(data_obs)[1]
  effect_true = simdata[[dataid]]$effect_simu
  ## recall posterior results
  beta00_pos = rbind(poster_naive_resm10[[dataid]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_naive_resm10_st2[[dataid]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_naive_resm10_st3[[dataid]]$beta00_pos[seq((burn+1),iter,thin),]) 
  beta10_pos = rbind(poster_naive_resm10[[dataid]]$beta10_pos[seq((burn+1),iter,thin),],
               poster_naive_resm10_st2[[dataid]]$beta10_pos[seq((burn+1),iter,thin),],
               poster_naive_resm10_st3[[dataid]]$beta10_pos[seq((burn+1),iter,thin),])
  beta11_pos = rbind(poster_naive_resm10[[dataid]]$beta11_pos[seq((burn+1),iter,thin),],
               poster_naive_resm10_st2[[dataid]]$beta11_pos[seq((burn+1),iter,thin),],
               poster_naive_resm10_st3[[dataid]]$beta11_pos[seq((burn+1),iter,thin),])                            
  
  alpha00_pos = rbind(poster_naive_resm10[[dataid]]$alpha00_pos[seq((burn+1),iter,thin),],
               poster_naive_resm10_st2[[dataid]]$alpha00_pos[seq((burn+1),iter,thin),],
               poster_naive_resm10_st3[[dataid]]$alpha00_pos[seq((burn+1),iter,thin),])
  alpha10_pos = rbind(poster_naive_resm10[[dataid]]$alpha10_pos[seq((burn+1),iter,thin),],
               poster_naive_resm10_st2[[dataid]]$alpha10_pos[seq((burn+1),iter,thin),],
               poster_naive_resm10_st3[[dataid]]$alpha10_pos[seq((burn+1),iter,thin),]) 
  alpha01_pos = rbind(poster_naive_resm10[[dataid]]$alpha01_pos[seq((burn+1),iter,thin),],
               poster_naive_resm10_st2[[dataid]]$alpha01_pos[seq((burn+1),iter,thin),],
               poster_naive_resm10_st3[[dataid]]$alpha01_pos[seq((burn+1),iter,thin),])
  alpha11_pos = rbind(poster_naive_resm10[[dataid]]$alpha11_pos[seq((burn+1),iter,thin),],
               poster_naive_resm10_st2[[dataid]]$alpha11_pos[seq((burn+1),iter,thin),],
               poster_naive_resm10_st3[[dataid]]$alpha11_pos[seq((burn+1),iter,thin),])
  
  beta00_ave = apply(beta00_pos,2,mean)
  beta10_ave = apply(beta10_pos,2,mean)
  beta11_ave = apply(beta11_pos,2,mean)
  alpha11_ave = apply(alpha11_pos,2,mean)
  alpha10_ave = apply(alpha10_pos,2,mean)
  alpha01_ave = apply(alpha01_pos,2,mean)
  alpha00_ave = apply(alpha00_pos,2,mean)
  
  data_G = cbind(data_obs$x1,data_obs$x2)
  mean_11 = data_G%*%beta11_ave
  mean_10 = data_G%*%beta10_ave
  mean_00 = data_G%*%beta00_ave
  p_11 = apply(mean_11,1,pnorm,mean=0,sd=1)
  p_10 = apply(mean_10,1,pnorm,mean=0,sd=1)
  p_00 = apply(mean_00,1,pnorm,mean=0,sd=1)
  
  pi_11_est = 1-p_11
  pi_10_est = p_11*(1-p_10)
  pi_00_est = p_11*p_10*(1-p_00)
  pi_01_est = p_11*p_10*p_00
  
  
  data_0_y = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N)))
  data_1_y = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N)))
  Ey_0_00 = data_0_y%*%alpha00_ave
  Ey_1_00 = data_1_y%*%alpha00_ave
  ace00 = (Ey_1_00-Ey_0_00)*pi_00_est
  
  Ey_0_10 = data_0_y%*%alpha10_ave
  Ey_1_10 = data_1_y%*%alpha10_ave
  ace10 = (Ey_1_10-Ey_0_10)*pi_10_est
  
  
  Ey_0_01 = data_0_y%*%alpha01_ave 
  Ey_1_01 = data_1_y%*%alpha01_ave 
  ace01 = (Ey_1_01-Ey_0_01)*pi_01_est
  
  
  Ey_0_11 = data_0_y%*%alpha11_ave 
  Ey_1_11 = data_1_y%*%alpha11_ave 
  ace11 = (Ey_1_11-Ey_0_11)*pi_11_est
  
  ave_ate = mean(ace00+ace01+ace10+ace11)
  
  
  ## estimated 
  ate_est = c()
 for (i in 1:((iter - burn)/thin*3)){
    
    data_G = cbind(data_obs$x1,data_obs$x2)
    mean_11 = data_G%*%beta11_pos[i,]
    mean_10 = data_G%*%beta10_pos[i,]
    mean_00 = data_G%*%beta00_pos[i,]
    p_11 = apply(mean_11,1,pnorm,mean=0,sd=1)
    p_10 = apply(mean_10,1,pnorm,mean=0,sd=1)
    p_00 = apply(mean_00,1,pnorm,mean=0,sd=1)
    
    pi_11_est = 1-p_11
    pi_10_est = p_11*(1-p_10)
    pi_00_est = p_11*p_10*(1-p_00)
    pi_01_est = p_11*p_10*p_00
    
    alpha00_est = alpha00_pos[i,]
    alpha01_est = alpha01_pos[i,]
    alpha10_est = alpha10_pos[i,]
    alpha11_est = alpha11_pos[i,]
    
    data_0_y = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N)))
    data_1_y = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N)))
    Ey_0_00 = data_0_y%*%alpha00_est
    Ey_1_00 = data_1_y%*%alpha00_est
    ace00 = (Ey_1_00-Ey_0_00)*pi_00_est
    
    Ey_0_10 = data_0_y%*%alpha10_est
    Ey_1_10 = data_1_y%*%alpha10_est
    ace10 = (Ey_1_10-Ey_0_10)*pi_10_est
    
    
    Ey_0_01 = data_0_y%*%alpha01_est 
    Ey_1_01 = data_1_y%*%alpha01_est 
    ace01 = (Ey_1_01-Ey_0_01)*pi_01_est
    
    
    Ey_0_11 = data_0_y%*%alpha11_est 
    Ey_1_11 = data_1_y%*%alpha11_est 
    ace11 = (Ey_1_11-Ey_0_11)*pi_11_est
    
    ate_est_j = mean(ace00+ace01+ace10+ace11)
    ate_est = c(ate_est,ate_est_j)
  }
  library(coda)
  ate_mcmc = as.mcmc(ate_est)
  hpd_ate_mcmc = HPDinterval(ate_mcmc, prob = 0.95)
  res_ate = as.data.frame(cbind(effect_true,mean(ate_est),sd(ate_est),hpd_ate_mcmc[1],hpd_ate_mcmc[2]))
  colnames(res_ate) = c("True","Mean","SD","2.5%","97.5%")
  print(dataid)
  return (list(ate_est=ate_est,ate_mcmc=ate_mcmc,hpd_ate_mcmc=hpd_ate_mcmc,res_ate=res_ate,ave_ate=ave_ate,effect_true=effect_true))
}
effect_naive_resm10_thin10_comb = mclapply(1:100,effect_naive_func,simdata_m10,mc.cores=20)
save(effect_naive_resm10_thin10_comb,file="effect_naive_resm10_thin10_comb.RData")