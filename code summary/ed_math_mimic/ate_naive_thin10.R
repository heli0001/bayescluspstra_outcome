load("simdata_m32_f.RData")
load("ed_mimic_naivem32f_class_posterior_c1.RData")
load("ed_mimic_naivem32f_class_posterior_c2.RData")
load("ed_mimic_naivem32f_class_posterior_c3.RData")
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)
library(survival)

iter=11000
burn=1000
thin=10
effect_naive_func = function(dataid,data_list,poster_res,poster_res_st2,poster_res_st3,iter=iter,burn=burn){
  set.seed(dataid)
  data_obs = data_list[[dataid]]$data_obs
  N = dim(data_obs)[1]
  ## recall posterior results
  beta00_pos = rbind(poster_res[[dataid]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_res_st2[[dataid]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_res_st3[[dataid]]$beta00_pos[seq((burn+1),iter,thin),]) 
  
  beta10_pos = rbind(poster_res[[dataid]]$beta10_pos[seq((burn+1),iter,thin),],
                 poster_res_st2[[dataid]]$beta10_pos[seq((burn+1),iter,thin),],
                 poster_res_st3[[dataid]]$beta10_pos[seq((burn+1),iter,thin),])
  
  beta11_pos = rbind(poster_res[[dataid]]$beta11_pos[seq((burn+1),iter,thin),],
                 poster_res_st2[[dataid]]$beta11_pos[seq((burn+1),iter,thin),],
                 poster_res_st3[[dataid]]$beta11_pos[seq((burn+1),iter,thin),])                            
  
  alpha00_pos = rbind(poster_res[[dataid]]$alpha00_pos[seq((burn+1),iter,thin),],
                  poster_res_st2[[dataid]]$alpha00_pos[seq((burn+1),iter,thin),],
                  poster_res_st3[[dataid]]$alpha00_pos[seq((burn+1),iter,thin),])
  
  alpha10_pos = rbind(poster_res[[dataid]]$alpha10_pos[seq((burn+1),iter,thin),],
                  poster_res_st2[[dataid]]$alpha10_pos[seq((burn+1),iter,thin),],
                  poster_res_st3[[dataid]]$alpha10_pos[seq((burn+1),iter,thin),]) 
  
  alpha01_pos = rbind(poster_res[[dataid]]$alpha01_pos[seq((burn+1),iter,thin),],
                  poster_res_st2[[dataid]]$alpha01_pos[seq((burn+1),iter,thin),],
                  poster_res_st3[[dataid]]$alpha01_pos[seq((burn+1),iter,thin),])
  
  alpha11_pos = rbind(poster_res[[dataid]]$alpha11_pos[seq((burn+1),iter,thin),],
                  poster_res_st2[[dataid]]$alpha11_pos[seq((burn+1),iter,thin),],
                  poster_res_st3[[dataid]]$alpha11_pos[seq((burn+1),iter,thin),])
  
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
  
  z1_ate = (Ey_1_00*pi_00_est + Ey_1_10*pi_10_est + Ey_1_01*pi_01_est + Ey_1_11*pi_11_est)
  z0_ate = (Ey_0_00*pi_00_est + Ey_0_10*pi_10_est + Ey_0_01*pi_01_est + Ey_0_11*pi_11_est)
  
  z1_1100_ate = (Ey_1_00*pi_00_est + Ey_1_11*pi_11_est)/(pi_00_est + pi_11_est)
  z0_1100_ate = (Ey_0_00*pi_00_est + Ey_0_11*pi_11_est)/(pi_00_est + pi_11_est)
  #ave_ate = exp(mean(z1_ate))-exp(mean(z0_ate))
  #ave_ate = exp(mean(z1_ate) - mean(z0_ate))
  ave_ate = mean(z1_ate) - mean(z0_ate)
  ave_1100_ate = mean(z1_1100_ate) - mean(z0_1100_ate)
  
  ## estimated 
  ate_est = c()
  ate_1100_est = c()
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
    
    z1_ate = (Ey_1_00*pi_00_est + Ey_1_10*pi_10_est + Ey_1_01*pi_01_est + Ey_1_11*pi_11_est)
    z0_ate = (Ey_0_00*pi_00_est + Ey_0_10*pi_10_est + Ey_0_01*pi_01_est + Ey_0_11*pi_11_est)
    
    z1_1100_ate = (Ey_1_00*pi_00_est + Ey_1_11*pi_11_est)/(pi_00_est + pi_11_est)
    z0_1100_ate = (Ey_0_00*pi_00_est + Ey_0_11*pi_11_est)/(pi_00_est + pi_11_est)
    
    ate_1100_est_j = mean(z1_1100_ate) - mean(z0_1100_ate)
    ate_1100_est = c(ate_1100_est, ate_1100_est_j)
    
    #ate_est[j] = exp(mean(z1_ate))-exp(mean(z0_ate))
    ate_est_j = (mean(z1_ate) - mean(z0_ate))
    ate_est = c(ate_est,ate_est_j)
  }
  library(coda)
  ate_mcmc = as.mcmc(ate_est)
  hpd_ate_mcmc = HPDinterval(ate_mcmc, prob = 0.95)
  res_ate = as.data.frame(cbind(mean(ate_est),sd(ate_est),hpd_ate_mcmc[1],hpd_ate_mcmc[2]))
  colnames(res_ate) = c("Mean","SD","2.5%","97.5%")
  
  ate_1100_mcmc = as.mcmc(ate_1100_est)
  hpd_ate_1100_mcmc = HPDinterval(ate_1100_mcmc, prob = 0.95)
  res_ate_1100 = as.data.frame(cbind(mean(ate_1100_est),sd(ate_1100_est),hpd_ate_1100_mcmc[1],hpd_ate_1100_mcmc[2]))
  colnames(res_ate_1100) = c("Mean","SD","2.5%","97.5%")
  
  print(dataid)
  return (list(ate_est=ate_est,ate_1100_est=ate_1100_est,ate_mcmc=ate_mcmc,hpd_ate_mcmc=hpd_ate_mcmc,res_ate=res_ate,ave_ate=ave_ate,ave_1100_ate=ave_1100_ate))
}
effect_naivem32f_class_logscale_thin10_comb = mclapply(1:100,effect_naive_func,simdata_m32_f,ed_mimic_naivem32f_class_posterior_c1,ed_mimic_naivem32f_class_posterior_c2,ed_mimic_naivem32f_class_posterior_c3,iter=iter,burn=burn,mc.cores=20)
save(effect_naivem32f_class_logscale_thin10_comb,file="effect_naivem32f_class_logscale_thin10_comb.RData")

