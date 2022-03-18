load("simdata_m10.RData")
load("poster_lassom10.RData")
load("poster_lassom10_st2.RData")
load("poster_lassom10_st3.RData")
library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)
library(stats)
library(rmutil) # generate laplace
library(statmod) # generate inverse gaussian
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

effect_func = function(dataid,simdata){
  set.seed(dataid)
  data_obs = simdata[[dataid]]$data_obs
  N = dim(data_obs)[1]
  effect_true = simdata[[dataid]]$effect_simu
  ## recall posterior results
  beta00_pos = rbind(poster_lassom10[[dataid]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_lassom10_st2[[dataid]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_lassom10_st3[[dataid]]$beta00_pos[seq((burn+1),iter,thin),]) 
  
  beta10_pos = rbind(poster_lassom10[[dataid]]$beta10_pos[seq((burn+1),iter,thin),],
                 poster_lassom10_st2[[dataid]]$beta10_pos[seq((burn+1),iter,thin),],
                 poster_lassom10_st3[[dataid]]$beta10_pos[seq((burn+1),iter,thin),])
  
  beta11_pos = rbind(poster_lassom10[[dataid]]$beta11_pos[seq((burn+1),iter,thin),],
                 poster_lassom10_st2[[dataid]]$beta11_pos[seq((burn+1),iter,thin),],
                 poster_lassom10_st3[[dataid]]$beta11_pos[seq((burn+1),iter,thin),])                            
  
  alpha00_pos = rbind(poster_lassom10[[dataid]]$alpha00_pos[seq((burn+1),iter,thin),],
                  poster_lassom10_st2[[dataid]]$alpha00_pos[seq((burn+1),iter,thin),],
                  poster_lassom10_st3[[dataid]]$alpha00_pos[seq((burn+1),iter,thin),])
  
  alpha10_pos = rbind(poster_lassom10[[dataid]]$alpha10_pos[seq((burn+1),iter,thin),],
                  poster_lassom10_st2[[dataid]]$alpha10_pos[seq((burn+1),iter,thin),],
                  poster_lassom10_st3[[dataid]]$alpha10_pos[seq((burn+1),iter,thin),]) 
  
  alpha01_pos = rbind(poster_lassom10[[dataid]]$alpha01_pos[seq((burn+1),iter,thin),],
                  poster_lassom10_st2[[dataid]]$alpha01_pos[seq((burn+1),iter,thin),],
                  poster_lassom10_st3[[dataid]]$alpha01_pos[seq((burn+1),iter,thin),])
  
  alpha11_pos = rbind(poster_lassom10[[dataid]]$alpha11_pos[seq((burn+1),iter,thin),],
                  poster_lassom10_st2[[dataid]]$alpha11_pos[seq((burn+1),iter,thin),],
                  poster_lassom10_st3[[dataid]]$alpha11_pos[seq((burn+1),iter,thin),])                                       
  
  zeta_00 = rbind(poster_lassom10[[dataid]]$zeta_00[seq((burn+1),iter,thin),],
              poster_lassom10_st2[[dataid]]$zeta_00[seq((burn+1),iter,thin),],
              poster_lassom10_st3[[dataid]]$zeta_00[seq((burn+1),iter,thin),])
  
  zeta_10 = rbind(poster_lassom10[[dataid]]$zeta_10[seq((burn+1),iter,thin),],
              poster_lassom10_st2[[dataid]]$zeta_10[seq((burn+1),iter,thin),],
              poster_lassom10_st3[[dataid]]$zeta_10[seq((burn+1),iter,thin),]) 
  
  zeta_11 = rbind(poster_lassom10[[dataid]]$zeta_11[seq((burn+1),iter,thin),],
              poster_lassom10_st2[[dataid]]$zeta_11[seq((burn+1),iter,thin),],
              poster_lassom10_st3[[dataid]]$zeta_11[seq((burn+1),iter,thin),]) 
                                         
  
  zeta_y_00 = rbind(poster_lassom10[[dataid]]$zeta_y_00[seq((burn+1),iter,thin),],
                poster_lassom10_st2[[dataid]]$zeta_y_00[seq((burn+1),iter,thin),],
                poster_lassom10_st3[[dataid]]$zeta_y_00[seq((burn+1),iter,thin),]) 
  
  zeta_y_10 = rbind(poster_lassom10[[dataid]]$zeta_y_10[seq((burn+1),iter,thin),],
                poster_lassom10_st2[[dataid]]$zeta_y_10[seq((burn+1),iter,thin),],
                poster_lassom10_st3[[dataid]]$zeta_y_10[seq((burn+1),iter,thin),]) 
  
  zeta_y_01 = rbind(poster_lassom10[[dataid]]$zeta_y_01[seq((burn+1),iter,thin),],
                poster_lassom10_st2[[dataid]]$zeta_y_01[seq((burn+1),iter,thin),],
                poster_lassom10_st3[[dataid]]$zeta_y_01[seq((burn+1),iter,thin),]) 
  
  zeta_y_11 = rbind(poster_lassom10[[dataid]]$zeta_y_11[seq((burn+1),iter,thin),],
                poster_lassom10_st2[[dataid]]$zeta_y_11[seq((burn+1),iter,thin),],
                poster_lassom10_st3[[dataid]]$zeta_y_11[seq((burn+1),iter,thin),])  
  
  
  ## dummy variable/indicator matrix for cluster
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  beta00_ave = apply(beta00_pos,2,mean)
  beta10_ave = apply(beta10_pos,2,mean)
  beta11_ave = apply(beta11_pos,2,mean)
  alpha11_ave = apply(alpha11_pos,2,mean)
  alpha10_ave = apply(alpha10_pos,2,mean)
  alpha01_ave = apply(alpha01_pos,2,mean)
  alpha00_ave = apply(alpha00_pos,2,mean)
  zeta_00_ave =  apply(zeta_00,2,mean)
  zeta_10_ave =  apply(zeta_10,2,mean)
  zeta_11_ave =  apply(zeta_11,2,mean)
  zeta_y_00_ave = apply(zeta_y_00,2,mean)
  zeta_y_10_ave = apply(zeta_y_10,2,mean)
  zeta_y_01_ave = apply(zeta_y_01,2,mean)
  zeta_y_11_ave = apply(zeta_y_11,2,mean)
  
  data_G = cbind(data_obs$x1,data_obs$x2,Istar%*%zeta_11_ave,Istar%*%zeta_10_ave,Istar%*%zeta_00_ave)
  mean_11 = data_G[,1:2]%*%beta11_ave+ data_G[,3]
  mean_10 = data_G[,1:2]%*%beta10_ave+ data_G[,4]
  mean_00 = data_G[,1:2]%*%beta00_ave+ data_G[,5]
  p_11 = apply(mean_11,1,pnorm,mean=0,sd=1)
  p_10 = apply(mean_10,1,pnorm,mean=0,sd=1)
  p_00 = apply(mean_00,1,pnorm,mean=0,sd=1)
  
  pi_11_est = 1-p_11
  pi_10_est = p_11*(1-p_10)
  pi_00_est = p_11*p_10*(1-p_00)
  pi_01_est = p_11*p_10*p_00
  
  data_0_y00 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_00_ave))
  data_1_y00 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_00_ave))
  Ey_0_00 = data_0_y00[,1:3]%*%alpha00_ave + data_0_y00[,4]
  Ey_1_00 = data_1_y00[,1:3]%*%alpha00_ave + data_1_y00[,4]
  ace00 = (Ey_1_00-Ey_0_00)*pi_00_est
  
  data_0_y10 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_10_ave))
  data_1_y10 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_10_ave))
  Ey_0_10 = data_0_y10[,1:3]%*%alpha10_ave + data_0_y10[,4]
  Ey_1_10 = data_1_y10[,1:3]%*%alpha10_ave + data_1_y10[,4]
  ace10 = (Ey_1_10-Ey_0_10)*pi_10_est
  
  data_0_y01 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_01_ave))
  data_1_y01 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_01_ave))
  Ey_0_01 = data_0_y01[,1:3]%*%alpha01_ave + data_0_y01[,4]
  Ey_1_01 = data_1_y01[,1:3]%*%alpha01_ave + data_1_y01[,4]
  ace01 = (Ey_1_01-Ey_0_01)*pi_01_est
  
  data_0_y11 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_11_ave))
  data_1_y11 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_11_ave))
  Ey_0_11 = data_0_y11[,1:3]%*%alpha11_ave + data_0_y11[,4]
  Ey_1_11 = data_1_y11[,1:3]%*%alpha11_ave + data_1_y11[,4]
  ace11 = (Ey_1_11-Ey_0_11)*pi_11_est
  
  ave_ate = mean(ace00+ace01+ace10+ace11)
  
  ## estimated 
  ate_est = c()
  for (i in 1:((iter - burn)/thin*3)){
    
    data_G = cbind(data_obs$x1,data_obs$x2,Istar%*%zeta_11[i,],Istar%*%zeta_10[i,],Istar%*%zeta_00[i,])
    mean_11 = data_G[,1:2]%*%beta11_pos[i,]+ data_G[,3]
    mean_10 = data_G[,1:2]%*%beta10_pos[i,]+ data_G[,4]
    mean_00 = data_G[,1:2]%*%beta00_pos[i,]+ data_G[,5]
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
    
    data_0_y00 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_00[i,]))
    data_1_y00 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_00[i,]))
    Ey_0_00 = data_0_y00[,1:3]%*%alpha00_est + data_0_y00[,4]
    Ey_1_00 = data_1_y00[,1:3]%*%alpha00_est + data_1_y00[,4]
    ace00 = (Ey_1_00-Ey_0_00)*pi_00_est
    
    data_0_y10 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_10[i,]))
    data_1_y10 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_10[i,]))
    Ey_0_10 = data_0_y10[,1:3]%*%alpha10_est + data_0_y10[,4]
    Ey_1_10 = data_1_y10[,1:3]%*%alpha10_est + data_1_y10[,4]
    ace10 = (Ey_1_10-Ey_0_10)*pi_10_est
    
    data_0_y01 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_01[i,]))
    data_1_y01 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_01[i,]))
    Ey_0_01 = data_0_y01[,1:3]%*%alpha01_est + data_0_y01[,4]
    Ey_1_01 = data_1_y01[,1:3]%*%alpha01_est + data_1_y01[,4]
    ace01 = (Ey_1_01-Ey_0_01)*pi_01_est
    
    data_0_y11 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(0,N),Istar%*%zeta_y_11[i,]))
    data_1_y11 = as.matrix(cbind(data_obs$x1,data_obs$x2,rep(1,N),Istar%*%zeta_y_11[i,]))
    Ey_0_11 = data_0_y11[,1:3]%*%alpha11_est + data_0_y11[,4]
    Ey_1_11 = data_1_y11[,1:3]%*%alpha11_est + data_1_y11[,4]
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
effect_lasso_resm10_thin10_comb = mclapply(1:100,effect_func,simdata_m10,mc.cores=20)
save(effect_lasso_resm10_thin10_comb,file="effect_lasso_resm10_thin10_comb.RData")
