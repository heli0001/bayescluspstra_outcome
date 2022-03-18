load("poster_nmig_resclass12_mcmc2.RData")
load("poster_nmig_resclass12_mcmc2_st2.RData")
load("poster_nmig_resclass12_mcmc2_st3.RData")
load("poster_nmig_resclass12_mcmc2_st4.RData")
load("poster_nmig_resclass12_mcmc2_st5.RData")
load("poster_nmig_resclass12_mcmc2_st6.RData")
load("poster_nmig_resclass12_mcmc2_st7.RData")
load("poster_nmig_resclass12_mcmc2_st8.RData")
load("poster_nmig_resclass12_mcmc2_st9.RData")
load("poster_nmig_resclass12_mcmc2_st10.RData")

load("poster_nmig_resclass34_mcmc2.RData")
load("poster_nmig_resclass34_mcmc2_st2.RData")
load("poster_nmig_resclass34_mcmc2_st3.RData")
load("poster_nmig_resclass34_mcmc2_st4.RData")
load("poster_nmig_resclass34_mcmc2_st5.RData")
load("poster_nmig_resclass34_mcmc2_st6.RData")
load("poster_nmig_resclass34_mcmc2_st7.RData")
load("poster_nmig_resclass34_mcmc2_st8.RData")
load("poster_nmig_resclass34_mcmc2_st9.RData")
load("poster_nmig_resclass34_mcmc2_st10.RData")


load("nmig12_class_posterior.RData")
load("nmig34_class_posterior.RData")
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

effect_nmigfunc = function(subject,dataobs,poster_res,poster_res_st2,poster_res_st3,poster_res_st4,
poster_res_st5, poster_res_st6, poster_res_st7, poster_res_st8, poster_res_st9, poster_res_st10, iter=iter,burn=burn){
  set.seed(subject)
  data_obs = dataobs[[subject]]$data_obs
  N = dim(data_obs)[1]
  ## recall posterior results
  beta00_pos = rbind(poster_res[[subject]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_res_st2[[subject]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_res_st3[[subject]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_res_st4[[subject]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_res_st5[[subject]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_res_st6[[subject]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_res_st7[[subject]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_res_st8[[subject]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_res_st9[[subject]]$beta00_pos[seq((burn+1),iter,thin),],
                 poster_res_st10[[subject]]$beta00_pos[seq((burn+1),iter,thin),])
                  
  beta10_pos = rbind(poster_res[[subject]]$beta10_pos[seq((burn+1),iter,thin),],
                   poster_res_st2[[subject]]$beta10_pos[seq((burn+1),iter,thin),],
                   poster_res_st3[[subject]]$beta10_pos[seq((burn+1),iter,thin),],
                   poster_res_st4[[subject]]$beta10_pos[seq((burn+1),iter,thin),],
                   poster_res_st5[[subject]]$beta10_pos[seq((burn+1),iter,thin),],
                   poster_res_st6[[subject]]$beta10_pos[seq((burn+1),iter,thin),],
                   poster_res_st7[[subject]]$beta10_pos[seq((burn+1),iter,thin),],
                   poster_res_st8[[subject]]$beta10_pos[seq((burn+1),iter,thin),],
                   poster_res_st9[[subject]]$beta10_pos[seq((burn+1),iter,thin),],
                   poster_res_st10[[subject]]$beta10_pos[seq((burn+1),iter,thin),])
                   
  beta11_pos = rbind(poster_res[[subject]]$beta11_pos[seq((burn+1),iter,thin),],
                   poster_res_st2[[subject]]$beta11_pos[seq((burn+1),iter,thin),],
                   poster_res_st3[[subject]]$beta11_pos[seq((burn+1),iter,thin),],
                   poster_res_st4[[subject]]$beta11_pos[seq((burn+1),iter,thin),],
                   poster_res_st5[[subject]]$beta11_pos[seq((burn+1),iter,thin),],
                   poster_res_st6[[subject]]$beta11_pos[seq((burn+1),iter,thin),],
                   poster_res_st7[[subject]]$beta11_pos[seq((burn+1),iter,thin),],
                   poster_res_st8[[subject]]$beta11_pos[seq((burn+1),iter,thin),],
                   poster_res_st9[[subject]]$beta11_pos[seq((burn+1),iter,thin),],
                   poster_res_st10[[subject]]$beta11_pos[seq((burn+1),iter,thin),])                          
  
  alpha00_pos = rbind(poster_res[[subject]]$alpha00_pos[seq((burn+1),iter,thin),],
                   poster_res_st2[[subject]]$alpha00_pos[seq((burn+1),iter,thin),],
                   poster_res_st3[[subject]]$alpha00_pos[seq((burn+1),iter,thin),],
                   poster_res_st4[[subject]]$alpha00_pos[seq((burn+1),iter,thin),],
                   poster_res_st5[[subject]]$alpha00_pos[seq((burn+1),iter,thin),],
                   poster_res_st6[[subject]]$alpha00_pos[seq((burn+1),iter,thin),],
                   poster_res_st7[[subject]]$alpha00_pos[seq((burn+1),iter,thin),],
                   poster_res_st8[[subject]]$alpha00_pos[seq((burn+1),iter,thin),],
                   poster_res_st9[[subject]]$alpha00_pos[seq((burn+1),iter,thin),],
                   poster_res_st10[[subject]]$alpha00_pos[seq((burn+1),iter,thin),])
                   
  alpha10_pos = rbind(poster_res[[subject]]$alpha10_pos[seq((burn+1),iter,thin),],
                   poster_res_st2[[subject]]$alpha10_pos[seq((burn+1),iter,thin),],
                   poster_res_st3[[subject]]$alpha10_pos[seq((burn+1),iter,thin),],
                   poster_res_st4[[subject]]$alpha10_pos[seq((burn+1),iter,thin),],
                   poster_res_st5[[subject]]$alpha10_pos[seq((burn+1),iter,thin),],
                   poster_res_st6[[subject]]$alpha10_pos[seq((burn+1),iter,thin),],
                   poster_res_st7[[subject]]$alpha10_pos[seq((burn+1),iter,thin),],
                   poster_res_st8[[subject]]$alpha10_pos[seq((burn+1),iter,thin),],
                   poster_res_st9[[subject]]$alpha10_pos[seq((burn+1),iter,thin),],
                   poster_res_st10[[subject]]$alpha10_pos[seq((burn+1),iter,thin),])
                   
  alpha01_pos = rbind(poster_res[[subject]]$alpha01_pos[seq((burn+1),iter,thin),],
                   poster_res_st2[[subject]]$alpha01_pos[seq((burn+1),iter,thin),],
                   poster_res_st3[[subject]]$alpha01_pos[seq((burn+1),iter,thin),],
                   poster_res_st4[[subject]]$alpha01_pos[seq((burn+1),iter,thin),],
                   poster_res_st5[[subject]]$alpha01_pos[seq((burn+1),iter,thin),],
                   poster_res_st6[[subject]]$alpha01_pos[seq((burn+1),iter,thin),],
                   poster_res_st7[[subject]]$alpha01_pos[seq((burn+1),iter,thin),],
                   poster_res_st8[[subject]]$alpha01_pos[seq((burn+1),iter,thin),],
                   poster_res_st9[[subject]]$alpha01_pos[seq((burn+1),iter,thin),],
                   poster_res_st10[[subject]]$alpha01_pos[seq((burn+1),iter,thin),])
                   
  alpha11_pos = rbind(poster_res[[subject]]$alpha11_pos[seq((burn+1),iter,thin),],
                   poster_res_st2[[subject]]$alpha11_pos[seq((burn+1),iter,thin),],
                   poster_res_st3[[subject]]$alpha11_pos[seq((burn+1),iter,thin),],
                   poster_res_st4[[subject]]$alpha11_pos[seq((burn+1),iter,thin),],
                   poster_res_st5[[subject]]$alpha11_pos[seq((burn+1),iter,thin),],
                   poster_res_st6[[subject]]$alpha11_pos[seq((burn+1),iter,thin),],
                   poster_res_st7[[subject]]$alpha11_pos[seq((burn+1),iter,thin),],
                   poster_res_st8[[subject]]$alpha11_pos[seq((burn+1),iter,thin),],
                   poster_res_st9[[subject]]$alpha11_pos[seq((burn+1),iter,thin),],
                   poster_res_st10[[subject]]$alpha11_pos[seq((burn+1),iter,thin),])                                     
  
  zeta_00 = rbind(poster_res[[subject]]$zeta_00[seq((burn+1),iter,thin),],
                   poster_res_st2[[subject]]$zeta_00[seq((burn+1),iter,thin),],
                   poster_res_st3[[subject]]$zeta_00[seq((burn+1),iter,thin),],
                   poster_res_st4[[subject]]$zeta_00[seq((burn+1),iter,thin),],
                   poster_res_st5[[subject]]$zeta_00[seq((burn+1),iter,thin),],
                   poster_res_st6[[subject]]$zeta_00[seq((burn+1),iter,thin),],
                   poster_res_st7[[subject]]$zeta_00[seq((burn+1),iter,thin),],
                   poster_res_st8[[subject]]$zeta_00[seq((burn+1),iter,thin),],
                   poster_res_st9[[subject]]$zeta_00[seq((burn+1),iter,thin),],
                   poster_res_st10[[subject]]$zeta_00[seq((burn+1),iter,thin),])
                   
  zeta_10 = rbind(poster_res[[subject]]$zeta_10[seq((burn+1),iter,thin),],
                   poster_res_st2[[subject]]$zeta_10[seq((burn+1),iter,thin),],
                   poster_res_st3[[subject]]$zeta_10[seq((burn+1),iter,thin),],
                   poster_res_st4[[subject]]$zeta_10[seq((burn+1),iter,thin),],
                   poster_res_st5[[subject]]$zeta_10[seq((burn+1),iter,thin),],
                   poster_res_st6[[subject]]$zeta_10[seq((burn+1),iter,thin),],
                   poster_res_st7[[subject]]$zeta_10[seq((burn+1),iter,thin),],
                   poster_res_st8[[subject]]$zeta_10[seq((burn+1),iter,thin),],
                   poster_res_st9[[subject]]$zeta_10[seq((burn+1),iter,thin),],
                   poster_res_st10[[subject]]$zeta_10[seq((burn+1),iter,thin),])
                   
  zeta_11 = rbind(poster_res[[subject]]$zeta_11[seq((burn+1),iter,thin),],
                   poster_res_st2[[subject]]$zeta_11[seq((burn+1),iter,thin),],
                   poster_res_st3[[subject]]$zeta_11[seq((burn+1),iter,thin),],
                   poster_res_st4[[subject]]$zeta_11[seq((burn+1),iter,thin),],
                   poster_res_st5[[subject]]$zeta_11[seq((burn+1),iter,thin),],
                   poster_res_st6[[subject]]$zeta_11[seq((burn+1),iter,thin),],
                   poster_res_st7[[subject]]$zeta_11[seq((burn+1),iter,thin),],
                   poster_res_st8[[subject]]$zeta_11[seq((burn+1),iter,thin),],
                   poster_res_st9[[subject]]$zeta_11[seq((burn+1),iter,thin),],
                   poster_res_st10[[subject]]$zeta_11[seq((burn+1),iter,thin),])
                                         
  zeta_y_00 = rbind(poster_res[[subject]]$zeta_y_00[seq((burn+1),iter,thin),],
                   poster_res_st2[[subject]]$zeta_y_00[seq((burn+1),iter,thin),],
                   poster_res_st3[[subject]]$zeta_y_00[seq((burn+1),iter,thin),],
                   poster_res_st4[[subject]]$zeta_y_00[seq((burn+1),iter,thin),],
                   poster_res_st5[[subject]]$zeta_y_00[seq((burn+1),iter,thin),],
                   poster_res_st6[[subject]]$zeta_y_00[seq((burn+1),iter,thin),],
                   poster_res_st7[[subject]]$zeta_y_00[seq((burn+1),iter,thin),],
                   poster_res_st8[[subject]]$zeta_y_00[seq((burn+1),iter,thin),],
                   poster_res_st9[[subject]]$zeta_y_00[seq((burn+1),iter,thin),],
                   poster_res_st10[[subject]]$zeta_y_00[seq((burn+1),iter,thin),])
                   
  zeta_y_10 = rbind(poster_res[[subject]]$zeta_y_10[seq((burn+1),iter,thin),],
                   poster_res_st2[[subject]]$zeta_y_10[seq((burn+1),iter,thin),],
                   poster_res_st3[[subject]]$zeta_y_10[seq((burn+1),iter,thin),],
                   poster_res_st4[[subject]]$zeta_y_10[seq((burn+1),iter,thin),],
                   poster_res_st5[[subject]]$zeta_y_10[seq((burn+1),iter,thin),],
                   poster_res_st6[[subject]]$zeta_y_10[seq((burn+1),iter,thin),],
                   poster_res_st7[[subject]]$zeta_y_10[seq((burn+1),iter,thin),],
                   poster_res_st8[[subject]]$zeta_y_10[seq((burn+1),iter,thin),],
                   poster_res_st9[[subject]]$zeta_y_10[seq((burn+1),iter,thin),],
                   poster_res_st10[[subject]]$zeta_y_10[seq((burn+1),iter,thin),])
                   
  zeta_y_01 = rbind(poster_res[[subject]]$zeta_y_01[seq((burn+1),iter,thin),],
                   poster_res_st2[[subject]]$zeta_y_01[seq((burn+1),iter,thin),],
                   poster_res_st3[[subject]]$zeta_y_01[seq((burn+1),iter,thin),],
                   poster_res_st4[[subject]]$zeta_y_01[seq((burn+1),iter,thin),],
                   poster_res_st5[[subject]]$zeta_y_01[seq((burn+1),iter,thin),],
                   poster_res_st6[[subject]]$zeta_y_01[seq((burn+1),iter,thin),],
                   poster_res_st7[[subject]]$zeta_y_01[seq((burn+1),iter,thin),],
                   poster_res_st8[[subject]]$zeta_y_01[seq((burn+1),iter,thin),],
                   poster_res_st9[[subject]]$zeta_y_01[seq((burn+1),iter,thin),],
                   poster_res_st10[[subject]]$zeta_y_01[seq((burn+1),iter,thin),])
                   
  zeta_y_11 = rbind(poster_res[[subject]]$zeta_y_11[seq((burn+1),iter,thin),],
                   poster_res_st2[[subject]]$zeta_y_11[seq((burn+1),iter,thin),],
                   poster_res_st3[[subject]]$zeta_y_11[seq((burn+1),iter,thin),],
                   poster_res_st4[[subject]]$zeta_y_11[seq((burn+1),iter,thin),],
                   poster_res_st5[[subject]]$zeta_y_11[seq((burn+1),iter,thin),],
                   poster_res_st6[[subject]]$zeta_y_11[seq((burn+1),iter,thin),],
                   poster_res_st7[[subject]]$zeta_y_11[seq((burn+1),iter,thin),],
                   poster_res_st8[[subject]]$zeta_y_11[seq((burn+1),iter,thin),],
                   poster_res_st9[[subject]]$zeta_y_11[seq((burn+1),iter,thin),],
                   poster_res_st10[[subject]]$zeta_y_11[seq((burn+1),iter,thin),])
  
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
  
  pi_11_est[which(pi_11_est==0)] = 1e-308
  pi_00_est[which(pi_00_est==0)] = 1e-308
  
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
  
  ## all strata
  #z1_ate = (Ey_1_00*pi_00_est + Ey_1_10*pi_10_est + Ey_1_01*pi_01_est + Ey_1_11*pi_11_est)
  #z0_ate = (Ey_0_00*pi_00_est + Ey_0_10*pi_10_est + Ey_0_01*pi_01_est + Ey_0_11*pi_11_est)
  
  ## 1 and 4 only
  z1_ate = (Ey_1_00*pi_00_est + Ey_1_11*pi_11_est)/(pi_00_est + pi_11_est)
  z0_ate = (Ey_0_00*pi_00_est + Ey_0_11*pi_11_est)/(pi_00_est + pi_11_est)
  
  
  #ave_ate = exp(mean(z1_ate))-exp(mean(z0_ate))
  #ave_ate = exp(mean(z1_ate) - mean(z0_ate))
  ave_ate = mean(z1_ate) - mean(z0_ate)
  
  ## estimated 
  ate_est = c()
  for (i in 1:((iter - burn)/thin*10)){
    
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
    
    pi_11_est[which(pi_11_est==0)] = 1e-308
    pi_00_est[which(pi_00_est==0)] = 1e-308
    
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
    
    ## all strata 
    #z1_ate = (Ey_1_00*pi_00_est + Ey_1_10*pi_10_est + Ey_1_01*pi_01_est + Ey_1_11*pi_11_est)
    #z0_ate = (Ey_0_00*pi_00_est + Ey_0_10*pi_10_est + Ey_0_01*pi_01_est + Ey_0_11*pi_11_est)
    
    ## 1 and 4 only
    z1_ate = (Ey_1_00*pi_00_est + Ey_1_11*pi_11_est)/(pi_00_est + pi_11_est)
    z0_ate = (Ey_0_00*pi_00_est + Ey_0_11*pi_11_est)/(pi_00_est + pi_11_est)
    
    ate_est_j = mean(z1_ate) - mean(z0_ate)
    #ate_est[j] = exp(mean(z1_ate) - mean(z0_ate))
    
    ate_est = c(ate_est,ate_est_j)
  }
  library(coda)
  ate_mcmc = as.mcmc(ate_est)
  hpd_ate_mcmc = HPDinterval(ate_mcmc, prob = 0.95)
  res_ate = as.data.frame(cbind(mean(ate_est),sd(ate_est),hpd_ate_mcmc[1],hpd_ate_mcmc[2]))
  colnames(res_ate) = c("Mean","SD","2.5%","97.5%")
  print(subject)
  return (list(ate_est=ate_est,ate_mcmc=ate_mcmc,hpd_ate_mcmc=hpd_ate_mcmc,res_ate=res_ate,ave_ate=ave_ate))
}
effect_nmig12res_cla_logscale_1100_thin10_mcmc2 = mclapply(1:4,effect_nmigfunc,nmig12_class_posterior,poster_nmig_resclass12_mcmc2,poster_nmig_resclass12_mcmc2_st2,poster_nmig_resclass12_mcmc2_st3,poster_nmig_resclass12_mcmc2_st4,poster_nmig_resclass12_mcmc2_st5,poster_nmig_resclass12_mcmc2_st6,poster_nmig_resclass12_mcmc2_st7,poster_nmig_resclass12_mcmc2_st8,poster_nmig_resclass12_mcmc2_st9,poster_nmig_resclass12_mcmc2_st10,iter=iter,burn=burn,mc.cores=4)
save(effect_nmig12res_cla_logscale_1100_thin10_mcmc2,file="effect_nmig12res_cla_logscale_1100_thin10_mcmc2.RData")

effect_nmig34res_cla_logscale_1100_thin10_mcmc2 = mclapply(1:4,effect_nmigfunc,nmig34_class_posterior,poster_nmig_resclass34_mcmc2,poster_nmig_resclass34_mcmc2_st2,poster_nmig_resclass34_mcmc2_st3,poster_nmig_resclass34_mcmc2_st4,poster_nmig_resclass34_mcmc2_st5,poster_nmig_resclass34_mcmc2_st6,poster_nmig_resclass34_mcmc2_st7,poster_nmig_resclass34_mcmc2_st8,poster_nmig_resclass34_mcmc2_st9,poster_nmig_resclass34_mcmc2_st10,iter=iter,burn=burn,mc.cores=4)
save(effect_nmig34res_cla_logscale_1100_thin10_mcmc2,file="effect_nmig34res_cla_logscale_1100_thin10_mcmc2.RData")
