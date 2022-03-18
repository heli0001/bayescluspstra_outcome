load("bwmimic_poster_nmig_res_mcmc2_single.RData")
#sub_comp = read.csv("more_cluster_induction_with_gesage.csv",header=T,fill=T,sep=",")
#n = matrix(table(sub_comp$cluster))
#sub_comp$id = unlist(apply(n,1,seq,from=1,by=1))
#sub_comp$s = ifelse(sub_comp$ld_indl=="N",0,1)
#sub_comp$log_y =log(sub_comp$dbwt)

## finalize observed dataset
#data_obs = sub_comp[,c(18,19,1,12,13,8,16,20,21)]
#colnames(data_obs) = c("clu","id","x1","x2","x3","x4","z","s_obs","y_obs")


load("simdata_bw_m99.RData")

############## MCMC #################
library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

iter=11000
burn=1000

effect_nmigfunc = function(dataid,poster_res,iter=iter,burn=burn){
  set.seed(dataid)
  data_obs = simdata_bw_m99[[dataid]]$data_obs
  N = dim(data_obs)[1]
  ## recall posterior results
  beta00_pos = poster_res[[dataid]]$beta00_pos
  beta10_pos = poster_res[[dataid]]$beta10_pos
  beta11_pos = poster_res[[dataid]]$beta11_pos
  alpha00_pos = poster_res[[dataid]]$alpha00_pos
  alpha10_pos = poster_res[[dataid]]$alpha10_pos
  alpha01_pos = poster_res[[dataid]]$alpha01_pos
  alpha11_pos = poster_res[[dataid]]$alpha11_pos
  zeta_00 = poster_res[[dataid]]$zeta_00
  zeta_10 = poster_res[[dataid]]$zeta_10
  zeta_11 = poster_res[[dataid]]$zeta_11
  zeta_y_00 = poster_res[[dataid]]$zeta_y_00
  zeta_y_10 = poster_res[[dataid]]$zeta_y_10
  zeta_y_01 = poster_res[[dataid]]$zeta_y_01
  zeta_y_11 = poster_res[[dataid]]$zeta_y_11
  
  ## dummy variable/indicator matrix for cluster
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  beta00_ave = apply(beta00_pos[(burn+1):iter,],2,mean)
  beta10_ave = apply(beta10_pos[(burn+1):iter,],2,mean)
  beta11_ave = apply(beta11_pos[(burn+1):iter,],2,mean)
  alpha11_ave = apply(alpha11_pos[(burn+1):iter,],2,mean)
  alpha10_ave = apply(alpha10_pos[(burn+1):iter,],2,mean)
  alpha01_ave = apply(alpha01_pos[(burn+1):iter,],2,mean)
  alpha00_ave = apply(alpha00_pos[(burn+1):iter,],2,mean)
  zeta_00_ave =  apply(zeta_00[(burn+1):iter,],2,mean)
  zeta_10_ave =  apply(zeta_10[(burn+1):iter,],2,mean)
  zeta_11_ave =  apply(zeta_11[(burn+1):iter,],2,mean)
  zeta_y_00_ave = apply(zeta_y_00[(burn+1):iter,],2,mean)
  zeta_y_10_ave = apply(zeta_y_10[(burn+1):iter,],2,mean)
  zeta_y_01_ave = apply(zeta_y_01[(burn+1):iter,],2,mean)
  zeta_y_11_ave = apply(zeta_y_11[(burn+1):iter,],2,mean)
  
  data_G = cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,Istar%*%zeta_11_ave,Istar%*%zeta_10_ave,Istar%*%zeta_00_ave)
  mean_11 = data_G[,1:4]%*%beta11_ave+ data_G[,5]
  mean_10 = data_G[,1:4]%*%beta10_ave+ data_G[,6]
  mean_00 = data_G[,1:4]%*%beta00_ave+ data_G[,7]
  p_11 = apply(mean_11,1,pnorm,mean=0,sd=1)
  p_10 = apply(mean_10,1,pnorm,mean=0,sd=1)
  p_00 = apply(mean_00,1,pnorm,mean=0,sd=1)
  
  pi_11_est = 1-p_11
  pi_10_est = p_11*(1-p_10)
  pi_00_est = p_11*p_10*(1-p_00)
  pi_01_est = p_11*p_10*p_00
  
  pi_11_est[which(pi_11_est==0)] = 1e-308
  pi_00_est[which(pi_00_est==0)] = 1e-308
  
  data_0_y00 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(0,N),Istar%*%zeta_y_00_ave))
  data_1_y00 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(1,N),Istar%*%zeta_y_00_ave))
  Ey_0_00 = (data_0_y00[,1:5]%*%alpha00_ave + data_0_y00[,6])
  Ey_1_00 = (data_1_y00[,1:5]%*%alpha00_ave + data_1_y00[,6])
  
  data_0_y10 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(0,N),Istar%*%zeta_y_10_ave))
  data_1_y10 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(1,N),Istar%*%zeta_y_10_ave))
  Ey_0_10 = (data_0_y10[,1:5]%*%alpha10_ave + data_0_y10[,6])
  Ey_1_10 = (data_1_y10[,1:5]%*%alpha10_ave + data_1_y10[,6])
  
  data_0_y01 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(0,N),Istar%*%zeta_y_01_ave))
  data_1_y01 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(1,N),Istar%*%zeta_y_01_ave))
  Ey_0_01 = (data_0_y01[,1:5]%*%alpha01_ave + data_0_y01[,6])
  Ey_1_01 = (data_1_y01[,1:5]%*%alpha01_ave + data_1_y01[,6])
  
  data_0_y11 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(0,N),Istar%*%zeta_y_11_ave))
  data_1_y11 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(1,N),Istar%*%zeta_y_11_ave))
  Ey_0_11 = (data_0_y11[,1:5]%*%alpha11_ave + data_0_y11[,6])
  Ey_1_11 = (data_1_y11[,1:5]%*%alpha11_ave + data_1_y11[,6])
  
  #z1_ate = (Ey_1_00*pi_00_est + Ey_1_10*pi_10_est + Ey_1_01*pi_01_est + Ey_1_11*pi_11_est)
  #z0_ate = (Ey_0_00*pi_00_est + Ey_0_10*pi_10_est + Ey_0_01*pi_01_est + Ey_0_11*pi_11_est)
  z1_ate = (Ey_1_00*pi_00_est + Ey_1_11*pi_11_est)/(pi_00_est + pi_11_est)
  z0_ate = (Ey_0_00*pi_00_est + Ey_0_11*pi_11_est)/(pi_00_est + pi_11_est)
  
  #ave_ate = exp(mean(z1_ate))-exp(mean(z0_ate))
  #ave_ate = exp(mean(z1_ate) - mean(z0_ate))
  ave_ate = mean(z1_ate) - mean(z0_ate)
  
  ave_ate_00 = mean(Ey_1_00 - Ey_0_00)
  ## estimated 
  ate_est = rep(-99,(iter-burn))
  ate_11_est = rep(-99,(iter-burn))
  ate_00_est = rep(-99,(iter-burn))
  
  for (i in (burn+1):iter){
    j=i-burn
    
    data_G = cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,Istar%*%zeta_11[i,],Istar%*%zeta_10[i,],Istar%*%zeta_00[i,])
    mean_11 = data_G[,1:4]%*%beta11_pos[i,]+ data_G[,5]
    mean_10 = data_G[,1:4]%*%beta10_pos[i,]+ data_G[,6]
    mean_00 = data_G[,1:4]%*%beta00_pos[i,]+ data_G[,7]
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
    
    data_0_y00 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(0,N),Istar%*%zeta_y_00[i,]))
    data_1_y00 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(1,N),Istar%*%zeta_y_00[i,]))
    Ey_0_00 = (data_0_y00[,1:5]%*%alpha00_est + data_0_y00[,6])
    Ey_1_00 = (data_1_y00[,1:5]%*%alpha00_est + data_1_y00[,6])
    
    
    data_0_y10 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(0,N),Istar%*%zeta_y_10[i,]))
    data_1_y10 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(1,N),Istar%*%zeta_y_10[i,]))
    Ey_0_10 = (data_0_y10[,1:5]%*%alpha10_est + data_0_y10[,6])
    Ey_1_10 = (data_1_y10[,1:5]%*%alpha10_est + data_1_y10[,6])
    
    data_0_y01 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(0,N),Istar%*%zeta_y_01[i,]))
    data_1_y01 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(1,N),Istar%*%zeta_y_01[i,]))
    Ey_0_01 = (data_0_y01[,1:5]%*%alpha01_est + data_0_y01[,6])
    Ey_1_01 = (data_1_y01[,1:5]%*%alpha01_est + data_1_y01[,6])
    
    data_0_y11 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(0,N),Istar%*%zeta_y_11[i,]))
    data_1_y11 = as.matrix(cbind(data_obs$x1,data_obs$x2,data_obs$x3,data_obs$x4,rep(1,N),Istar%*%zeta_y_11[i,]))
    Ey_0_11 = (data_0_y11[,1:5]%*%alpha11_est + data_0_y11[,6])
    Ey_1_11 = (data_1_y11[,1:5]%*%alpha11_est + data_1_y11[,6])
    
    #z1_ate = (Ey_1_00*pi_00_est + Ey_1_10*pi_10_est + Ey_1_01*pi_01_est + Ey_1_11*pi_11_est)
    #z0_ate = (Ey_0_00*pi_00_est + Ey_0_10*pi_10_est + Ey_0_01*pi_01_est + Ey_0_11*pi_11_est)
    z1_ate = (Ey_1_00*pi_00_est + Ey_1_11*pi_11_est)/(pi_00_est + pi_11_est)
    z0_ate = (Ey_0_00*pi_00_est + Ey_0_11*pi_11_est)/(pi_00_est + pi_11_est)
    
    ate_11_est[j] = mean(Ey_1_11 - Ey_0_11)
    ate_00_est[j] = mean(Ey_1_00 - Ey_0_00)
    
    #ate_est[j] = exp(mean(z1_ate))-exp(mean(z0_ate))
    #ate_est[j] = exp(mean(z1_ate) - mean(z0_ate))
    ate_est[j] = mean(z1_ate) - mean(z0_ate)
  }
  library(coda)
  ate_mcmc = as.mcmc(ate_est)
  hpd_ate_mcmc = HPDinterval(ate_mcmc, prob = 0.95)
  res_ate = as.data.frame(cbind(mean(ate_est),sd(ate_est),hpd_ate_mcmc[1],hpd_ate_mcmc[2]))
  colnames(res_ate) = c("Mean","SD","2.5%","97.5%")
  
  ate00_mcmc = as.mcmc(ate_00_est)
  hpd_ate00_mcmc = HPDinterval(ate00_mcmc, prob = 0.95)
  res_ate00 = as.data.frame(cbind(mean(ate_00_est),sd(ate_00_est),hpd_ate00_mcmc[1],hpd_ate00_mcmc[2]))
  colnames(res_ate00) = c("Mean","SD","2.5%","97.5%")
  
  print(dataid)
  return (list(ate_11_est=ate_11_est,ate_00_est=ate_00_est,ate_est=ate_est,ate_mcmc=ate_mcmc,hpd_ate_mcmc=hpd_ate_mcmc,res_ate=res_ate,ave_ate=ave_ate,ate00_mcmc=ate00_mcmc,hpd_ate00_mcmc=hpd_ate00_mcmc,res_ate00=res_ate00,ave_ate_00=ave_ate_00))
}

effect_nmigres_1100 = mclapply(1:100,effect_nmigfunc,bwmimic_poster_nmig_res_mcmc2_single,iter=iter,burn=burn,mc.cores=20)
save(effect_nmigres_1100,file="effect_nmigres_1100.RData")
