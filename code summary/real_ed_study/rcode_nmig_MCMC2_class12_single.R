load("nmig12_class_posterior.RData")
load("get_pinclu_12_class_single.RData") # first MCMC result
get_pinclu_first = get_pinclu_12_class_single

library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

iter=11000
burn=1000

###########################################################################
######################### MCMC gibb sample ################################
###########################################################################

nmigfunc = function(dataid,iter=100,burn=0){
  set.seed(dataid)
  data_obs = nmig12_class_posterior[[dataid]]$data_obs
  N = dim(data_obs)[1]
  m = length(unique(data_obs$clu))
  n = as.vector(table(data_obs$clu))
  
  ### selected fixed prob from first MCMC
  p_inc11 = get_pinclu_first[[dataid]]$p_inc11
  p_inc10 = get_pinclu_first[[dataid]]$p_inc10
  p_inc00 = get_pinclu_first[[dataid]]$p_inc00
  p_inc_y11 = get_pinclu_first[[dataid]]$p_inc_y11
  p_inc_y10 = get_pinclu_first[[dataid]]$p_inc_y10
  p_inc_y00 = get_pinclu_first[[dataid]]$p_inc_y00
  p_inc_y01 = get_pinclu_first[[dataid]]$p_inc_y01
  
  
  ########### MCMC posterior samples
  # principal strata
  beta00_pos = matrix(-99,nrow=iter,ncol=2)
  beta10_pos = matrix(-99,nrow=iter,ncol=2)
  beta11_pos = matrix(-99,nrow=iter,ncol=2)
  var_beta00 = rep(-99,iter)
  var_beta10 = rep(-99,iter)
  var_beta11 = rep(-99,iter)
  
  # outcomes
  alpha00_pos = matrix(-99,nrow=iter,ncol=3)
  alpha10_pos = matrix(-99,nrow=iter,ncol=3)
  alpha01_pos = matrix(-99,nrow=iter,ncol=3)
  alpha11_pos = matrix(-99,nrow=iter,ncol=3)
  var_alpha11 = rep(-99,iter)
  var_alpha10 = rep(-99,iter)
  var_alpha00 = rep(-99,iter)
  var_alpha01 = rep(-99,iter)
  
  sigsq00_pos = rep(-99,iter)
  sigsq01_pos = rep(-99,iter)
  sigsq10_pos = rep(-99,iter)
  sigsq11_pos = rep(-99,iter)
  
  ## hyper cluster parameters
  zeta_11 = matrix(-99,nrow=iter,ncol=m)
  zeta_10 = matrix(-99,nrow=iter,ncol=m)
  zeta_00 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_11 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_10 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_01 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_00 = matrix(-99,nrow=iter,ncol=m)
  
  ## cluster-variable indicator parameters
  delta11 = matrix(p_inc11, nrow=iter, ncol=m, byrow=TRUE)
  delta10 = matrix(p_inc10, nrow=iter, ncol=m, byrow=TRUE)
  delta00 = matrix(p_inc00, nrow=iter, ncol=m, byrow=TRUE)
  delta_y11 = matrix(p_inc_y11, nrow=iter, ncol=m, byrow=TRUE)
  delta_y10 = matrix(p_inc_y10, nrow=iter, ncol=m, byrow=TRUE)
  delta_y01 = matrix(p_inc_y01, nrow=iter, ncol=m, byrow=TRUE)
  delta_y00 = matrix(p_inc_y00, nrow=iter, ncol=m, byrow=TRUE)
  
  ## probablity inclustion (use fixed from 1 MCMC)
  prob_inc11 = matrix(p_inc11, nrow=iter, ncol=m, byrow=TRUE)
  prob_inc10 = matrix(p_inc10, nrow=iter, ncol=m, byrow=TRUE)
  prob_inc00 = matrix(p_inc00, nrow=iter, ncol=m, byrow=TRUE)
  prob_inc_y11 = matrix(p_inc_y11, nrow=iter, ncol=m, byrow=TRUE)
  prob_inc_y10 = matrix(p_inc_y10, nrow=iter, ncol=m, byrow=TRUE)
  prob_inc_y00 = matrix(p_inc_y00, nrow=iter, ncol=m, byrow=TRUE)
  prob_inc_y01 = matrix(p_inc_y01, nrow=iter, ncol=m, byrow=TRUE)
  
  ## variance of spike and slab distribution
  phi1 = rep(-99,iter) # for principal strata
  phi2 = rep(-99,iter) # for outcome
  
  ## starting value
  # principal strata
  pstra_m = glm(s_obs~-1+x1+x2+as.factor(clu),family=binomial(link = "probit"),data=data_obs)
  beta00_pos[1,] = coef(pstra_m)[1:2]
  beta10_pos[1,] = coef(pstra_m)[1:2]
  beta11_pos[1,] = coef(pstra_m)[1:2]
  var_beta00[1] = 100
  var_beta10[1] = 100
  var_beta11[1] = 100
  
  # outcomes
  out_m = lm(y_obs~-1+x1+x2+z+as.factor(clu),data=data_obs)
  alpha00_pos[1,] = coef(out_m)[1:3]
  alpha10_pos[1,] = coef(out_m)[1:3]
  alpha01_pos[1,] = coef(out_m)[1:3]
  alpha11_pos[1,] = coef(out_m)[1:3]
  var_alpha11[1] = 100
  var_alpha10[1] = 100
  var_alpha00[1] = 100
  var_alpha01[1] = 100
  
  sigsq00_pos[1] = sum(out_m$residuals^2)/(N-3)
  sigsq01_pos[1] = sum(out_m$residuals^2)/(N-3)
  sigsq10_pos[1] = sum(out_m$residuals^2)/(N-3)
  sigsq11_pos[1] = sum(out_m$residuals^2)/(N-3)
  
  ## hyper lassonaive parameters
  zeta_11[1,] = rep(min(coef(pstra_m)[3:(m+2)],na.rm=T),m)
  zeta_10[1,] = rep(min(coef(pstra_m)[3:(m+2)],na.rm=T),m)
  zeta_00[1,] = rep(min(coef(pstra_m)[3:(m+2)],na.rm=T),m)
  zeta_y_11[1,] = rep(min(coef(out_m)[4:(m+3)],na.rm=T),m)
  zeta_y_10[1,] = rep(min(coef(out_m)[4:(m+3)],na.rm=T),m)
  zeta_y_01[1,] = rep(min(coef(out_m)[4:(m+3)],na.rm=T),m)
  zeta_y_00[1,] = rep(min(coef(out_m)[4:(m+3)],na.rm=T),m)

  
  ## cluster-variable indicator parameters (no need)
  # delta11[1,] = c(rep(0,m0),rep(1,(m-m0)))
  # delta10[1,] = c(rep(0,m0),rep(1,(m-m0)))
  # delta00[1,] = c(rep(0,m0),rep(1,(m-m0)))
  # delta_y11[1,] = c(rep(0,m0),rep(1,(m-m0)))
  # delta_y10[1,] = c(rep(0,m0),rep(1,(m-m0)))
  # delta_y01[1,] = c(rep(0,m0),rep(1,(m-m0)))
  # delta_y00[1,] = c(rep(0,m0),rep(1,(m-m0)))
  
  ## probablity inclustion (no need)
  # prob_inc11[1,] = rep(0.5,m)
  # prob_inc10[1,] = rep(0.5,m)
  # prob_inc00[1,] = rep(0.5,m)
  # prob_inc_y11[1,] = rep(0.5,m)
  # prob_inc_y10[1,] = rep(0.5,m)
  # prob_inc_y01[1,] = rep(0.5,m)
  # prob_inc_y00[1,] = rep(0.5,m)
  
  phi1[1] = 100
  phi2[1] = 100
  
  ## dummy variable/indicator matrix for cluster
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  for (k in 2:iter){
  # Step 1: update principal strata parameters
  # 1-0 : sample G
  data_G = cbind(data_obs$x1,data_obs$x2,Istar%*%zeta_11[k-1,],Istar%*%zeta_10[k-1,],Istar%*%zeta_00[k-1,])
  mean_11 = data_G[,1:2]%*%beta11_pos[k-1,]+ data_G[,3]
  mean_10 = data_G[,1:2]%*%beta10_pos[k-1,]+ data_G[,4]
  mean_00 = data_G[,1:2]%*%beta00_pos[k-1,]+ data_G[,5]
  p_11 = apply(mean_11,1,pnorm,mean=0,sd=1)
  p_10 = apply(mean_10,1,pnorm,mean=0,sd=1)
  p_00 = apply(mean_00,1,pnorm,mean=0,sd=1)
  
  pi_11 = 1-p_11
  pi_10 = p_11*(1-p_10)
  pi_00 = p_11*p_10*(1-p_00)
  pi_01 = p_11*p_10*p_00
  
  fi_func = function(data_yg,alpha_pos,sigsq){
    # data_yg includes y_obs,x1,x2,zeta_yg
    yobs = data_yg[1]
    data_stra = data_yg[2:3]
    fi0 = dnorm(yobs,mean=data_stra%*%alpha_pos[1:2]+data_yg[4],sd=sqrt(sigsq))
    fi1 = dnorm(yobs,mean=data_stra%*%alpha_pos[1:2]+alpha_pos[3]+data_yg[4],sd=sqrt(sigsq))
    return(c(fi0,fi1))
  }
  fi_0_00 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_00[k-1,])),1,fi_func,alpha00_pos[k-1,],sigsq00_pos[k-1])[1,] 
  fi_0_01 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_01[k-1,])),1,fi_func,alpha01_pos[k-1,],sigsq01_pos[k-1])[1,] 
  fi_0_10 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_10[k-1,])),1,fi_func,alpha10_pos[k-1,],sigsq10_pos[k-1])[1,] 
  fi_0_11 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_11[k-1,])),1,fi_func,alpha11_pos[k-1,],sigsq11_pos[k-1])[1,] 
  
  fi_1_00 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_00[k-1,])),1,fi_func,alpha00_pos[k-1,],sigsq00_pos[k-1])[2,] 
  fi_1_01 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_01[k-1,])),1,fi_func,alpha01_pos[k-1,],sigsq01_pos[k-1])[2,] 
  fi_1_10 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_10[k-1,])),1,fi_func,alpha10_pos[k-1,],sigsq10_pos[k-1])[2,] 
  fi_1_11 = apply(as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2,Istar%*%zeta_y_11[k-1,])),1,fi_func,alpha11_pos[k-1,],sigsq11_pos[k-1])[2,] 
  
  p_me = matrix(0,nrow=N,ncol=4) # intermediate calculation
  p_s = matrix(0,nrow=N,ncol=4)      # 4 strata probabilities for each unit
  for (i in 1:N){
    if (data_obs$z[i]==0 &data_obs$s_obs[i]==0) {
      p_me[i,4] = pi_01[i]*fi_0_01[i]
      p_me[i,3] = pi_00[i]*fi_0_00[i]
    }else if (data_obs$z[i]==0 &data_obs$s_obs[i]==1){
      p_me[i,2] = pi_10[i]*fi_0_10[i]
      p_me[i,1] = pi_11[i]*fi_0_11[i]
    }else if (data_obs$z[i]==1 &data_obs$s_obs[i]==0){
      p_me[i,3] = pi_00[i]*fi_1_00[i]
      p_me[i,2] = pi_10[i]*fi_1_10[i]
    }else if (data_obs$z[i]==1 &data_obs$s_obs[i]==1){
      p_me[i,4] = pi_01[i]*fi_1_01[i]
      p_me[i,1] = pi_11[i]*fi_1_11[i]
    }
    
    for (j in 1:4){
      p_s[i,j] = p_me[i,j]/sum(p_me[i,])
    }
  }
  
  strata_num = rep(-99,N)
  for (i in 1:N){
    if (data_obs$z[i]==1&data_obs$s_obs[i]==1){
      strata_num[i] = sample(c(1,4),size=1,replace=TRUE,prob=p_s[i,c(1,4)])
    }else if (data_obs$z[i]==1&data_obs$s_obs[i]==0){
      strata_num[i] = sample(2:3,size=1,replace=TRUE,prob=p_s[i,2:3])
    }else if (data_obs$z[i]==0&data_obs$s_obs[i]==1){
      strata_num[i] = sample(1:2,size=1,replace=TRUE,prob=p_s[i,1:2]) 
    }else if (data_obs$z[i]==0&data_obs$s_obs[i]==0){
      strata_num[i] = sample(3:4,size=1,replace=TRUE,prob=p_s[i,3:4])
    }
  }  
  cl_g_11 = Istar%*%zeta_11[k-1,]
  cl_g_10 = Istar%*%zeta_10[k-1,]
  cl_g_00 = Istar%*%zeta_00[k-1,]
  data_obs2 = as.data.frame(cbind(data_obs,strata_num,cl_g_11,cl_g_10,cl_g_00))
  
  # 1-1: sample latent variable G*
  for (i in 1:N){
    if (data_obs2$strata_num[i]==1){
      data_obs2$gstar_11[i] = rtruncnorm(1,a = -Inf,b=0,mean=mean_11[i],sd=1)
    }else if (data_obs2$strata_num[i]!=1){data_obs2$gstar_11[i] = rtruncnorm(1,a =0,b=Inf,mean=mean_11[i],sd=1)}
    if (data_obs2$strata_num[i]==2){
      data_obs2$gstar_10[i] = rtruncnorm(1,a = -Inf,b=0,mean=mean_10[i],sd=1)
    }else if (data_obs2$strata_num[i]!=2){data_obs2$gstar_10[i] = rtruncnorm(1,a =0,b=Inf,mean=mean_10[i],sd=1)}
    if (data_obs2$strata_num[i]==3){
      data_obs2$gstar_00[i] = rtruncnorm(1,a = -Inf,b=0,mean=mean_00[i],sd=1)
    }else if (data_obs2$strata_num[i]!=3){data_obs2$gstar_00[i] = rtruncnorm(1,a =0,b=Inf,mean=mean_00[i],sd=1)}
  }
  # 1-2: update beta_11,beta_10,beta_00
  #(a) beta_11
  data_pristr = cbind(data_obs2$x1,data_obs2$x2)
  beta11_Sigma = solve(1/var_beta11[k-1]*diag(2)+t(data_pristr)%*%data_pristr)
  beta11_mu = beta11_Sigma%*%(t(data_pristr)%*%(data_obs2$gstar_11-data_obs2$cl_g_11))
  beta11_pos[k,] = mvrnorm(1,beta11_mu,beta11_Sigma)
  
  
  #(b) beta_10
  datano11 = data_obs2[data_obs2$strata_num!=1,]
  data_pristr = cbind(datano11$x1,datano11$x2)
  beta10_Sigma = solve(1/var_beta10[k-1]*diag(2)+t(data_pristr)%*%data_pristr)
  beta10_mu = beta10_Sigma%*%(t(data_pristr)%*%(datano11$gstar_10-datano11$cl_g_10))
  beta10_pos[k,] = mvrnorm(1,beta10_mu,beta10_Sigma)
  
  #(c) beta_00
  datano1110 = datano11[datano11$strata_num!=2,]
  data_pristr = cbind(datano1110$x1,datano1110$x2)
  beta00_Sigma = solve(1/var_beta00[k-1]*diag(2)+t(data_pristr)%*%data_pristr)
  beta00_mu = beta00_Sigma%*%(t(data_pristr)%*%(datano1110$gstar_00-datano1110$cl_g_00))
  beta00_pos[k,] = mvrnorm(1,beta00_mu,beta00_Sigma)
  
  # 1-3 var_beta
  var_beta00[k] = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(beta00_pos[k,])%*%beta00_pos[k,]/2)
  var_beta10[k] = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(beta10_pos[k,])%*%beta10_pos[k,]/2)
  var_beta11[k] = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(beta11_pos[k,])%*%beta11_pos[k,]/2)
  
  # 1-4 zeta_g and delta_g, 1,2,...m 
  #(a-1) zeta_11
  ## check variable inclusion
  ind_0 = which(delta11[k-1,]==0) # zero
  ind_n0 = which(delta11[k-1,]!=0) # nonzero
  rzeta11 = rep(1,m)
  rzeta11[ind_0] = 1/0.00025 ## specify small r
  C = Istar
  sigma_zeta11 = diag(rzeta11*1/phi1[k-1],m,m) # specify phi as variance
  
  data_cl_all = data_obs2
  data_cl = cbind(data_cl_all$x1,data_cl_all$x2) #x1,x2
  cl_Sigma11 = solve(sigma_zeta11+t(C)%*%C)
  cl_mu11 = cl_Sigma11%*%t(C)%*%(data_cl_all$gstar_11-data_cl%*%beta11_pos[k,])
  zeta_11[k,] = mvrnorm(1,cl_mu11,cl_Sigma11)
  
  # #(a-2) delta_11
  # pslab = dnorm(zeta_11[k,],0,sqrt(phi1[k-1]))*prob_inc11[k-1,]
  # pspike = dnorm(zeta_11[k,],0,sqrt(0.00025*phi1[k-1]))*(1-prob_inc11[k-1,])
  # delta11[k,] = rbinom(m,1,pslab/(pspike+pslab))
  
  #(b-1) zeta_10
  data_Istar = as.data.frame(cbind(strata_num,Istar))
  
  data_cl_all = data_obs2[data_obs2$strata_num!=1,]
  Istarno11 = as.matrix(data_Istar[data_Istar$strata_num!=1,-1]) ## only left the 10,00,01 data
  data_cl = cbind(data_cl_all$x1,data_cl_all$x2)  #x1,x2
  
  ind_0 = which(delta10[k-1,]==0) # zero
  ind_n0 = which(delta10[k-1,]!=0) # nonzero
  rzeta10 = rep(1,m)
  rzeta10[ind_0] = 1/0.00025 ## specify small r
  C = Istarno11
  sigma_zeta10 = diag(rzeta10*1/phi1[k-1],m,m)
  
  cl_Sigma10 = solve(sigma_zeta10+t(C)%*%C)
  cl_mu10 = cl_Sigma10%*%t(C)%*%(data_cl_all$gstar_10-data_cl%*%beta10_pos[k,])
  zeta_10[k,] = mvrnorm(1,cl_mu10,cl_Sigma10)
  
  # #(b-1) delta10
  # pslab = dnorm(zeta_10[k,],0,sqrt(phi1[k-1]))*prob_inc10[k-1,]
  # pspike = dnorm(zeta_10[k,],0,sqrt(0.00025*phi1[k-1]))*(1-prob_inc10[k-1,])
  # delta10[k,] = rbinom(m,1,pslab/(pspike+pslab))
  
  #(c-1) zeta_00
  data_cl_all = data_obs2[data_obs2$strata_num==3|data_obs2$strata_num==4,]
  Istarno1110 = as.matrix(data_Istar[data_Istar$strata_num==3|data_Istar$strata_num==4,-1]) ## only left the 00,01 data
  data_cl = cbind(data_cl_all$x1,data_cl_all$x2)  #x1,x2
  
  ind_0 = which(delta00[k-1,]==0) # zero
  ind_n0 = which(delta00[k-1,]!=0) # nonzero
  rzeta00 = rep(1,m)
  rzeta00[ind_0] = 1/0.00025 ## specify small r
  C = Istarno1110
  sigma_zeta00 = diag(rzeta00*1/phi1[k-1],m,m) # specify 10 as variance
  
  
  cl_Sigma00 = solve(sigma_zeta00+t(C)%*%C)
  cl_mu00 = cl_Sigma00%*%t(C)%*%(data_cl_all$gstar_00-data_cl%*%beta00_pos[k,])
  zeta_00[k,] = mvrnorm(1,cl_mu00,cl_Sigma00)
  
  # #(c-1) delta00
  # pslab = dnorm(zeta_00[k,],0,sqrt(phi1[k-1]))*prob_inc00[k-1,]
  # pspike = dnorm(zeta_00[k,],0,sqrt(0.00025*phi1[k-1]))*(1-prob_inc00[k-1,])
  # delta00[k,] = rbinom(m,1,pslab/(pspike+pslab))
  # 
  # # 1-5 update prob_incs 
  # prob_inc11[k,] = rbeta(m,1+delta11[k,],1+(1-delta11[k,]))
  # prob_inc10[k,] = rbeta(m,1+delta10[k,],1+(1-delta10[k,]))
  # prob_inc00[k,] = rbeta(m,1+delta00[k,],1+(1-delta00[k,]))
  
  ## step 2 : update outcome parameters
  # 2-0 : add data information 
  clu_y_11 = Istar%*%zeta_y_11[k-1,]
  clu_y_10 = Istar%*%zeta_y_10[k-1,]
  clu_y_00 = Istar%*%zeta_y_00[k-1,]
  clu_y_01 = Istar%*%zeta_y_01[k-1,]
  data_obs3 =as.data.frame(cbind(data_obs2,clu_y_11,clu_y_10,clu_y_00,clu_y_01))
  
  # 2-1 : alpha_g
  update_alpha_func = function(data,sigsq_alpha,sig){
    if (dim(data)[1]==1){
      var_alpha = solve(1/sigsq_alpha*diag(3)+matrix(data[,1:3])%*%data[,1:3]/sig)
      mu_alpha = var_alpha%*%(matrix(data[,1:3])%*%(data[,5]-data[,4])/sig)
      return (mvrnorm(1,mu=mu_alpha,Sigma=var_alpha))
    } else if (dim(data)[1]>1){
      var_alpha = solve(1/sigsq_alpha*diag(3)+t(data[,1:3])%*%data[,1:3]/sig)
      mu_alpha = var_alpha%*%(t(data[,1:3])%*%(data[,5]-data[,4])/sig)
      return (mvrnorm(1,mu=mu_alpha,Sigma=var_alpha))
    }
  }
  
  # (a) : alpha_11
  data_11 = as.matrix(data_obs3[data_obs3$strata_num==1,c(3,4,5,15,7)]) #x1,x2,z,zeta_y11,y_obs
  if (dim(data_11)[1]==0){
    alpha11_pos[k,] = alpha11_pos[k-1,]
  }else{
    alpha11_pos[k,] = update_alpha_func(data_11,var_alpha11[k-1],sigsq11_pos[k-1])
  }
  
  # (b) : alpha_10
  data_10 = as.matrix(data_obs3[data_obs3$strata_num==2,c(3,4,5,16,7)]) #x1,x2,z,zeta_y10,y_obs
  if (dim(data_10)[1]==0){
    alpha10_pos[k,] = alpha10_pos[k-1,]
  }else{
    alpha10_pos[k,] = update_alpha_func(data_10,var_alpha10[k-1],sigsq10_pos[k-1])
  }
  # (c) : alpha_00
  data_00 = as.matrix(data_obs3[data_obs3$strata_num==3,c(3,4,5,17,7)]) #x1,x2,z,zeta_y00,y_obs
  if (dim(data_00)[1]==0){
    alpha00_pos[k,] = alpha00_pos[k-1,]
  }else{
    alpha00_pos[k,] = update_alpha_func(data_00,var_alpha00[k-1],sigsq00_pos[k-1])
  }
  # (d) : alpha_01
  data_01 = as.matrix(data_obs3[data_obs3$strata_num==4,c(3,4,5,18,7)]) #x1,x2,z,zeta_y01,y_obs
  if (dim(data_01)[1]==0){
    alpha01_pos[k,] = alpha01_pos[k-1,]
  }else{
    alpha01_pos[k,] = update_alpha_func(data_01,var_alpha01[k-1],sigsq01_pos[k-1])
  }
  
  # 2-2: var_alpha_g
  var_alpha11[k] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(alpha11_pos[k,])%*%alpha11_pos[k,]/2)
  var_alpha10[k] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(alpha10_pos[k,])%*%alpha10_pos[k,]/2)
  var_alpha01[k] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(alpha01_pos[k,])%*%alpha01_pos[k,]/2)
  var_alpha00[k] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(alpha00_pos[k,])%*%alpha00_pos[k,]/2)
  
  # 2-3: sigsq_g
  ry11 = data_11[,5]-data_11[,4]-data_11[,1:3]%*%alpha11_pos[k,]
  shape_sigsq11 = (0.002+nrow(ry11))/2
  scale_sigsq11 = (0.002*1 + t(ry11)%*%ry11)/2
  sigsq11_pos[k] = rinvgamma(1,shape=shape_sigsq11,rate=scale_sigsq11)
  
  ry10 = data_10[,5]-data_10[,4]-data_10[,1:3]%*%alpha10_pos[k,]
  shape_sigsq10 = (0.002+nrow(ry10))/2
  scale_sigsq10 = (0.002*1 + t(ry10)%*%ry10)/2
  sigsq10_pos[k] = rinvgamma(1,shape=shape_sigsq10,rate=scale_sigsq10)
  
  ry00 = data_00[,5]-data_00[,4]-data_00[,1:3]%*%alpha00_pos[k,]
  shape_sigsq00 = (0.002+nrow(ry00))/2
  scale_sigsq00 = (0.002*1 + t(ry00)%*%ry00)/2
  sigsq00_pos[k] = rinvgamma(1,shape=shape_sigsq00,rate=scale_sigsq00)
  
  ry01 = data_01[,5]-data_01[,4]-data_01[,1:3]%*%alpha01_pos[k,]
  shape_sigsq01 = (0.002+nrow(ry01))/2
  scale_sigsq01 = (0.002*1 + t(ry01)%*%ry01)/2
  sigsq01_pos[k] = rinvgamma(1,shape=shape_sigsq01,rate=scale_sigsq01)
  
  # 2-4: zeta_yg 1,2,...m
  # (a-1) zeta_y11
  data_cl_11 = as.matrix(data_obs3[data_obs3$strata_num==1,c(3,4,5,7)])
  Istar11 = as.matrix(data_Istar[data_Istar$strata_num==1,-1])
  
  ind_0 = which(delta_y11[k-1,]==0) # zero
  ind_n0 = which(delta_y11[k-1,]!=0) # nonzero
  rzeta_y11 = rep(1,m)
  rzeta_y11[ind_0] = 1/0.00025 ## specify small r
  C = Istar11
  sigma_zetay11 = diag(rzeta_y11*1/phi2[k-1],m,m) 
  
  var_zetay11 = solve(sigma_zetay11+t(C)%*%C/sigsq11_pos[k])
  mu_zetay11 = var_zetay11%*%(t(C)%*%(data_cl_11[,4]-data_cl_11[,1:3]%*%alpha11_pos[k,])/sigsq11_pos[k])
  zeta_y_11[k,] = mvrnorm(1,mu_zetay11,var_zetay11)
  
  # # (a-2) delta_y11
  # pslab = dnorm(zeta_y_11[k,],0,sqrt(phi2[k-1]))*prob_inc_y11[k-1,]
  # pspike = dnorm(zeta_y_11[k,],0,sqrt(0.00025*phi2[k-1]))*(1-prob_inc_y11[k-1,])
  # delta_y11[k,] = rbinom(m,1,pslab/(pspike+pslab))
  
  # (b-1) zeta_y10
  data_cl_10 = as.matrix(data_obs3[data_obs3$strata_num==2,c(3,4,5,7)])
  Istar10 = as.matrix(data_Istar[data_Istar$strata_num==2,-1])
  
  ind_0 = which(delta_y10[k-1,]==0) # zero
  ind_n0 = which(delta_y10[k-1,]!=0) # nonzero
  rzeta_y10 = rep(1,m)
  rzeta_y10[ind_0] = 1/0.00025 ## specify small r
  C = Istar10
  sigma_zetay10 = diag(rzeta_y10*1/phi2[k-1],m,m) # specify phi as variance
  
  var_zetay10 = solve(sigma_zetay10+t(C)%*%C/sigsq10_pos[k])
  mu_zetay10 = var_zetay10%*%(t(C)%*%(data_cl_10[,4]-data_cl_10[,1:3]%*%alpha10_pos[k,])/sigsq10_pos[k])
  zeta_y_10[k,] = mvrnorm(1,mu_zetay10,var_zetay10)
  
  # # (b-2) delta_y10
  # pslab = dnorm(zeta_y_10[k,],0,sqrt(phi2[k-1]))*prob_inc_y10[k-1,]
  # pspike = dnorm(zeta_y_10[k,],0,sqrt(0.00025*phi2[k-1]))*(1-prob_inc_y10[k-1,])
  # delta_y10[k,] = rbinom(m,1,pslab/(pspike+pslab))
  
  
  # (c-1) zeta_y00
  data_cl_00 = as.matrix(data_obs3[data_obs3$strata_num==3,c(3,4,5,7)])
  Istar00 = as.matrix(data_Istar[data_Istar$strata_num==3,-1])
  
  ind_0 = which(delta_y00[k-1,]==0) # zero
  ind_n0 = which(delta_y00[k-1,]!=0) # nonzero
  rzeta_y00 = rep(1,m)
  rzeta_y00[ind_0] = 1/0.00025 ## specify small r
  C = Istar00
  sigma_zetay00 = diag(rzeta_y00*1/phi2[k-1],m,m) # specify phi as variance
  
  var_zetay00 = solve(sigma_zetay00+t(C)%*%C/sigsq00_pos[k])
  mu_zetay00 = var_zetay00%*%(t(C)%*%(data_cl_00[,4]-data_cl_00[,1:3]%*%alpha00_pos[k,])/sigsq00_pos[k])
  zeta_y_00[k,] = mvrnorm(1,mu_zetay00,var_zetay00)
  
  # # (c-2) delta_y00
  # pslab = dnorm(zeta_y_00[k,],0,sqrt(phi2[k-1]))*prob_inc_y00[k-1,]
  # pspike = dnorm(zeta_y_00[k,],0,sqrt(0.00025*phi2[k-1]))*(1-prob_inc_y00[k-1,])
  # delta_y00[k,] = rbinom(m,1,pslab/(pspike+pslab))
  
  # (d-1) zeta_y01
  data_cl_01 = as.matrix(data_obs3[data_obs3$strata_num==4,c(3,4,5,7)])
  Istar01 = as.matrix(data_Istar[data_Istar$strata_num==4,-1])
  
  ind_0 = which(delta_y01[k-1,]==0) # zero
  ind_n0 = which(delta_y01[k-1,]!=0) # nonzero
  rzeta_y01 = rep(1,m)
  rzeta_y01[ind_0] = 1/0.00025 ## specify small r
  C = Istar01
  sigma_zetay01 = diag(rzeta_y01*1/phi2[k-1],m,m) # specify phi as variance
  
  
  var_zetay01 = solve(sigma_zetay01+t(C)%*%C/sigsq01_pos[k])
  mu_zetay01 = var_zetay01%*%(t(C)%*%(data_cl_01[,4]-data_cl_01[,1:3]%*%alpha01_pos[k,])/sigsq01_pos[k])
  zeta_y_01[k,] = mvrnorm(1,mu_zetay01,var_zetay01)
  
  # # (d-2) delta_y01
  # pslab = dnorm(zeta_y_01[k,],0,sqrt(phi2[k-1]))*prob_inc_y01[k-1,]
  # pspike = dnorm(zeta_y_01[k,],0,sqrt(0.00025*phi2[k-1]))*(1-prob_inc_y01[k-1,])
  # delta_y01[k,] = rbinom(m,1,pslab/(pspike+pslab))
  # 
  # # 2-5 : update prob_incs
  # prob_inc_y11[k,] = rbeta(m,1+delta_y11[k,],1+(1-delta_y11[k,]))
  # prob_inc_y10[k,] = rbeta(m,1+delta_y10[k,],1+(1-delta_y10[k,]))
  # prob_inc_y00[k,] = rbeta(m,1+delta_y00[k,],1+(1-delta_y00[k,]))
  # prob_inc_y01[k,] = rbeta(m,1+delta_y01[k,],1+(1-delta_y01[k,]))
  
  # step 3 : update hyperparameters slab-spike variance
  ## 3-1 for prin stra sigma1 
  b1_phi = 0
  for (i in 1:m){
    b1_phi = b1_phi+zeta_00[k,i]^2*rzeta00[i]+zeta_10[k,i]^2*rzeta10[i]+zeta_11[k,i]^2*rzeta11[i]
  }
  phi1[k] = rinvgamma(1,shape=5+3*m/2,rate=50+b1_phi/2)
  
  ## 3-2 for outcome sigma2
  b2_phi = 0
  for (i in 1:m){
    b2_phi = b2_phi+zeta_y_00[k,i]^2*rzeta_y00[i]+zeta_y_10[k,i]^2*rzeta_y10[i]+
      zeta_y_01[k,i]^2*rzeta_y01[i]+zeta_y_11[k,i]^2*rzeta_y11[i]
  }
  phi2[k] = rinvgamma(1,shape=5+4*m/2,rate=50+b2_phi/2)
  
  print(k)
}
  print(dataid)
  return(list(beta00_pos=beta00_pos,beta10_pos=beta10_pos,beta11_pos=beta11_pos,alpha00_pos=alpha00_pos,alpha10_pos=alpha10_pos,alpha01_pos=alpha01_pos,alpha11_pos=alpha11_pos,zeta_11=zeta_11,zeta_10=zeta_10,zeta_00=zeta_00,zeta_y_11=zeta_y_11,zeta_y_10=zeta_y_10,zeta_y_01=zeta_y_01,zeta_y_00=zeta_y_00))
}
poster_nmig_resclass12_mcmc2_single = mclapply(1:4,nmigfunc,iter=iter,burn=burn,mc.cores=20)
save(poster_nmig_resclass12_mcmc2_single,file="poster_nmig_resclass12_mcmc2_single.RData")
