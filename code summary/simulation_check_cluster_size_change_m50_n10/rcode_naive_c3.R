load("simdata_m50.RData")
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

###########################################################################
######################### MCMC gibb sample ################################
###########################################################################

naivefunc = function(dataid,simdata,iter=100,burn=0){
  set.seed(dataid)
  data_obs = simdata[[dataid]]$data_obs
  N = dim(data_obs)[1]
  m = simdata[[dataid]]$m
  n = simdata[[dataid]]$n
  ui = simdata[[dataid]]$ui
  
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
  
  ## starting value
  # v2:+0.1, v3:-1
  # principal strata
  beta00_pos[1,] = beta00_t[1:2]-1
  beta10_pos[1,] = beta10_t[1:2]-1
  beta11_pos[1,] = beta11_t[1:2]-1
  var_beta00[1] = 100-1
  var_beta10[1] = 100-1
  var_beta11[1] = 100-1
  
  # outcomes
  alpha00_pos[1,] = alpha00_t[1:3]-1
  alpha10_pos[1,] = alpha10_t[1:3]-1
  alpha01_pos[1,] = alpha01_t[1:3]-1
  alpha11_pos[1,] = alpha11_t[1:3]-1
  var_alpha11[1] = 100-1
  var_alpha10[1] = 100-1
  var_alpha00[1] = 100-1
  var_alpha01[1] = 100-1
  
  sigsq00_pos[1] = sigsq00_t-1
  sigsq01_pos[1] = sigsq01_t-1
  sigsq10_pos[1] = sigsq10_t-1
  sigsq11_pos[1] = sigsq11_t-1
  
  for (k in 2:iter){
    # Step 2: update principal strata parameters
    # 2-0 : sample G
    data_G = cbind(data_obs$x1,data_obs$x2)
    mean_11 = data_G%*%beta11_pos[k-1,]
    mean_10 = data_G%*%beta10_pos[k-1,]
    mean_00 = data_G%*%beta00_pos[k-1,]
    p_11 = apply(mean_11,1,pnorm,mean=0,sd=1)
    p_10 = apply(mean_10,1,pnorm,mean=0,sd=1)
    p_00 = apply(mean_00,1,pnorm,mean=0,sd=1)
    
    pi_11 = 1-p_11
    pi_10 = p_11*(1-p_10)
    pi_00 = p_11*p_10*(1-p_00)
    pi_01 = p_11*p_10*p_00
    
    fi_func = function(data_yg,alpha_pos,sigsq){
      # data_yg includes y_obs,x1,x2
      yobs = data_yg[1]
      data_stra = data_yg[2:3]
      fi0 = dnorm(yobs,mean=data_stra%*%alpha_pos[1:2],sd=sqrt(sigsq))
      fi1 = dnorm(yobs,mean=data_stra%*%alpha_pos[1:2]+alpha_pos[3],sd=sqrt(sigsq))
      return(c(fi0,fi1))
    }
    datay = as.matrix(cbind(data_obs$y_obs,data_obs$x1,data_obs$x2))
    
    fi_0_00 = apply(datay,1,fi_func,alpha00_pos[k-1,],sigsq00_pos[k-1])[1,] 
    fi_0_01 = apply(datay,1,fi_func,alpha01_pos[k-1,],sigsq01_pos[k-1])[1,] 
    fi_0_10 = apply(datay,1,fi_func,alpha10_pos[k-1,],sigsq10_pos[k-1])[1,] 
    fi_0_11 = apply(datay,1,fi_func,alpha11_pos[k-1,],sigsq11_pos[k-1])[1,] 
    
    fi_1_00 = apply(datay,1,fi_func,alpha00_pos[k-1,],sigsq00_pos[k-1])[2,] 
    fi_1_01 = apply(datay,1,fi_func,alpha01_pos[k-1,],sigsq01_pos[k-1])[2,] 
    fi_1_10 = apply(datay,1,fi_func,alpha10_pos[k-1,],sigsq10_pos[k-1])[2,] 
    fi_1_11 = apply(datay,1,fi_func,alpha11_pos[k-1,],sigsq11_pos[k-1])[2,] 
    
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
    
    data_obs2 = as.data.frame(cbind(data_obs,strata_num))
    
    # 2-1: sample latent variable G*
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
    
    # 2-2: update beta_11,beta_10,beta_00
    #(a) beta_11
    data_pristr = cbind(data_obs2$x1,data_obs2$x2)
    beta11_Sigma = solve(1/var_beta11[k-1]*diag(2)+t(data_pristr)%*%data_pristr)
    beta11_mu = beta11_Sigma%*%(t(data_pristr)%*%(data_obs2$gstar_11))
    beta11_pos[k,] = mvrnorm(1,beta11_mu,beta11_Sigma)
    
    #(b) beta_10
    datano11 = data_obs2[data_obs2$strata_num!=1,]
    data_pristr = cbind(datano11$x1,datano11$x2)
    beta10_Sigma = solve(1/var_beta10[k-1]*diag(2)+t(data_pristr)%*%data_pristr)
    beta10_mu = beta10_Sigma%*%(t(data_pristr)%*%(datano11$gstar_10))
    beta10_pos[k,] = mvrnorm(1,beta10_mu,beta10_Sigma)
    
    #(c) beta_00
    datano1110 = datano11[datano11$strata_num!=2,]
    data_pristr = cbind(datano1110$x1,datano1110$x2)
    beta00_Sigma = solve(1/var_beta00[k-1]*diag(2)+t(data_pristr)%*%data_pristr)
    beta00_mu = beta00_Sigma%*%(t(data_pristr)%*%(datano1110$gstar_00))
    beta00_pos[k,] = mvrnorm(1,beta00_mu,beta00_Sigma)
    
    # 2-3 var_beta
    var_beta00[k] = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(beta00_pos[k,])%*%beta00_pos[k,]/2)
    var_beta10[k] = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(beta10_pos[k,])%*%beta10_pos[k,]/2)
    var_beta11[k] = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(beta11_pos[k,])%*%beta11_pos[k,]/2)
    
    data_obs3 = data_obs2
    
    ## step 3 : update outcome parameters
    # 3-1 : alpha_g
    update_alpha_func = function(data,sigsq_alpha,sig){
      if (dim(data)[1]==1){
        var_alpha = solve(1/sigsq_alpha*diag(3)+matrix(data[,1:3])%*%data[,1:3]/sig)
        mu_alpha = var_alpha%*%(matrix(data[,1:3])%*%data[,4]/sig)
        return (mvrnorm(1,mu=mu_alpha,Sigma=var_alpha))
      } else if (dim(data)[1]>1){
        var_alpha = solve(1/sigsq_alpha*diag(3)+t(data[,1:3])%*%data[,1:3]/sig)
        mu_alpha = var_alpha%*%(t(data[,1:3])%*%data[,4]/sig)
        return (mvrnorm(1,mu=mu_alpha,Sigma=var_alpha))
      }
    }
    # (a) : alpha_11
    data_11 = as.matrix(data_obs3[data_obs3$strata_num==1,c(3,4,5,7)]) #x1,x2,z,y_obs
    if (dim(data_11)[1]==0){
      alpha11_pos[k,] = alpha11_pos[k-1,]
    }else{
      alpha11_pos[k,] = update_alpha_func(data_11,var_alpha11[k-1],sigsq11_pos[k-1])
    }
    # (b) : alpha_10
    data_10 = as.matrix(data_obs3[data_obs3$strata_num==2,c(3,4,5,7)]) #x1,x2,z,y_obs
    if (dim(data_10)[1]==0){
      alpha10_pos[k,] = alpha10_pos[k-1,]
    }else{
      alpha10_pos[k,] = update_alpha_func(data_10,var_alpha10[k-1],sigsq10_pos[k-1])
    }
    # (c) : alpha_00
    data_00 = as.matrix(data_obs3[data_obs3$strata_num==3,c(3,4,5,7)]) #x1,x2,z,y_obs
    if (dim(data_00)[1]==0){
      alpha00_pos[k,] = alpha00_pos[k-1,]
    }else{
      alpha00_pos[k,] = update_alpha_func(data_00,var_alpha00[k-1],sigsq00_pos[k-1])
    }
    # (d) : alpha_01
    data_01 = as.matrix(data_obs3[data_obs3$strata_num==4,c(3,4,5,7)]) #x1,x2,z,y_obs
    if (dim(data_01)[1]==0){
      alpha01_pos[k,] = alpha01_pos[k-1,]
    }else{
      alpha01_pos[k,] = update_alpha_func(data_01,var_alpha01[k-1],sigsq01_pos[k-1])
    }
    
    # 3-2: var_alpha_g
    var_alpha11[k] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(alpha11_pos[k,])%*%alpha11_pos[k,]/2)
    var_alpha10[k] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(alpha10_pos[k,])%*%alpha10_pos[k,]/2)
    var_alpha01[k] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(alpha01_pos[k,])%*%alpha01_pos[k,]/2)
    var_alpha00[k] = rinvgamma(1,shape=0.001+4/2,rate=0.001+t(alpha00_pos[k,])%*%alpha00_pos[k,]/2)
    
    # 3-3: sigsq_g
    ry11 = data_11[,4]-data_11[,1:3]%*%alpha11_pos[k,]
    shape_sigsq11 = (0.002+nrow(ry11))/2
    scale_sigsq11 = (0.002*1 + t(ry11)%*%ry11)/2
    sigsq11_pos[k] = rinvgamma(1,shape=shape_sigsq11,rate=scale_sigsq11)
    
    ry10 = data_10[,4]-data_10[,1:3]%*%alpha10_pos[k,]
    shape_sigsq10 = (0.002+nrow(ry10))/2
    scale_sigsq10 = (0.002*1 + t(ry10)%*%ry10)/2
    sigsq10_pos[k] = rinvgamma(1,shape=shape_sigsq10,rate=scale_sigsq10)
    
    ry00 = data_00[,4]-data_00[,1:3]%*%alpha00_pos[k,]
    shape_sigsq00 = (0.002+nrow(ry00))/2
    scale_sigsq00 = (0.002*1 + t(ry00)%*%ry00)/2
    sigsq00_pos[k] = rinvgamma(1,shape=shape_sigsq00,rate=scale_sigsq00)
    
    ry01 = data_01[,4]-data_01[,1:3]%*%alpha01_pos[k,]
    shape_sigsq01 = (0.002+nrow(ry01))/2
    scale_sigsq01 = (0.002*1 + t(ry01)%*%ry01)/2
    sigsq01_pos[k] = rinvgamma(1,shape=shape_sigsq01,rate=scale_sigsq01)
    
    #print(k)
    
  }
  print(dataid)
  return(list(beta00_pos=beta00_pos,beta10_pos=beta10_pos,beta11_pos=beta11_pos,alpha00_pos=alpha00_pos,alpha10_pos=alpha10_pos,alpha01_pos=alpha01_pos,alpha11_pos=alpha11_pos,var_alpha11=var_alpha11,sigsq00_pos=sigsq00_pos,sigsq01_pos=sigsq01_pos,sigsq10_pos=sigsq10_pos,sigsq11_pos=sigsq11_pos))
}
poster_naive_resm50_st3 = mclapply(1:100,naivefunc,simdata_m50,iter=iter,burn=burn,mc.cores = 20)
save(poster_naive_resm50_st3,file="poster_naive_resm50_st3.RData")
