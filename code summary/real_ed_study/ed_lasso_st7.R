load("data_class_as_clu_list.RData")
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)
library(survival)
library(stats)
library(rmutil) # generate laplace
library(statmod) 
library(MatchIt)

iter=11000
burn=1000

##################MCMC#################################
poster12_lasso_func = function(subject,data_list,iter=10000,burn=1000){
  set.seed(subject)
  data_sub = data_list[[subject]]
  ###################### Co_Ed F==1 vc F only==2 ###################################
  data_obs=data_sub[data_sub$group_rcode==1|data_sub$group_rcode==2,]
  data_obs$z = ifelse(data_obs$group_rcode==1,0,1) ## Co_Ed F =0 vc F only=1
  N = dim(data_obs)[1]
  m = length(unique(data_obs$clu))
  n = as.vector(table(data_obs$clu))
  data_obs$clu = rep(1:m,n)
  data_obs$id = unlist(apply(matrix(n),1,seq,from=1,by=1))
  data_obs = data_obs[,c(1:4,8,6,7)]
  colnames(data_obs) = c("clu","id","x1","x2","z","s_obs","y_obs")
  data_obs$x1 = scale(data_obs$x1)
  data_obs$x2 = scale(data_obs$x2)
  data_obs$y_obs = scale(data_obs$y_obs)
  
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
  
  ## hyper lasso parameters
  # for prin strata
  lambdasq1 = rep(-99,iter)
  lasso_sigsq1 = rep(-99,iter)
  # for outcome
  lambdasq2 = rep(-99,iter)
  lasso_sigsq2 = rep(-99,iter)
  
  zeta_11 = matrix(-99,nrow=iter,ncol=m)
  zeta_10 = matrix(-99,nrow=iter,ncol=m)
  zeta_00 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_11 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_10 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_01 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_00 = matrix(-99,nrow=iter,ncol=m)
  
  ## starting value
  # principal strata
  pstra_m = glm(s_obs~-1+x1+x2+as.factor(clu),family=binomial(link = "probit"),data=data_obs)
  beta00_pos[1,] = coef(pstra_m)[1:2]*0.99
  beta10_pos[1,] = coef(pstra_m)[1:2]*0.99
  beta11_pos[1,] = coef(pstra_m)[1:2]*0.99
  var_beta00[1] = 100
  var_beta10[1] = 100
  var_beta11[1] = 100
  
  # outcomes
  out_m = lm(y_obs~-1+x1+x2+z+as.factor(clu),data=data_obs)
  alpha00_pos[1,] = coef(out_m)[1:3]*0.99
  alpha10_pos[1,] = coef(out_m)[1:3]*0.99
  alpha01_pos[1,] = coef(out_m)[1:3]*0.99
  alpha11_pos[1,] = coef(out_m)[1:3]*0.99
  var_alpha11[1] = 100
  var_alpha10[1] = 100
  var_alpha00[1] = 100
  var_alpha01[1] = 100
  
  sigsq00_pos[1] = sum(out_m$residuals^2)/(N-3)
  sigsq01_pos[1] = sum(out_m$residuals^2)/(N-3)
  sigsq10_pos[1] = sum(out_m$residuals^2)/(N-3)
  sigsq11_pos[1] = sum(out_m$residuals^2)/(N-3)
  
  ## hyper lasso parameters
  lambdasq1[1] = rgamma(1,shape=1,rate=0.1)
  lasso_sigsq1[1] = 100
  lambdasq2[1] = rgamma(1,shape=1,rate=0.1)
  lasso_sigsq2[1] = 100
  #zeta_11[1,] = rep(coef(pstra_m)[3],m)
  #zeta_10[1,] = rep(coef(pstra_m)[3],m)
  #zeta_00[1,] = rep(coef(pstra_m)[3],m)
  zeta_11[1,] = rep(min(coef(pstra_m)[3:(m+2)],na.rm=T),m)
  zeta_10[1,] = rep(min(coef(pstra_m)[3:(m+2)],na.rm=T),m)
  zeta_00[1,] = rep(min(coef(pstra_m)[3:(m+2)],na.rm=T),m)
  zeta_y_11[1,] = rep(min(coef(out_m)[4:(m+3)],na.rm=T),m)
  zeta_y_10[1,] = rep(min(coef(out_m)[4:(m+3)],na.rm=T),m)
  zeta_y_01[1,] = rep(min(coef(out_m)[4:(m+3)],na.rm=T),m)
  zeta_y_00[1,] = rep(min(coef(out_m)[4:(m+3)],na.rm=T),m)
  
  ## dummy variable/indicator matrix for cluster
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  for (k in 2:iter){
    # Step 2: update principal strata parameters
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
    
    pi_11[which(pi_11==0)] = 1e-308
    pi_10[which(pi_10==0)] = 1e-308
    pi_01[which(pi_01==0)] = 1e-308
    pi_00[which(pi_00==0)] = 1e-308
    
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
    
    # 2-3 zeta_g, 1,2,...m
    #(a) zeta_11
    tau_11 = rep(-99,m)
    for (cl in 1:m){
      # generate tau_11
      mean_tau11 = sqrt(lambdasq1[k-1]*lasso_sigsq1[k-1]/(zeta_11[k-1,cl]^2))
      invtau11 = rinvgauss(1,mean=mean_tau11,dispersion = 1/lambdasq1[k-1])
      tau_11[cl] = 1/invtau11
      
      # update 
      data_cl_all = data_obs2[data_obs2$clu==cl,]
      n_cl = dim(data_cl_all)[1]
      data_cl = cbind(data_cl_all$x1,data_cl_all$x2) #x1,x2
      cl_Sigma11 = 1/(1/(lasso_sigsq1[k-1]*tau_11[cl])+n_cl)
      cl_mu11 = cl_Sigma11*(sum(data_cl_all$gstar_11-data_cl%*%beta11_pos[k,]))
      zeta_11[k,cl] = rnorm(1,mean=cl_mu11,sd=sqrt(cl_Sigma11))
    }
    #(b) zeta_10
    tau_10 = rep(-99,m)
    for (cl in 1:m){
      # generate tau_10
      mean_tau10 = sqrt(lambdasq1[k-1]*lasso_sigsq1[k-1]/(zeta_10[k-1,cl]^2))
      invtau10 = rinvgauss(1,mean=mean_tau10,dispersion = 1/lambdasq1[k-1])
      tau_10[cl] = 1/invtau10
      # update
      data_cl_all = data_obs2[data_obs2$clu==cl&data_obs2$strata_num!=1,]
      n_cl = dim(data_cl_all)[1]
      data_cl = cbind(data_cl_all$x1,data_cl_all$x2)  #x1,x2
      cl_Sigma10 = 1/(1/(lasso_sigsq1[k-1]*tau_10[cl])+n_cl)
      cl_mu10 = cl_Sigma10*(sum(data_cl_all$gstar_10-data_cl%*%beta10_pos[k,]))
      zeta_10[k,cl] = rnorm(1,mean=cl_mu10,sd=sqrt(cl_Sigma10))
    }
    
    #(c) zeta_00
    tau_00 = rep(-99,m)
    for (cl in 1:m){
      # generate tau_00
      mean_tau00 = sqrt(lambdasq1[k-1]*lasso_sigsq1[k-1]/(zeta_00[k-1,cl]^2))
      invtau00 = rinvgauss(1,mean=mean_tau00,dispersion = 1/lambdasq1[k-1])
      tau_00[cl] = 1/invtau00
      # update
      data_cl_all = data_obs2[data_obs2$clu==cl&data_obs2$strata_num==3|data_obs2$clu==cl&data_obs2$strata_num==4,]
      n_cl = dim(data_cl_all)[1]
      data_cl = cbind(data_cl_all$x1,data_cl_all$x2)  #x1,x2
      cl_Sigma00 = 1/(1/(lasso_sigsq1[k-1]*tau_00[cl])+n_cl)
      cl_mu00 = cl_Sigma00*(sum(data_cl_all$gstar_00-data_cl%*%beta00_pos[k,]))
      zeta_00[k,cl] = rnorm(1,mean=cl_mu00,sd=sqrt(cl_Sigma00))
    }
    
    ## step 3 : update outcome parameters
    # 3-0 : add data information 
    clu_y_11 = Istar%*%zeta_y_11[k-1,]
    clu_y_10 = Istar%*%zeta_y_10[k-1,]
    clu_y_00 = Istar%*%zeta_y_00[k-1,]
    clu_y_01 = Istar%*%zeta_y_01[k-1,]
    data_obs3 =as.data.frame(cbind(data_obs2,clu_y_11,clu_y_10,clu_y_00,clu_y_01))
    
    # 3-1 : alpha_g
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
    
    # 3-4: zeta_yg 1,2,...m
    # (a) zeta_y11
    tau_y11 = rep(-99,m)
    for(cl in 1:m){
      # generate tau_y11
      mean_tauy11 = sqrt(lambdasq2[k-1]*lasso_sigsq2[k-1]/(zeta_y_11[k-1,cl]^2))
      invtauy11 = rinvgauss(1,mean=mean_tauy11,dispersion = 1/lambdasq2[k-1])
      tau_y11[cl] = 1/invtauy11
      # update
      data_cl_11 = as.matrix(data_obs3[data_obs3$clu==cl&data_obs3$strata_num==1,c(3,4,5,7)])
      n_cl_11 = dim(data_cl_11)[1]
      var_zetay11 = 1/(1/(lasso_sigsq2[k-1]*tau_y11[cl])+n_cl_11/sigsq11_pos[k])
      mu_zetay11 = var_zetay11*(rep(1,n_cl_11)%*%(data_cl_11[,4]-data_cl_11[,1:3]%*%alpha11_pos[k,])/sigsq11_pos[k])
      zeta_y_11[k,cl] = rnorm(1,mu_zetay11,sqrt(var_zetay11))
    }
    # (b) zeta_y10
    tau_y10 = rep(-99,m)
    for(cl in 1:m){
      # generate tau_y10
      mean_tauy10 = sqrt(lambdasq2[k-1]*lasso_sigsq2[k-1]/(zeta_y_10[k-1,cl]^2))
      invtauy10 = rinvgauss(1,mean=mean_tauy10,dispersion = 1/lambdasq2[k-1])
      tau_y10[cl] = 1/invtauy10
      # update
      data_cl_10 = as.matrix(data_obs3[data_obs3$clu==cl&data_obs3$strata_num==2,c(3,4,5,7)])
      n_cl_10 = dim(data_cl_10)[1]
      var_zetay10 = 1/(1/(lasso_sigsq2[k-1]*tau_y10[cl])+n_cl_10/sigsq10_pos[k])
      mu_zetay10 = var_zetay10*(rep(1,n_cl_10)%*%(data_cl_10[,4]-data_cl_10[,1:3]%*%alpha10_pos[k,])/sigsq10_pos[k])
      zeta_y_10[k,cl] = rnorm(1,mu_zetay10,sqrt(var_zetay10))
    }
    # (c) zeta_y00
    tau_y00 = rep(-99,m)
    for(cl in 1:m){
      # generate tau_y00
      mean_tauy00 = sqrt(lambdasq2[k-1]*lasso_sigsq2[k-1]/(zeta_y_00[k-1,cl]^2))
      invtauy00 = rinvgauss(1,mean=mean_tauy00,dispersion = 1/lambdasq2[k-1])
      tau_y00[cl] = 1/invtauy00
      # update
      data_cl_00 = as.matrix(data_obs3[data_obs3$clu==cl&data_obs3$strata_num==3,c(3,4,5,7)])
      n_cl_00 = dim(data_cl_00)[1]
      var_zetay00 = 1/(1/(lasso_sigsq2[k-1]*tau_y00[cl])+n_cl_00/sigsq00_pos[k])
      mu_zetay00 = var_zetay00*(rep(1,n_cl_00)%*%(data_cl_00[,4]-data_cl_00[,1:3]%*%alpha00_pos[k,])/sigsq00_pos[k])
      zeta_y_00[k,cl] = rnorm(1,mu_zetay00,sqrt(var_zetay00))
    }
    # (d) zeta_y01
    tau_y01 = rep(-99,m)
    for(cl in 1:m){
      # generate tau_y01
      mean_tauy01 = sqrt(lambdasq2[k-1]*lasso_sigsq2[k-1]/(zeta_y_01[k-1,cl]^2))
      invtauy01 = rinvgauss(1,mean=mean_tauy01,dispersion = 1/lambdasq2[k-1])
      tau_y01[cl] = 1/invtauy01
      # update
      data_cl_01 = as.matrix(data_obs3[data_obs3$clu==cl&data_obs3$strata_num==4,c(3,4,5,7)])
      n_cl_01 = dim(data_cl_01)[1]
      var_zetay01 = 1/(1/(lasso_sigsq2[k-1]*tau_y01[cl])+n_cl_01/sigsq01_pos[k])
      mu_zetay01 = var_zetay01*(rep(1,n_cl_01)%*%(data_cl_01[,4]-data_cl_01[,1:3]%*%alpha01_pos[k,])/sigsq01_pos[k])
      zeta_y_01[k,cl] = rnorm(1,mu_zetay01,sqrt(var_zetay01))
    }
    # step 4 : update hyperparameters
    # 4-1 lambdasq1 and 2
    lambda_rate1 = sum(tau_00)+sum(tau_10)+sum(tau_11)
    lambdasq1[k] = rgamma(1,shape=3*m+1,rate=0.1+lambda_rate1/2)
    
    lambda_rate2 = sum(tau_y00)+sum(tau_y10)+sum(tau_y11)+sum(tau_y01)
    lambdasq2[k] = rgamma(1,shape=4*m+1,rate=0.1+lambda_rate2/2)
    
    # 4-2 lasso_sigsq1 and 2
    lasso_rate1 =0
    lasso_rate2 =0
    for (i in 1:m){
      lasso_rate1 = lasso_rate1+zeta_00[k,i]^2/tau_00[i]+zeta_10[k,i]^2/tau_10[i]+zeta_11[k,i]^2/tau_11[i]
      lasso_rate2 = lasso_rate2+zeta_y_00[k,i]^2/tau_y00[i]+zeta_y_01[k,i]^2/tau_y01[i]+zeta_y_10[k,i]^2/tau_y10[i]+zeta_y_11[k,i]^2/tau_y11[i]
    }
    lasso_sigsq1[k] = rinvgamma(1,shape=3*m/2,rate=lasso_rate1/2)
    lasso_sigsq2[k] = rinvgamma(1,shape=4*m/2,rate=lasso_rate2/2)
    
    #print(k)
  }
  print(subject)
  return(list(data_obs=data_obs,beta00_pos=beta00_pos,beta10_pos=beta10_pos,beta11_pos=beta11_pos,var_beta00=var_beta00,
              var_beta10=var_beta10,var_beta11=var_beta11,alpha00_pos=alpha00_pos,alpha10_pos=alpha10_pos,alpha01_pos=alpha01_pos,
              alpha11_pos=alpha11_pos,var_alpha11=var_alpha11,var_alpha10=var_alpha10,var_alpha00=var_alpha00,var_alpha01=var_alpha01,
              sigsq00_pos=sigsq00_pos,sigsq01_pos=sigsq01_pos,sigsq10_pos=sigsq10_pos,sigsq11_pos=sigsq11_pos,lambdasq1=lambdasq1,lambdasq2=lambdasq2,
              lasso_sigsq1=lasso_sigsq1,lasso_sigsq2=lasso_sigsq2,
              zeta_11=zeta_11,zeta_10=zeta_10,zeta_00=zeta_00,zeta_y_11=zeta_y_11,zeta_y_10=zeta_y_10,zeta_y_01=zeta_y_01,zeta_y_00=zeta_y_00))
  
  
  
}


poster34_lasso_func = function(subject,data_list,iter=11000,burn=1000){
  set.seed(subject)
  data_sub = data_list[[subject]]
  ###################### Co_Ed M ==3 vc M only==4 ###################################
  data_obs=data_sub[data_sub$group_rcode==3|data_sub$group_rcode==4,]
  data_obs$z = ifelse(data_obs$group_rcode==3,0,1) ## Co_Ed M =0 vc M only=1
  N = dim(data_obs)[1]
  m = length(unique(data_obs$clu))
  n = as.vector(table(data_obs$clu))
  data_obs$clu = rep(1:m,n)
  data_obs$id = unlist(apply(matrix(n),1,seq,from=1,by=1))
  data_obs = data_obs[,c(1:4,8,6,7)]
  colnames(data_obs) = c("clu","id","x1","x2","z","s_obs","y_obs")
  data_obs$x1 = scale(data_obs$x1)
  data_obs$x2 = scale(data_obs$x2)
  data_obs$y_obs = scale(data_obs$y_obs)
  
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
  
  ## hyper lasso parameters
  # for prin strata
  lambdasq1 = rep(-99,iter)
  lasso_sigsq1 = rep(-99,iter)
  # for outcome
  lambdasq2 = rep(-99,iter)
  lasso_sigsq2 = rep(-99,iter)
  
  zeta_11 = matrix(-99,nrow=iter,ncol=m)
  zeta_10 = matrix(-99,nrow=iter,ncol=m)
  zeta_00 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_11 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_10 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_01 = matrix(-99,nrow=iter,ncol=m)
  zeta_y_00 = matrix(-99,nrow=iter,ncol=m)
  
  ## starting value
  # principal strata
  pstra_m = glm(s_obs~-1+x1+x2+as.factor(clu),family=binomial(link = "probit"),data=data_obs)
  #beta00_pos[1,] = coef(pstra_m)[1:2]
  #beta10_pos[1,] = coef(pstra_m)[1:2]
  #beta11_pos[1,] = coef(pstra_m)[1:2]
  beta00_pos[1,] = coef(pstra_m)[1:2]*0.99
  beta10_pos[1,] = coef(pstra_m)[1:2]*0.99
  beta11_pos[1,] = coef(pstra_m)[1:2]*0.99
  var_beta00[1] = 100
  var_beta10[1] = 100
  var_beta11[1] = 100
  
  # outcomes
  out_m = lm(y_obs~-1+x1+x2+z+as.factor(clu),data=data_obs)
  alpha00_pos[1,] = coef(out_m)[1:3]*0.99
  alpha10_pos[1,] = coef(out_m)[1:3]*0.99
  alpha01_pos[1,] = coef(out_m)[1:3]*0.99
  alpha11_pos[1,] = coef(out_m)[1:3]*0.99
  var_alpha11[1] = 100
  var_alpha10[1] = 100
  var_alpha00[1] = 100
  var_alpha01[1] = 100
  
  sigsq00_pos[1] = sum(out_m$residuals^2)/(N-3)
  sigsq01_pos[1] = sum(out_m$residuals^2)/(N-3)
  sigsq10_pos[1] = sum(out_m$residuals^2)/(N-3)
  sigsq11_pos[1] = sum(out_m$residuals^2)/(N-3)
  
  ## hyper lasso parameters
  lambdasq1[1] = rgamma(1,shape=1,rate=0.1)
  lasso_sigsq1[1] = 100
  lambdasq2[1] = rgamma(1,shape=1,rate=0.1)
  lasso_sigsq2[1] = 100
  #zeta_11[1,] = rep(coef(pstra_m)[5],m)
  #zeta_10[1,] = rep(coef(pstra_m)[5],m)
  #zeta_00[1,] = rep(coef(pstra_m)[5],m)
  zeta_11[1,] = rep(min(coef(pstra_m)[3:(m+2)],na.rm=T),m)
  zeta_10[1,] = rep(min(coef(pstra_m)[3:(m+2)],na.rm=T),m)
  zeta_00[1,] = rep(min(coef(pstra_m)[3:(m+2)],na.rm=T),m)
  zeta_y_11[1,] = rep(min(coef(out_m)[4:(m+3)],na.rm=T),m)
  zeta_y_10[1,] = rep(min(coef(out_m)[4:(m+3)],na.rm=T),m)
  zeta_y_01[1,] = rep(min(coef(out_m)[4:(m+3)],na.rm=T),m)
  zeta_y_00[1,] = rep(min(coef(out_m)[4:(m+3)],na.rm=T),m)
  
  ## dummy variable/indicator matrix for cluster
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  for (k in 2:iter){
    # Step 2: update principal strata parameters
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
    
    pi_11[which(pi_11==0)] = 1e-308
    pi_10[which(pi_10==0)] = 1e-308
    pi_01[which(pi_01==0)] = 1e-308
    pi_00[which(pi_00==0)] = 1e-308
    
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
    
    # 2-3 zeta_g, 1,2,...m
    #(a) zeta_11
    tau_11 = rep(-99,m)
    for (cl in 1:m){
      # generate tau_11
      mean_tau11 = sqrt(lambdasq1[k-1]*lasso_sigsq1[k-1]/(zeta_11[k-1,cl]^2))
      invtau11 = rinvgauss(1,mean=mean_tau11,dispersion = 1/lambdasq1[k-1])
      tau_11[cl] = 1/invtau11
      
      # update 
      data_cl_all = data_obs2[data_obs2$clu==cl,]
      n_cl = dim(data_cl_all)[1]
      data_cl = cbind(data_cl_all$x1,data_cl_all$x2) #x1,x2
      cl_Sigma11 = 1/(1/(lasso_sigsq1[k-1]*tau_11[cl])+n_cl)
      cl_mu11 = cl_Sigma11*(sum(data_cl_all$gstar_11-data_cl%*%beta11_pos[k,]))
      zeta_11[k,cl] = rnorm(1,mean=cl_mu11,sd=sqrt(cl_Sigma11))
    }
    #(b) zeta_10
    tau_10 = rep(-99,m)
    for (cl in 1:m){
      # generate tau_10
      mean_tau10 = sqrt(lambdasq1[k-1]*lasso_sigsq1[k-1]/(zeta_10[k-1,cl]^2))
      invtau10 = rinvgauss(1,mean=mean_tau10,dispersion = 1/lambdasq1[k-1])
      tau_10[cl] = 1/invtau10
      # update
      data_cl_all = data_obs2[data_obs2$clu==cl&data_obs2$strata_num!=1,]
      n_cl = dim(data_cl_all)[1]
      data_cl = cbind(data_cl_all$x1,data_cl_all$x2)  #x1,x2
      cl_Sigma10 = 1/(1/(lasso_sigsq1[k-1]*tau_10[cl])+n_cl)
      cl_mu10 = cl_Sigma10*(sum(data_cl_all$gstar_10-data_cl%*%beta10_pos[k,]))
      zeta_10[k,cl] = rnorm(1,mean=cl_mu10,sd=sqrt(cl_Sigma10))
    }
    
    #(c) zeta_00
    tau_00 = rep(-99,m)
    for (cl in 1:m){
      # generate tau_00
      mean_tau00 = sqrt(lambdasq1[k-1]*lasso_sigsq1[k-1]/(zeta_00[k-1,cl]^2))
      invtau00 = rinvgauss(1,mean=mean_tau00,dispersion = 1/lambdasq1[k-1])
      tau_00[cl] = 1/invtau00
      # update
      data_cl_all = data_obs2[data_obs2$clu==cl&data_obs2$strata_num==3|data_obs2$clu==cl&data_obs2$strata_num==4,]
      n_cl = dim(data_cl_all)[1]
      data_cl = cbind(data_cl_all$x1,data_cl_all$x2)  #x1,x2
      cl_Sigma00 = 1/(1/(lasso_sigsq1[k-1]*tau_00[cl])+n_cl)
      cl_mu00 = cl_Sigma00*(sum(data_cl_all$gstar_00-data_cl%*%beta00_pos[k,]))
      zeta_00[k,cl] = rnorm(1,mean=cl_mu00,sd=sqrt(cl_Sigma00))
    }
    
    ## step 3 : update outcome parameters
    # 3-0 : add data information 
    clu_y_11 = Istar%*%zeta_y_11[k-1,]
    clu_y_10 = Istar%*%zeta_y_10[k-1,]
    clu_y_00 = Istar%*%zeta_y_00[k-1,]
    clu_y_01 = Istar%*%zeta_y_01[k-1,]
    data_obs3 =as.data.frame(cbind(data_obs2,clu_y_11,clu_y_10,clu_y_00,clu_y_01))
    
    # 3-1 : alpha_g
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
    
    # 3-4: zeta_yg 1,2,...m
    # (a) zeta_y11
    tau_y11 = rep(-99,m)
    for(cl in 1:m){
      # generate tau_y11
      mean_tauy11 = sqrt(lambdasq2[k-1]*lasso_sigsq2[k-1]/(zeta_y_11[k-1,cl]^2))
      invtauy11 = rinvgauss(1,mean=mean_tauy11,dispersion = 1/lambdasq2[k-1])
      tau_y11[cl] = 1/invtauy11
      # update
      data_cl_11 = as.matrix(data_obs3[data_obs3$clu==cl&data_obs3$strata_num==1,c(3,4,5,7)])
      n_cl_11 = dim(data_cl_11)[1]
      var_zetay11 = 1/(1/(lasso_sigsq2[k-1]*tau_y11[cl])+n_cl_11/sigsq11_pos[k])
      mu_zetay11 = var_zetay11*(rep(1,n_cl_11)%*%(data_cl_11[,4]-data_cl_11[,1:3]%*%alpha11_pos[k,])/sigsq11_pos[k])
      zeta_y_11[k,cl] = rnorm(1,mu_zetay11,sqrt(var_zetay11))
    }
    # (b) zeta_y10
    tau_y10 = rep(-99,m)
    for(cl in 1:m){
      # generate tau_y10
      mean_tauy10 = sqrt(lambdasq2[k-1]*lasso_sigsq2[k-1]/(zeta_y_10[k-1,cl]^2))
      invtauy10 = rinvgauss(1,mean=mean_tauy10,dispersion = 1/lambdasq2[k-1])
      tau_y10[cl] = 1/invtauy10
      # update
      data_cl_10 = as.matrix(data_obs3[data_obs3$clu==cl&data_obs3$strata_num==2,c(3,4,5,7)])
      n_cl_10 = dim(data_cl_10)[1]
      var_zetay10 = 1/(1/(lasso_sigsq2[k-1]*tau_y10[cl])+n_cl_10/sigsq10_pos[k])
      mu_zetay10 = var_zetay10*(rep(1,n_cl_10)%*%(data_cl_10[,4]-data_cl_10[,1:3]%*%alpha10_pos[k,])/sigsq10_pos[k])
      zeta_y_10[k,cl] = rnorm(1,mu_zetay10,sqrt(var_zetay10))
    }
    # (c) zeta_y00
    tau_y00 = rep(-99,m)
    for(cl in 1:m){
      # generate tau_y00
      mean_tauy00 = sqrt(lambdasq2[k-1]*lasso_sigsq2[k-1]/(zeta_y_00[k-1,cl]^2))
      invtauy00 = rinvgauss(1,mean=mean_tauy00,dispersion = 1/lambdasq2[k-1])
      tau_y00[cl] = 1/invtauy00
      # update
      data_cl_00 = as.matrix(data_obs3[data_obs3$clu==cl&data_obs3$strata_num==3,c(3,4,5,7)])
      n_cl_00 = dim(data_cl_00)[1]
      var_zetay00 = 1/(1/(lasso_sigsq2[k-1]*tau_y00[cl])+n_cl_00/sigsq00_pos[k])
      mu_zetay00 = var_zetay00*(rep(1,n_cl_00)%*%(data_cl_00[,4]-data_cl_00[,1:3]%*%alpha00_pos[k,])/sigsq00_pos[k])
      zeta_y_00[k,cl] = rnorm(1,mu_zetay00,sqrt(var_zetay00))
    }
    # (d) zeta_y01
    tau_y01 = rep(-99,m)
    for(cl in 1:m){
      # generate tau_y01
      mean_tauy01 = sqrt(lambdasq2[k-1]*lasso_sigsq2[k-1]/(zeta_y_01[k-1,cl]^2))
      invtauy01 = rinvgauss(1,mean=mean_tauy01,dispersion = 1/lambdasq2[k-1])
      tau_y01[cl] = 1/invtauy01
      # update
      data_cl_01 = as.matrix(data_obs3[data_obs3$clu==cl&data_obs3$strata_num==4,c(3,4,5,7)])
      n_cl_01 = dim(data_cl_01)[1]
      var_zetay01 = 1/(1/(lasso_sigsq2[k-1]*tau_y01[cl])+n_cl_01/sigsq01_pos[k])
      mu_zetay01 = var_zetay01*(rep(1,n_cl_01)%*%(data_cl_01[,4]-data_cl_01[,1:3]%*%alpha01_pos[k,])/sigsq01_pos[k])
      zeta_y_01[k,cl] = rnorm(1,mu_zetay01,sqrt(var_zetay01))
    }
    # step 4 : update hyperparameters
    # 4-1 lambdasq1 and 2
    lambda_rate1 = sum(tau_00)+sum(tau_10)+sum(tau_11)
    lambdasq1[k] = rgamma(1,shape=3*m+1,rate=0.1+lambda_rate1/2)
    
    lambda_rate2 = sum(tau_y00)+sum(tau_y10)+sum(tau_y11)+sum(tau_y01)
    lambdasq2[k] = rgamma(1,shape=4*m+1,rate=0.1+lambda_rate2/2)
    
    # 4-2 lasso_sigsq1 and 2
    lasso_rate1 =0
    lasso_rate2 =0
    for (i in 1:m){
      lasso_rate1 = lasso_rate1+zeta_00[k,i]^2/tau_00[i]+zeta_10[k,i]^2/tau_10[i]+zeta_11[k,i]^2/tau_11[i]
      lasso_rate2 = lasso_rate2+zeta_y_00[k,i]^2/tau_y00[i]+zeta_y_01[k,i]^2/tau_y01[i]+zeta_y_10[k,i]^2/tau_y10[i]+zeta_y_11[k,i]^2/tau_y11[i]
    }
    lasso_sigsq1[k] = rinvgamma(1,shape=3*m/2,rate=lasso_rate1/2)
    lasso_sigsq2[k] = rinvgamma(1,shape=4*m/2,rate=lasso_rate2/2)
    
    #print(k)
  }
  print(subject)
  return(list(data_obs=data_obs,beta00_pos=beta00_pos,beta10_pos=beta10_pos,beta11_pos=beta11_pos,var_beta00=var_beta00,
              var_beta10=var_beta10,var_beta11=var_beta11,alpha00_pos=alpha00_pos,alpha10_pos=alpha10_pos,alpha01_pos=alpha01_pos,
              alpha11_pos=alpha11_pos,var_alpha11=var_alpha11,var_alpha10=var_alpha10,var_alpha00=var_alpha00,var_alpha01=var_alpha01,
              sigsq00_pos=sigsq00_pos,sigsq01_pos=sigsq01_pos,sigsq10_pos=sigsq10_pos,sigsq11_pos=sigsq11_pos,lambdasq1=lambdasq1,lambdasq2=lambdasq2,
              lasso_sigsq1=lasso_sigsq1,lasso_sigsq2=lasso_sigsq2,
              zeta_11=zeta_11,zeta_10=zeta_10,zeta_00=zeta_00,zeta_y_11=zeta_y_11,zeta_y_10=zeta_y_10,zeta_y_01=zeta_y_01,zeta_y_00=zeta_y_00))
  
  
  
}
lasso12_class_posterior_st7 = mclapply(1:4,poster12_lasso_func,data_class_as_clu_list,iter=iter,burn=burn,mc.cores=4)
save(lasso12_class_posterior_st7,file="lasso12_class_posterior_st7.RData")

lasso34_class_posterior_st7 = mclapply(1:4,poster34_lasso_func,data_class_as_clu_list,iter=iter,burn=burn,mc.cores=4)
save(lasso34_class_posterior_st7,file="lasso34_class_posterior_st7.RData")

