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

data_gene = function(dataid,m=50,m0=30){
  set.seed(dataid)
  #n = as.matrix(sample(6:15,m,replace=T,prob=rep(1/10,10)))
  #id = unlist(apply(n,1,seq,from=1,by=1))
  n = as.matrix(sample(rep(10,10),m,replace=T,prob=rep(1/10,10)))
  id = rep(c(1,2),m)
  clu = rep(1:m,n)
  data = as.data.frame(cbind(clu,id))
  N = dim(data)[1]
  #data$x1 = rnorm(N,0,1)
  #data$x2 = sample(c(-1,0,1),N,replace=T,prob=rep(1/3,3))
  X = rmvnorm(N,mean=rep(0,2),sigma = matrix(c(1,0.5,0.5,1),nrow=2,ncol=2))
  data$x1 = X[,1]
  data$x2 = X[,2]
  ui = rnorm(m,3,1)
  ui_0 = rnorm(m0,0,0.1) 
  ui[1:m0] =ui_0 ### some cluster effect close to 0
  data$U = rep(ui,n)
  data$z = rbinom(N,1,0.5)
  # principal strata
  ######## where x includes x1,x2,U
  prinstra_data = as.matrix(cbind(data$x1,data$x2,data$U))
  gstar_11 = prinstra_data%*%beta11_t+rnorm(N,0,1)
  gstar_10 = prinstra_data%*%beta10_t+rnorm(N,0,1)
  gstar_00 = prinstra_data%*%beta00_t+rnorm(N,0,1)
  G = ifelse(gstar_11<=0,1,ifelse(gstar_10<=0,2,ifelse(gstar_00<=0,3,4))) # principal membership
  data$G = G
  
  data$s0 = rep(-99,N)
  data$s1 = rep(-99,N)
  for (i in 1:N){
    if (data$G[i]==1){
      data$s0[i] = 1
      data$s1[i] = 1
    }else if (data$G[i]==2){
      data$s0[i] = 1
      data$s1[i] = 0
    }else if (data$G[i]==3){
      data$s0[i] = 0
      data$s1[i] = 0
    }else if (data$G[i]==4){
      data$s0[i] = 0
      data$s1[i] = 1
    }
  }
  data$s_obs = data$z*data$s1+(1-data$z)*data$s0
  
  # outcomes
  ########################
  y1=rep(-99,N)
  y0=rep(-99,N)
  data_y_0 = cbind(data$x1,data$x2,rep(0,N),data$U)
  data_y_1 = cbind(data$x1,data$x2,rep(1,N),data$U)
  for (i in 1:N){
    if (data$G[i]==3){
      y0[i] = rnorm(1,mean = data_y_0[i,]%*%alpha00_t,sd = sqrt(sigsq00_t))
      y1[i] = rnorm(1,mean = data_y_1[i,]%*%alpha00_t,sd = sqrt(sigsq00_t))
    }else if (data$G[i]==2){
      y0[i] = rnorm(1,mean = data_y_0[i,]%*%alpha10_t,sd = sqrt(sigsq10_t))
      y1[i] = rnorm(1,mean = data_y_1[i,]%*%alpha10_t,sd = sqrt(sigsq10_t))
    }else if (data$G[i]==4){
      y0[i] = rnorm(1,mean = data_y_0[i,]%*%alpha01_t,sd = sqrt(sigsq01_t))
      y1[i] = rnorm(1,mean = data_y_1[i,]%*%alpha01_t,sd = sqrt(sigsq01_t))
    }else if (data$G[i]==1){
      y0[i] = rnorm(1,mean = data_y_0[i,]%*%alpha11_t,sd = sqrt(sigsq11_t))
      y1[i] = rnorm(1,mean = data_y_1[i,]%*%alpha11_t,sd = sqrt(sigsq11_t))
    }
  }
  y_obs = data$z*y1 + (1-data$z)*y0
  data = as.data.frame(cbind(data,y0,y1,y_obs))
  effect_simu = mean(data$y1-data$y0)
  ##### if change the data structure, double check the column numbers
  data_obs = data[,c(1:4,6,10,13)] # include int
  return (list(data=data,data_obs=data_obs,m=m,n=n,ui=ui,effect_simu=effect_simu))
}
simdata_m50= mclapply(1:100,data_gene,m=50,m0=30,mc.cores=20)
save(simdata_m50,file="simdata_m50.RData")