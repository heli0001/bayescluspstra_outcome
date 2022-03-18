library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

load("data_class_as_clu_list.RData")
load("nmig12_class_posterior.RData")
load("poster_nmig_resclass12_mcmc2.RData")


iter=11000
burn=1000
thin=10
## recall posterior results
#################### for principal strata 
beta00_pos = poster_nmig_resclass12_mcmc2[[2]]$beta00_pos[1001:11000,]
beta10_pos = poster_nmig_resclass12_mcmc2[[2]]$beta10_pos[1001:11000,]
beta11_pos = poster_nmig_resclass12_mcmc2[[2]]$beta11_pos[1001:11000,]                          

#################### for outcomes  
alpha00_pos = poster_nmig_resclass12_mcmc2[[2]]$alpha00_pos[1001:11000,]
alpha01_pos = poster_nmig_resclass12_mcmc2[[2]]$alpha01_pos[1001:11000,]
alpha11_pos = poster_nmig_resclass12_mcmc2[[2]]$alpha11_pos[1001:11000,] 
alpha10_pos = poster_nmig_resclass12_mcmc2[[2]]$alpha10_pos[1001:11000,]

################# cluster effects
zeta_00 = poster_nmig_resclass12_mcmc2[[2]]$zeta_00[1001:11000,]
zeta_10 = poster_nmig_resclass12_mcmc2[[2]]$zeta_10[1001:11000,]
zeta_11 = poster_nmig_resclass12_mcmc2[[2]]$zeta_11[1001:11000,]
zeta_y_00 = poster_nmig_resclass12_mcmc2[[2]]$zeta_y_00[1001:11000,]
zeta_y_10 = poster_nmig_resclass12_mcmc2[[2]]$zeta_y_10[1001:11000,]
zeta_y_01 = poster_nmig_resclass12_mcmc2[[2]]$zeta_y_01[1001:11000,]
zeta_y_11 = poster_nmig_resclass12_mcmc2[[2]]$zeta_y_11[1001:11000,]
sigsq00_t = mean(c(nmig12_class_posterior[[2]]$sigsq00_pos[1001:11000]))
sigsq01_t = mean(c(nmig12_class_posterior[[2]]$sigsq00_pos[1001:11000]))
sigsq10_t = mean(c(nmig12_class_posterior[[2]]$sigsq10_pos[1001:11000]))
sigsq11_t = mean(c(nmig12_class_posterior[[2]]$sigsq11_pos[1001:11000]))

beta00_t = apply(beta00_pos,2,mean)
beta10_t = apply(beta10_pos,2,mean)
beta11_t = apply(beta11_pos,2,mean)
alpha11_t = apply(alpha11_pos,2,mean)
alpha10_t = apply(alpha10_pos,2,mean)
alpha01_t = apply(alpha01_pos,2,mean)
alpha00_t = apply(alpha00_pos,2,mean)
zeta_00_t =  apply(zeta_00,2,mean)
zeta_10_t =  apply(zeta_10,2,mean)
zeta_11_t =  apply(zeta_11,2,mean)
zeta_y_00_t = apply(zeta_y_00,2,mean)
zeta_y_10_t = apply(zeta_y_10,2,mean)
zeta_y_01_t = apply(zeta_y_01,2,mean)
zeta_y_11_t = apply(zeta_y_11,2,mean)



data_math = data_class_as_clu_list$data_math_comp
data_math_f=data_math[data_math$group_rcode==1|data_math$group_rcode==2,]
data_math_f$z = ifelse(data_math_f$group_rcode==1,0,1) ## Co_Ed F =0 vc F only=1
N_ori = dim(data_math_f)[1]
m = length(unique(data_math_f$clu))
n_ori = as.vector(table(data_math_f$clu))
probs = matrix(prop.table(table(data_math_f$clu,data_math_f$Ethnicity_rcode),margin=1),nrow=32,ncol=5)
probs_z = matrix(prop.table(table(data_math_f$clu,data_math_f$z),margin=1),nrow=32,ncol=2)

data_gene = function(dataid,m,n){
  set.seed(dataid)
  id = unlist(apply(as.matrix(n),1,seq,from=1,by=1))
  clu = rep(1:m,n)
  data = as.data.frame(cbind(clu,id))
  N = dim(data)[1]
  
  ## dummy variable/indicator matrix for cluster
  cluind = data$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  x1 = rnorm(N,mean(data_math_f$MathSS.2011_log),sd(data_math_f$MathSS.2011_log))
  data$x1 = scale(x1)
  x2=NULL
  for (cl in 1:m){
    gen =  sample(1:5,n[cl],replace=T,prob=probs[cl,])
    x2=c(x2,gen)
  }
  data$x2 = scale(x2)
  z = NULL
  for (cl in 1:m){
    gen_z = sample(0:1,n[cl],replace=T,prob = probs_z[cl,])
    z=c(z,gen_z)
  }
  data$z = z
  # principal strata
  ######## where x includes x1,x2,U
  prinstra_data = as.matrix(cbind(data$x1,data$x2))
  gstar_11 = prinstra_data%*%beta11_t+Istar%*%zeta_11_t+rnorm(N,0,1)
  gstar_10 = prinstra_data%*%beta10_t+Istar%*%zeta_10_t+rnorm(N,0,1)
  gstar_00 = prinstra_data%*%beta00_t+Istar%*%zeta_00_t+rnorm(N,0,1)
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
  data_y_0 = cbind(data$x1,data$x2,rep(0,N))
  data_y_1 = cbind(data$x1,data$x2,rep(1,N))
  for (i in 1:N){
    cl = data[i,]$clu
    if (data$G[i]==3){
      y0[i] = rnorm(1,mean = data_y_0[i,]%*%alpha00_t+zeta_y_00_t[cl],sd = sqrt(sigsq00_t))
      y1[i] = rnorm(1,mean = data_y_1[i,]%*%alpha00_t+zeta_y_00_t[cl],sd = sqrt(sigsq00_t))
    }else if (data$G[i]==2){
      y0[i] = rnorm(1,mean = data_y_0[i,]%*%alpha10_t+zeta_y_10_t[cl],sd = sqrt(sigsq10_t))
      y1[i] = rnorm(1,mean = data_y_1[i,]%*%alpha10_t+zeta_y_10_t[cl],sd = sqrt(sigsq10_t))
    }else if (data$G[i]==4){
      y0[i] = rnorm(1,mean = data_y_0[i,]%*%alpha01_t+zeta_y_01_t[cl],sd = sqrt(sigsq01_t))
      y1[i] = rnorm(1,mean = data_y_1[i,]%*%alpha01_t+zeta_y_01_t[cl],sd = sqrt(sigsq01_t))
    }else if (data$G[i]==1){
      y0[i] = rnorm(1,mean = data_y_0[i,]%*%alpha11_t+zeta_y_11_t[cl],sd = sqrt(sigsq11_t))
      y1[i] = rnorm(1,mean = data_y_1[i,]%*%alpha11_t+zeta_y_11_t[cl],sd = sqrt(sigsq11_t))
    }
  }
  y_obs = data$z*y1 + (1-data$z)*y0
  data = as.data.frame(cbind(data,y0,y1,y_obs))
  effect_simu = mean(data$y1-data$y0)
  ##### if change the data structure, double check the column numbers
  data_obs = data[,c(1:4,5,9,12)] 
  return (list(data=data,data_obs=data_obs,m=m,n=n,effect_simu=effect_simu))
}
simdata_m32_f= mclapply(1:100,data_gene,m=32,n=n_ori,mc.cores=20)

true_1100 = rep(-99,100)
for (i in 1:100){ 
  data = simdata_m32_f[[i]]$data
  data1100 = data[(data$G==1)|(data$G==3),]
  true_1100[i] = mean(data1100$y1-data1100$y0)
  }
print("mean(true_1100)")
print(mean(true_1100))
print("sd(true_1100)")
print(sd(true_1100))


save(simdata_m32_f,file="simdata_m32_f.RData")