########################## birth wt data mimic ######################
# this file is for simulation mimic based on birthwt data application. The total size is N = 1068, m=99, cluster = age & bmi category, each has 1 to 32, 
# dataset: from "more_cluster_induction_with_gesage.csv"
# X: three other covariate, marital status, Mother's Education, race/ethnicity, since all are categorial value, based on the probalily in each cluster to randomly generate over 100 replicate, another one is gestational weeks age
# beta, alpha: use the data application nmig parameter from ("poster_nmig_res") result to generate the data
# Z: treatment: marital health status (in terms of smoking behaviour, hypertension, hypertension elapcemia), healthy =0 vc poor =1, based on the probalily in each cluster to randomly generate over 100 replicate 
# m and ni in each cluster is the same as data application for 100 replicates
# S: post variable is whether or nor induction of labor 


library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

load("poster_nmig_agebmi_with_gesage_res_v4.RData")
load("poster_nmig_resagebmi_mcmc2_v4.RData")

#################### for principal strata 
## 1. fixed effect
beta00_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$beta00_pos[1001:11000,],2,mean)#x1,x2,x3,x4
beta10_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$beta10_pos[1001:11000,],2,mean)
beta11_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$beta11_pos[1001:11000,],2,mean)
## cluster effect
zeta_11_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$zeta_11[1001:11000,],2,mean)
zeta_10_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$zeta_10[1001:11000,],2,mean)
zeta_00_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$zeta_00[1001:11000,],2,mean)


################### for outcomes 
#
## 1. fixed effect
alpha00_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$alpha00_pos[1001:11000,],2,mean) #x1,x2,x3,x4,z
alpha10_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$alpha10_pos[1001:11000,],2,mean)
alpha01_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$alpha01_pos[1001:11000,],2,mean)
alpha11_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$alpha11_pos[1001:11000,],2,mean)
sigsq00_t = mean(poster_nmig_agebmi_with_gesage_res_v4[[1]]$sigsq00_pos[1001:11000])
sigsq01_t = mean(poster_nmig_agebmi_with_gesage_res_v4[[1]]$sigsq01_pos[1001:11000])
sigsq10_t = mean(poster_nmig_agebmi_with_gesage_res_v4[[1]]$sigsq10_pos[1001:11000])
sigsq11_t = mean(poster_nmig_agebmi_with_gesage_res_v4[[1]]$sigsq11_pos[1001:11000])
## 2. cluster effect
zeta_y_11_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$zeta_y_11[1001:11000,],2,mean)
zeta_y_10_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$zeta_y_10[1001:11000,],2,mean)
zeta_y_01_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$zeta_y_01[1001:11000,],2,mean)
zeta_y_00_t = apply(poster_nmig_resagebmi_mcmc2_v4[[1]]$zeta_y_00[1001:11000,],2,mean)

sub_comp = read.csv("more_cluster_induction_with_gesage.csv",header=T,fill=T,sep=",")
N_ori = dim(sub_comp)[1]
m = length(unique(sub_comp$cluster))
n_ori = as.vector(table(sub_comp$cluster))
# Marital Status, married vs not-married
probs_x1 = matrix(prop.table(table(sub_comp$clu,sub_comp$dmar),margin=1),nrow=99,ncol=2) 
# Mother's Education 1 8th grade or less, 9 levels
probs_x2 = matrix(prop.table(table(sub_comp$clu,sub_comp$meduc),margin=1),nrow=99,ncol=9)
# Mother's Race Recode 6, values = 10,20,30,40,41,51,61
probs_x3 = matrix(prop.table(table(sub_comp$clu,sub_comp$mrace6),margin=1),nrow=99,ncol=7)
# Gestational weeks


probs_z = matrix(prop.table(table(sub_comp$clu,sub_comp$z),margin=1),nrow=99,ncol=2)

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
  
  x1 = NULL
  for (cl in 1:m){
    gen_x1 =  sample(1:2,n[cl],replace=T,prob=probs_x1[cl,])
    x1=c(x1,gen_x1)
  }
  data$x1 = x1
  x2=NULL
  for (cl in 1:m){
    gen_x2 =  sample(1:9,n[cl],replace=T,prob=probs_x2[cl,])
    x2=c(x2,gen_x2)
  }
  data$x2 = x2
  x3=NULL
  for (cl in 1:m){
    gen_x3 =  sample(c(10,20,30,40,41,51,61),n[cl],replace=T,prob=probs_x3[cl,])
    x3=c(x3,gen_x3)
  }
  data$x3 = x3
  
  x4 = rnorm(N,mean(sub_comp$combgest),sd(sub_comp$combgest)) # conti  use mean and sd
  data$x4 = x4
  #data$x4 = sub_comp$combgest
  data$x5 = sub_comp$mager
  data$x6 = sub_comp$bmi_r
  z = NULL
  for (cl in 1:m){
    gen_z = sample(0:1,n[cl],replace=T,prob = probs_z[cl,])
    z=c(z,gen_z)
  }
  data$z = z
  # principal strata
  ######## where x includes x1,x2,x3,x4,U
  prinstra_data = as.matrix(cbind(data$x1,data$x2,data$x3,data$x4))
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
  data_y_0 = cbind(data$x1,data$x2,data$x3,data$x4,rep(0,N))
  data_y_1 = cbind(data$x1,data$x2,data$x3,data$x4,rep(1,N))
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
  data_obs = data[,c(1:6,9,13,16)] 
  data_obs_naive_linear = data[,c(1:8,9,13,16)] 
  return (list(data=data,data_obs=data_obs,data_obs_naive_linear=data_obs_naive_linear,m=m,n=n,effect_simu=effect_simu))
}
simdata_bw_m99= mclapply(1:100,data_gene,m=m,n=n_ori,mc.cores=20)
true_00 = rep(-99,100)
true_1100 = rep(-99,100)
for (i in 1:100){ 
  data = simdata_bw_m99[[i]]$data
  data1100 = data[(data$G==1)|(data$G==3),]
  true_1100[i] = mean(data1100$y1-data1100$y0)
  
  data = simdata_bw_m99[[i]]$data
  data00 = data[data$G==3,]
  true_00[i] = mean(data00$y1-data00$y0)
}
print("mean(true_1100)")
print(mean(true_1100))
print("sd(true_1100)")
print(sd(true_1100))

print("mean(true_00)")
print(mean(true_00))
print("sd(true_00)")
print(sd(true_00))

save(simdata_bw_m99,file="simdata_bw_m99.RData")
