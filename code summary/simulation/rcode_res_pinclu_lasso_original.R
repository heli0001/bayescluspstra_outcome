
load("simdata_m50.RData")
#load("poster_lassom50_st3_probinclu.RData") # first MCMC result
load("poster_lassom50.RData")


library(survival)
library(coda) # for HPD interval 
library(data.table)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

iter=11000
burn=1000


get_pinclu_res = function(dataid,simdata,iter,burn){
####### recall result from first MCMC 
  pstra11_all = poster_lassom50[[dataid]]$zeta_11[(burn+1):iter,]
  pstra11 = HPDinterval(as.mcmc(pstra11_all), prob = 0.95)
  
  pstra10_all = poster_lassom50[[dataid]]$zeta_10[(burn+1):iter,]
  pstra10 = HPDinterval(as.mcmc(pstra10_all), prob = 0.95)

  pstra00_all = poster_lassom50[[dataid]]$zeta_00[(burn+1):iter,]
  pstra00 = HPDinterval(as.mcmc(pstra00_all), prob = 0.95)

  pstra_y11_all = poster_lassom50[[dataid]]$zeta_y_11[(burn+1):iter,] 
  pstra_y11 = HPDinterval(as.mcmc(pstra_y11_all), prob = 0.95)
  
  pstra_y10_all = poster_lassom50[[dataid]]$zeta_y_10[(burn+1):iter,] 
  pstra_y10 = HPDinterval(as.mcmc(pstra_y10_all), prob = 0.95)
  
  pstra_y00_all = poster_lassom50[[dataid]]$zeta_y_00[(burn+1):iter,]
  pstra_y00 = HPDinterval(as.mcmc(pstra_y00_all), prob = 0.95)
  
  pstra_y01_all = poster_lassom50[[dataid]]$zeta_y_01[(burn+1):iter,]
  pstra_y01 = HPDinterval(as.mcmc(pstra_y01_all), prob = 0.95)
  
  
  ### if contain 0, then not include  
  p_inc11 = ifelse(between(0,pstra11[,1],pstra11[,2]),0,1)
  p_inc10 = ifelse(between(0,pstra10[,1],pstra11[,2]),0,1)
  p_inc00 = ifelse(between(0,pstra00[,1],pstra11[,2]),0,1)
  p_inc_y11 = ifelse(between(0,pstra_y11[,1],pstra11[,2]),0,1)
  p_inc_y10 = ifelse(between(0,pstra_y10[,1],pstra11[,2]),0,1)
  p_inc_y00 = ifelse(between(0,pstra_y00[,1],pstra11[,2]),0,1)
  p_inc_y01 = ifelse(between(0,pstra_y01[,1],pstra11[,2]),0,1)

  
  print(dataid)
  
  return(list(p_inc11 = p_inc11, p_inc10 = p_inc10, p_inc00 = p_inc00, p_inc_y11 = p_inc_y11, p_inc_y10 = p_inc_y10, p_inc_y00 = p_inc_y00, p_inc_y01 = p_inc_y01))
}
get_pinclu_resm50_lasso = mclapply(1:100,get_pinclu_res,simdata_m50,iter=iter,burn=burn,mc.cores=20)
save(get_pinclu_resm50_lasso,file="get_pinclu_resm50_lasso.RData")
