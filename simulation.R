#------------------------------------------------------------------------
# Filename: simulation.R
# Description:
# This R script contains functions to simulate true parameters and data.
#
# Contents of this file:
# 1. generatePrm(N, I, itemmean) 
#   - Simulate true model parameters for sample size N, test length I 
#   - itemmean can take on c(-2.197,4.394,4,3) or c(-2.197,4.394,4,5) 
#     for different item omission proportions.
# 2. simulateData(persons, pattern, items,prob_solv,Q)
#   - Simulate data based on the JCDM-RT-NA model
#
# Dependencies:
#  - mvtnorm: For multivariate normal distributions
#------------------------------------------------------------------------
library(mvtnorm)
generatePrm<-function(N, I, itemmean){
  set.seed(22724)
  K<-4
  sigma_person=matrix(c(1,-0.25,-0.7,
                        -0.25,0.16,0.2,
                        -0.7,0.2,3),nrow=3,ncol=3)
  persons<-rmvnorm(N,mean = rep(0,3),sigma_person)
  persons[,1]<-scale(persons[,1])
  persons[,2]<-persons[,2]-mean(persons[,2])
  persons[,3]<-persons[,3]-mean(persons[,3])
  theta<-persons[,1]
  tau<-persons[,2]
  psai<-persons[,3]
  
  lamda_K<-c(1,1,1,1)
  lamda_0<-c(1,0.5,-0.5,-1)
  exp.attribute = prob.attribute = matrix(NA,N,K)
  for(n in 1:N){
    for(k in 1:K){
      exp.attribute[n,k]=exp(lamda_K[k]*theta[n]+lamda_0[k])
      prob.attribute[n,k]=exp.attribute[n,k]/(1+exp.attribute[n,k])
    }}
  pattern=matrix(NA,N,K)
  for(n in 1:N){pattern[n,]=rbinom(K,1,prob.attribute[n,])}
  
  sigma_item=matrix(c(1,-0.8,-0.25,1,
                      -0.8,1,0.15,-0.25,
                      -0.25,0.15,0.25,-0.5,
                      1,-0.25,-0.5,5),nrow=4,ncol=4)
  items<-rmvnorm(I,mean = itemmean,sigma=sigma_item)
  
  return(list(persons = persons,
              pattern = pattern,
              items = items))
}

simulateData<-function(persons, pattern, items,prob_solv,Q){
  #--------------------------------------------
  #input: previously generated model parameters,
  #proportion of solution behavior, Q matrix
  #--------------------------------------------
  N<-nrow(persons)
  I<-nrow(items)
  K<-ncol(Q)
  theta<-persons[,1]
  tau<-persons[,2]
  psai<-persons[,3]
  
  beta<-items[,1]
  delta<-items[,2]
  zeta<-items[,3]
  gamma<-items[,4]
  
  D<-matrix(NA, nrow = N, ncol = I)
  
  for (n in 1:N){
    for (i in 1:I){
      D[n,i]<-rbinom(1,1,prob_solv)
    }
  }
  
  #omission probability for engaged people
  
  prob.omit1<-matrix(NA,nrow = N, ncol = I)
  linear1<-matrix(NA,nrow = N, ncol = I)
  for (n in 1:N){
    for (i in 1:I){
      linear1[n,i] <- psai[n]-gamma[i]
      prob.omit1[n,i]<-exp(linear1[n,i])/(1+exp(linear1[n,i]))
    }
  }
  
  # omission binary indicators
  d<-matrix(NA,nrow = N, ncol = I)
  for (n in 1:N){
    for(i in 1:I){
      if(D[n,i]==1){d[n,i]=rbinom(1,1,prob.omit1[n,i])}
      if(D[n,i]==0){d[n,i]= rbinom(1,1,0.1)}
    }
  }

  #generate response data
  exp=prob=matrix(NA,N,I)
  w=array(NA,c(N,I,K))
  for(n in 1:N){
    for(m in 1:I){
      for(k in 1:K){w[n,m,k]= pattern[n,k]^Q[m,k]}#item require and person master
      exp[n,m]=exp(beta[m]+delta[m]*prod(w[n,m,]))
      prob[n,m]=exp[n,m]/(1+exp[n,m])
    }
  }
 
  Score=matrix(NA,N,I)
  for (n in 1:N){
    for(m in 1:I){
      if(D[n,m]==1&d[n,m]==0){Score[n,m]=rbinom(1,1,prob[n,m])}
      if(D[n,m]==1&d[n,m]==1){Score[n,m]=NA}
      if(D[n,m]==0&d[n,m]==0){Score[n,m]=rbinom(1,1,0.2)}
      if(D[n,m]==0&d[n,m]==1){Score[n,m]=NA}
    }
  }
  
  #generate RT data
  muD<-2
  sigmaD<-1.4^2
  sigmaE<-0.25
  logT<-matrix(NA,nrow = N, ncol = I)  

  for (n in 1:N){
    for (i in 1:I){
      if (D[n,i]==1)
      {logT[n,i]<-rnorm(1, zeta[i]-tau[n],sqrt(sigmaE))}
      if (D[n,i]==0)
      {logT[n,i]<-rnorm(1, muD,sqrt(sigmaD))}
    }
  }
 return(list(d=d,Score=Score, logT=logT))
  
}
