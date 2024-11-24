#------------------------------------------------------------------------
# Filename: estimation.R
# Description:
# This R script uses R2jags to estimate the data using the JCDM-RT-NA model,
# and two competing models: RT-DINA-RG and JRT-DINA 
# Dependencies:
#  - R2jags: For fully Bayesian estimation
#------------------------------------------------------------------------

library(R2jags)
#JCDM-RT-DINA
jags.inits<- NULL
jags.data1<- list(N=N, I=I, d=d, logT=logT, Score=Score,K=K,totna=d,Q=Q)
jags.parameters1 <- c("tau","zeta", "alpha","sigmaE","sigmaD","beta","delta","sigma_item","item_mu",
                      "g", "muC","gamma0", "pvalue_0","pvalue_rt","pvalue_d","theta","psai","sigma_person","lambda_K","lambda_0","romit","D","Rho","prob_solv")
t1 <- Sys.time()
sim1<- jags(data=jags.data1, inits=jags.inits, parameters.to.save=jags.parameters1,
            model.file = "JCDMRGNA.txt", n.chains = 2,n.iter = 60000,n.thin = 1,n.burnin = 30000, DIC = TRUE)
end.time1 <- Sys.time()-t1

#RT-DINA-RG
jags.inits<- NULL
jags.data2<- list(N=N, I=I,logT=logT, Score=Score,K=K,totna=d,Q=Q)
jags.parameters2 <- c("tau","zeta", "alpha","sigmaE","sigmaD", "beta","delta","muC","D","g","prob_solv","Rho",
                      "pvalue_0","pvalue_rt","theta","sigma_person","lambda_K","lambda_0","sigma_item")
t2 <- Sys.time()
sim2<- jags(data=jags.data2, inits=jags.inits, parameters.to.save=jags.parameters2,
            model.file = "JRTRG.txt", n.chains = 2,n.iter = 60000,n.thin = 1,n.burnin = 30000, DIC = TRUE)
end.time2 <- Sys.time()-t2

#JRT-DINA
jags.inits<- NULL
jags.data1<- list(N=N, I=I, logT=logT, Score=Score,K=K,Q=Q)
jags.parameters1 <- c("tau","zeta", "alpha","Sigma_omiga","beta","delta",
                      "pvalue_0","pvalue_rt","theta","sigma_theta","lambda_K","lambda_0","sigma_item")
t3 <- Sys.time()
sim3<- jags(data=jags.data3, inits=jags.inits, parameters.to.save=jags.parameters3,
            model.file = "JRT_DINA.txt", n.chains = 2,n.iter = 60000,n.thin = 1,n.burnin = 30000, DIC = TRUE)
end.time3 <- Sys.time()-t3
