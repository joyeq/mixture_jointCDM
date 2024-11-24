load("empiricalData.RData")

jags.inits<- NULL
jags.data1<- list(N=N, I=I, d=d, logT=Timelog, Score=Y,K=K,totna=totna,nr=nr,Q=Q)
jags.parameters1 <- c("tau","zeta", "alpha","sigmaE","sigmaD","beta","delta", "pvalue_0_item","pvalue_rt_item","pvalue_d_item","sigma_item","item_mu",
                      "g", "muC","gamma0", "pvalue_0","pvalue_rt","pvalue_d","theta","psai","sigma_person","lambda_K","lambda_0","romit","D","Rho","prob_solv")

sim1<- jags(data=jags.data1, inits=jags.inits, parameters.to.save=jags.parameters1,
            model.file = "JCDMRTNA.txt", n.chains = 2,n.iter = 60000,n.burnin = 30000, DIC = TRUE)

result_new<-sim1$BUGSoutput$summary


#Plot posterior probabilities of solution behavior
alpha_new_data<-as.data.frame(matrix(result_new[grep("alpha",result_new$X),6],nrow=N,ncol=K))%>%mutate(., pattern = paste(V1,V2,V3,V4)) 
D_new<-as.data.frame(matrix(result_new[grep("D",result_new$X),6][-(N*I+1)],nrow=N,ncol=I))
P_new<-as.data.frame(matrix(result_new[grep("D",result_new$X),2][-(N*I+1)],nrow=N,ncol=I))

Prob_data1<-(P_new[d==1])
Prob_data2<-(P_new[d==0])
label1<-rep(1, length(Prob_data1))
label2<-rep(0, length(Prob_data2))

appender <- function(string) TeX(paste("$P(\\Delta_{ij} = 1|D_{ij} = $", string,")")) 
Prob_data<-cbind.data.frame(c(Prob_data1,Prob_data2),c(label1,label2))
colnames(Prob_data)<-c("value","label")

ggplot( data=Prob_data,aes(x=value)) +
  geom_histogram() +
  theme_bw()+
  xlab(NULL) +
  ylab("Frequency") +
  facet_wrap(~label,scales="free_y",labeller = as_labeller(appender, default = label_parsed))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        legend.position="bottom",
        legend.direction = "horizontal", 
        strip.text = element_text(size = 20),
        axis.title=element_text(size=20),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5,size = 16),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5,size = 16))
