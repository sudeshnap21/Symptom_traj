ip<-getProfiles(t="Time", y=c("EPDS", "PSS", "Sleephrs", "Fatigue", "jFatigue", "IL6_mean", "IL10_mean", "TNF_mean", "IL6_del", "IL10_del", "TNF_del" ), id="id", data=a1)

a1_time=a1[a1$Time == 1, ]
corr(a1_time$EPDS, a1_time$IL6)


mod<-GLMM_MCMC(y=a1[,c("EPDS","PSS", "Sleephrs", "Fatigue")], dist=c("gaussian","gaussian", "gaussian", "binomial(logit)" ), id=a1[,"id"], x=list(EPDS="empty", PSS="empty", sleephrs_avg="empty", Fatigue=a1[,"Weeks"]),z=list(EPDS=a1[,"Weeks"],PSS=a1[,"Weeks"], Sleephrs=a1[,"Weeks"], Fatigue="empty"), random.intercept=rep(TRUE,4), prior.b=list(Kmax=2), nMCMC=c(burn=500,keep=2000, thin=50, info=200), parallel=FALSE)


Devs1 <- Devs2 <- list()
PED <- data.frame()
for (K in 1:4) {
   set.seed(1234567)
   modK <- GLMM_MCMC(y=a1[,c("EPDS","PSS", "Sleephrs", "Fatigue")], 
                     dist=c("gaussian","gaussian", "gaussian", "binomial(logit)" ), 
                     id=a1[,"id"], x=list(EPDS=a1[,c("depression_hx", "TNF_mean", "TNF_del","IL10")], PSS="empty", 
                                          sleephrs_avg=a1[,c("IL6_mean", "IL6_del", "IL10")], Fatigue=a1[,"Weeks"]),
                     z=list(EPDS=a1[,"Weeks"],PSS=a1[,"Weeks"], Sleephrs=a1[,"Weeks"], Fatigue="empty"), 
                     random.intercept=rep(TRUE,4), prior.b=list(Kmax=K), nMCMC=c(burn=500,keep=2000, thin=50, info=200), parallel=FALSE)
   Devs1[[K]] <- modK$Deviance1
   Devs2[[K]] <- modK$Deviance2
   PED <- rbind(PED, modK$PED)
   colnames(PED) <- names(modK$PED)
}
print(PED, digits=6)
set.seed(1234567)
mod<-GLMM_MCMC(y=a1[,c("EPDS","PSS", "Sleephrs", "Fatigue")], dist=c("gaussian","gaussian", "gaussian", "binomial(logit)" ), id=a1[,"id"], x=list(EPDS=a1[,c("depression_hx", "TNF_mean", "TNF_del","IL10")], PSS="empty", sleephrs_avg=a1[,c("IL6_mean", "IL6_del", "IL10")], Fatigue=a1[,"Weeks"]),z=list(EPDS=a1[,"Weeks"],PSS=a1[,"Weeks"], Sleephrs=a1[,"Weeks"], Fatigue="empty"), random.intercept=rep(TRUE,4), prior.b=list(Kmax=2), nMCMC=c(burn=500,keep=2500, thin=50, info=200), parallel=FALSE)

mod1<-GLMM_MCMC(y=a1[,c("EPDS","PSS", "Sleephrs", "Fatigue")], dist=c("gaussian","gaussian", "gaussian", "binomial(logit)" ), id=a1[,"id"], x=list(EPDS=a1[,c("depression_hx", "TNF_mean", "TNF_del","IL10")], PSS="empty", sleephrs_avg=a1[,c("IL6_mean", "IL6_del", "IL10")], Fatigue=a1[,"Weeks"]),z=list(EPDS=a1[,"Weeks"],PSS=a1[,"Weeks"], Sleephrs=a1[,"Weeks"], Fatigue="empty"), random.intercept=rep(TRUE,4), prior.b=list(Kmax=1), nMCMC=c(burn=500,keep=2000, thin=50, info=200), parallel=FALSE)

