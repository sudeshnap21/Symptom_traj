#Load Package: mixAK and dependencies***/
library(colorspace)
library(lme4)
library(Matrix)
library(mixAK)

setwd(file.path("T:\\Paper\\LCMM\\paper1"))
data<- read.csv(file="data_final2.csv",header=TRUE,sep=",")
a1=data.frame(data)
colnames(a1)[1]="id"

a1$jFatigue=jitter(a1$Fatigue, amount=0.05)

ip<-getProfiles(t="Time", y=c("EPDS", "PSS", "Sleephrs", "Fatigue", "jFatigue", "IL6_mean", "IL10_mean", "TNF_mean", "IL6_del", "IL10_del", "TNF_del" ), id="id", data=a1)

plotProfiles(ip=ip, data=a1,var="jFatigue", tvar="Time", highlight=c(10,60),ylab="Fatigue(yes/no)", xlab="Time", lines=FALSE, points=TRUE)
plotProfiles(ip=ip, data=a1,var="EPDS", tvar="Time", highlight=c(10, 60),ylab="Total EPDS score", xlab="Time", lwd=2, lwd.highlight=5)

#trial 1: Symptoms model
set.seed(1234567)
mod<-GLMM_MCMC(y=a1[,c("EPDS","PSS", "Sleephrs", "Fatigue")], dist=c("gaussian","gaussian", "gaussian", "binomial(logit)" ), id=a1[,"id"], x=list(EPDS="empty", PSS="empty", sleephrs_avg="empty", Fatigue=a1[,"Weeks"]),z=list(EPDS=a1[,"Weeks"],PSS=a1[,"Weeks"], Sleephrs=a1[,"Weeks"], Fatigue="empty"), random.intercept=rep(TRUE,4), prior.b=list(Kmax=2), nMCMC=c(burn=200,keep=2000, thin=50, info=200), parallel=FALSE)

#Trial 2: Symptoms model+clinical covariates
set.seed(1234567)
mod<-GLMM_MCMC(y=a1[,c("EPDS","PSS", "Sleephrs", "Fatigue")], dist=c("gaussian","gaussian", "gaussian", "binomial(logit)" ), id=a1[,"id"], x=list(EPDS=a1[,c("depression_hx","depression_famhx")], PSS="empty", sleephrs_avg="empty", Fatigue=a1[,"Weeks"]),z=list(EPDS=a1[,"Weeks"],PSS=a1[,"Weeks"], Sleephrs=a1[,"Weeks"], Fatigue="empty"), random.intercept=rep(TRUE,4), prior.b=list(Kmax=2), nMCMC=c(burn=200,keep=2000, thin=50, info=200), parallel=FALSE)

#Trial 3: symptoms model+ clinical covariates +biomarkers
set.seed(1234567)
mod<-GLMM_MCMC(y=a1[,c("EPDS","PSS", "Sleephrs", "Fatigue")], dist=c("gaussian","gaussian", "gaussian", "binomial(logit)" ), id=a1[,"id"], x=list(EPDS=a1[,c("depression_hx", "TNF_mean", "TNF_del","IL10")], PSS="empty", sleephrs_avg=a1[,c("IL6_mean", "IL6_del")], Fatigue=a1[,"Weeks"]),z=list(EPDS=a1[,"Weeks"],PSS=a1[,"Weeks"], Sleephrs=a1[,"Weeks"], Fatigue="empty"), random.intercept=rep(TRUE,4), prior.b=list(Kmax=2), nMCMC=c(burn=500,keep=2000, thin=50, info=200), parallel=FALSE)


#analysing the output

library("coda")
mod <- NMixRelabel(mod, type = "stephens", keep.comp.prob = TRUE)
print(mod[[1]]$Deviance[1:10], digits = 9)

muSamp1<-NMixChainComp(mod[[1]], relabel=TRUE, param="mu_b")
print(muSamp1[1:3,])
tracePlots(mod, param = "Deviance")
tracePlots(mod, param = "mu_b", relabel = TRUE)
print(mod)

#saving the r object: read it back as load(file="model.Rdata")#

save(mod, file="model.Rdata")

#Mixture analysis
NMixSummComp(mod[[1]])
summary(mcmc(muSamp1))


#clustering specific code
delta<-2
tpred=seq(0,24,by=delta)
Npred=length(tpred)
mepds=cbind(rep(0,Npred),rep(0,Npred),rep(30,Npred))
tage=rep(30,Npred)
fit <- fitted(mod[[1]], x = list(mepds,tage, tage, tpred), z = list(tpred, tpred, tpred, "empty"), glmer = TRUE)

fit <- fitted(mod[[1]], x = list("empty","empty", "empty", tpred), z = list(tpred, tpred, tpred, "empty"), glmer = TRUE)
names(fit) <- c("EPDS", "PSS", "Sleephrs",  "Fatigue")
print(fit[["EPDS"]][1:3, ], digits = 5)

K <- mod[[1]]$prior.b$Kmax
clCOL <- c("blue", "red")
plotProfiles(ip = ip, data = a1, var = "EPDS", tvar = "Time", points=TRUE ,col = "azure3", xlab = "Time", ylab = "Total EPDS", lwd=2)
for (k in 1:K) lines(tpred, fit[["EPDS"]][, k], col = clCOL[k], lwd = 5)


K <- mod[[1]]$prior.b$Kmax
clCOL <- c("blue", "red")
plotProfiles(ip = ip, data = a1, var = "PSS", tvar = "Time", points=TRUE ,col = "azure3", xlab = "Time", ylab = "PSS", lwd=2)
for (k in 1:K) lines(tpred, fit[["PSS"]][, k], col = clCOL[k], lwd = 5)

K <- mod[[1]]$prior.b$Kmax
clCOL <- c("blue", "red")
plotProfiles(ip = ip, data = a1, var = "Sleephrs", tvar = "Time", points=TRUE ,col = "azure3", xlab = "Time", ylab = "Sleep hours", lwd=2)
for (k in 1:K) lines(tpred, fit[["Sleephrs"]][, k], col = clCOL[k], lwd = 5)


K <- mod[[1]]$prior.b$Kmax
clCOL <- c("blue", "red")
plotProfiles(ip = ip, data = a1, var = "jFatigue", tvar = "Time", points=TRUE ,col = "azure3", xlab = "Time", ylab = "Fatigue (yes/no)", lines=FALSE)
for (k in 1:K) lines(tpred, fit[["Fatigue"]][, k], col = clCOL[k], lwd = 5)

#median post prob

print(mod[[1]]$quant.comp.prob[["50%"]][1:3, ])
groupMed <- apply(mod[[1]]$quant.comp.prob[["50%"]], 1, which.max)
pMed <- apply(mod[[1]]$quant.comp.prob[["50%"]], 1, max)
table(groupMed)

TAB <- table(a1$id)
a1$groupMed <- factor(rep(groupMed, TAB))
GCOL <-c("cornflowerblue", "coral") 
names(GCOL) <- levels(a1$groupMed)

par(cex.axis=1.5, cex.lab=1.5, cex.main=1.2, cex.sub=1)
ip<-getProfiles(t="Time", y=c("EPDS", "PSS", "Sleephrs", "jFatigue", "groupMed"), id="id", data=a1)

plotProfiles(ip = ip, data = a1, var = "EPDS", tvar = "Time", gvar = "groupMed", col = GCOL, auto.layout = FALSE, xlab = "Time (months)",ylab = "EPDS", lwd=2)
plotProfiles(ip = ip, data = a1, var = "PSS", tvar = "Time", gvar = "groupMed", col = GCOL, auto.layout = FALSE, xlab = "Time (months)",ylab = "PSS", lwd=2)
plotProfiles(ip = ip, data = a1, var = "Fatigue", tvar = "Time", gvar = "groupMed", col = GCOL, auto.layout = FALSE, xlab = "Time",ylab = "Fatigue")
plotProfiles(ip = ip, data = a1, var = "Sleephrs", tvar = "Time", gvar = "groupMed", col = GCOL, auto.layout = FALSE, xlab = "Time",ylab = "Average sleep hours", lwd=2)



a1_time=a1[a1$Time == 1, ]
by(a1_time$sleephrs_avg,a1_time$groupMed, sd, na.rm=T)
t.test(a1_time$sleephrs_avg~a1_time$groupMed)

by(a1_time$age,a1_time$groupMed, sd, na.rm=T)
t.test(a1_time$age~a1_time$groupMed)

tbl = table(a1_time$depression_hx, a1_time$groupMed)
chisq.test(tbl)


####mean Table by clusters#####
  
  a1_sub=a1[a1$Time == 7, ]
by(a1_sub$EPDS,a1_sub$groupMed, mean, na.rm=T)
by(a1_sub$EPDS,a1_sub$groupMed, sd, na.rm=T)

by(a1_sub$PSS,a1_sub$groupMed, mean, na.rm=T)
by(a1_sub$PSS,a1_sub$groupMed, sd, na.rm=T)

by(a1_sub$Sleephrs,a1_sub$groupMed, mean, na.rm=T)
by(a1_sub$Sleephrs,a1_sub$groupMed, sd, na.rm=T)

by(a1_sub$Fatigue,a1_sub$groupMed, mean, na.rm=T)
by(a1_sub$Fatigue,a1_sub$groupMed, sd, na.rm=T)
