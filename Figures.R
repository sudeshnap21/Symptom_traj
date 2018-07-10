#Paper : Figure 1

setwd(file.path("T:\\Paper\\LCMM\\paper1"))
data<- read.csv(file="data_final2.csv",header=TRUE,sep=",")
a1=data.frame(data)
colnames(a1)[1]="id"

a1$jFatigue=jitter(a1$Fatigue, amount=0.075)
ip<-getProfiles(t="Weeks", y=c("EPDS", "PSS", "Sleephrs", "jFatigue", "IL6", "IL10", "TNF", "groupMed"), id="id", data=a1)

delta<-0.1
tpred=seq(0,24,by=delta)
Npred=length(tpred)

mu_TNF_mean=rep(mean(a1$TNF_mean, na.rm=T),Npred)
mu_TNF_del=rep(0,Npred)
mu_IL6_mean=rep(mean(a1$IL6_mean, na.rm=T),Npred)
mu_IL6_del=rep(0,Npred)
mu_IL10=rep(mean(a1$IL10, na.rm=T), Npred)
mepds=cbind(rep(0,Npred),mu_TNF_mean, mu_TNF_del,mu_IL10)
msleephrs=cbind(mu_IL6_mean, mu_IL6_del, mu_IL10)

fit <- fitted(mod[[1]], x = list(mepds,"empty", msleephrs, tpred), z = list(tpred, tpred, tpred, "empty"), glmer = TRUE)
names(fit) <- c("EPDS", "PSS", "Sleephrs", "Fatigue")

ip<-getProfiles(t="Weeks", y=c("EPDS", "PSS", "Sleephrs", "IL6","IL10", "TNF", "Fatigue", "jFatigue"), id="id", data=a1)

GCOL=c("pink", "skyblue")
clCOL=c("red", "blue")

plotProfiles(ip = ip, data = a1, var = "PSS", tvar = "Weeks", gvar = "groupMed", col = GCOL, auto.layout = FALSE, xlab = "Weeks",ylab = "PSS", lwd=2)
for (k in 1:K) lines(tpred, fit[["PSS"]][, k], col = clCOL[k], lwd = 3, lty=1)

plotProfiles(ip = ip, data = a1, var = "jFatigue", tvar = "Weeks", gvar = "groupMed", col = GCOL, auto.layout = FALSE, xlab = "Weeks",ylab = "Fatigue", lwd=2,lty=1, points=TRUE, lty=2)
for (k in 1:K) lines(tpred, fit[["Fatigue"]][, k], col = clCOL[k], lwd = 3, lty=1)

lotProfiles(ip = ip, data = a1, var = "Sleephrs", tvar = "Weeks", gvar = "groupMed", col = GCOL, auto.layout = FALSE, xlab = "Weeks",ylab = "Sleep hours", lwd=2)
for (k in 1:K) lines(tpred, fit[["Sleephrs"]][, k], col = clCOL[k], lwd = 3, lty=1)

legend(15, 12, legend=c("cluster 1", "cluster 2"),
       col=clCOL, lty=1, cex=1)
#Add error bars: arrows(x-sdev, y, x+sdev, y, length=0.05, angle=90, code=3)#

#95%CI of fitted trajectory
LCL95=fitted(mod[[1]], x = list(mepds,"empty", msleephrs, tpred), z = list(tpred, tpred, tpred, "empty"), statistic="2.5%")
UCL95=fitted(mod[[1]], x = list(mepds,"empty", msleephrs, tpred), z = list(tpred, tpred, tpred, "empty"), statistic="97.5%")

GCOL=c("pink", "skyblue")
clCOL=c("red", "blue")
plotProfiles(ip = ip, data = a1, var = "EPDS", tvar = "Weeks", gvar = "groupMed", col = GCOL, auto.layout = FALSE, xlab = "Weeks",ylab = "EPDS", lwd=2)
for (k in 1:K) lines(tpred, fit[["EPDS"]][, k], col = clCOL[k], lwd = 3, lty=1)
for (k in 1:K) lines(tpred, LCL95[["EPDS"]][, k], col = clCOL[k], lwd = 3, lty=1)
for (k in 1:K) lines(tpred, UCL95[["EPDS"]][, k], col = clCOL[k], lwd = 3, lty=1)
