

library(blma)

source("./Common/BayesFactors.R")
source("./Common/BayesFactorsGonN.R",chdir=TRUE)
source("./Common/BayesFactorsZS.R")


plot_all_BayesFactors <- function(n,R2,vp) {

	vBF1 <- c(); for (j in 1:length(vp)) { vBF1[j] <- liang_g2(n,p_gamma=vp[j],R2); }
	vBF2 <- c(); for (j in 1:length(vp)) { vBF2[j] <- liang_g_n_quad(n,p_gamma=vp[j],R2); }
	vBF3 <- c(); for (j in 1:length(vp)) { vBF3[j] <- robust_bayarri2(n,p_gamma=vp[j],R2); }
	vBF4 <- c(); for (j in 1:length(vp)) { vBF4[j] <- ZE(n,p_gamma=vp[j],R2); }
	vBF5 <- c(); for (j in 1:length(vp)) { vBF5[j] <- BIC(n,p_gamma=vp[j],R2); }
	vBF6 <- c(); for (j in 1:length(vp)) { vBF6[j] <- log_BF_Zellner_Siow_quad(n,p_gamma=vp[j],R2); }

	vBF1[!is.finite(vBF1)] <- NA
	vBF2[!is.finite(vBF2)] <- NA
	vBF3[!is.finite(vBF3)] <- NA
	vBF4[!is.finite(vBF4)] <- NA
	vBF5[!is.finite(vBF5)] <- NA
	vBF6[!is.finite(vBF6)] <- NA

	xlim <- range(vp)
	ylim <- range(c(vBF1,vBF2,vBF3,vBF4,vBF5),na.rm=TRUE)

	plot(NA,type="n",xlim=xlim,ylim=ylim,xlab="p",ylab="log BF",
		main=paste("log BF - n = ",n,", R-sq=",R2),cex.lab=1.5,cex.main=1.5)

	lines(vp,vBF1,lwd=2,col="black")
	lines(vp,vBF2,lwd=2,col="blue")
	lines(vp,vBF3,lwd=2,col="red")
	lines(vp,vBF4,lwd=2,col="green")
	lines(vp,vBF5,lwd=2,col="orange")
	lines(vp,vBF6,lwd=2,col="purple")

	lines(vp,vp*0,col="gray",lwd=2)
}

pdf("BayesFactors.pdf",width=8,height=8)

par(mfrow=c(3,3))

plot_all_BayesFactors(n=50,R2=0.1,vp=1:20)
plot_all_BayesFactors(n=200,R2=0.1,vp=1:20)
plot_all_BayesFactors(n=1000,R2=0.1,vp=1:20)

plot_all_BayesFactors(n=50,R2=0.5,vp=1:20)
plot_all_BayesFactors(n=200,R2=0.5,vp=1:20)
plot_all_BayesFactors(n=1000,R2=0.5,vp=1:20)

plot_all_BayesFactors(n=50,R2=0.9,vp=1:20)
plot_all_BayesFactors(n=200,R2=0.9,vp=1:20)
plot_all_BayesFactors(n=1000,R2=0.9,vp=1:20)

dev.off()

