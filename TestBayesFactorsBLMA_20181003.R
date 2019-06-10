

library(blma)
library(blmatest)
library(pracma)
library(BAS)

####################################################################################################

mutual_diff <- function(x) {
	N <- length(x)
	D <- matrix(NA,N,N)
	for (i in 1:N) {
		for (j in 1:N) {
			D[i,j] = x[i] - x[j]
		}
	}
	return(D)
}

####################################################################################################

check_hyper_g <- function(n,p_gamma,R2)
{
	vals_hyper_g <- c(
		liang_g1(n,p_gamma,R2),  # hyper-g via eq. 10 in package blma
		liang_g2(n,p_gamma,R2),  # hyper-g via eq. 11 in package blma
		log.BF.g.naive(R2, n, p_gamma), # hyper-g via eq. 10 in BayesFactors.R
		log.BF.g.safe1(R2, n, p_gamma), # hyper-g via eq. 11 in BayesFactors.R with log.hyperg_2F1
		log.BF.g.safe2(R2, n, p_gamma) # hyper-g via eq. 11 in BayesFactors.R
	)

	D <- mutual_diff(vals_hyper_g)
	return(D)
}

####################################################################################################

# hyper-g/n

check_hyper_g_n <- function(n,p_gamma,R2)
{
	a <- 3
	quadrule = gaussLegendre(n=1000,a=0,b=1)

	vals_hyper_g_n <- c(
		liang_g_n_appell(n,p_gamma,R2),                 # hyper-g/n via eq. 14 in package blma using gsl for appell
		liang_g_n_quad(n,p_gamma,R2),                   # hyper-g/n via eq. 13 in package blma with Gauss-Legendre
		liang_g_n_approx(n,p_gamma,R2),                 # hyper-g/n via old appell fn approx.
		log.BF.g.on.n.naive(R2, n, p_gamma, flag=1)$vals,      # hyper-g/n via eq. 14 in BayesFactors.R using gsl for appel
		log.BF.g.on.n.quad(R2, n, p_gamma, a, quadrule),       # hyper-g/n via eq. 13 in BayesFactors.R with Gauss-Legendre
		log.BF.g.on.n.Laplace(R2, n, p_gamma)                  # hyper-g/n via eq. 13 in BayesFactors.R with Laplace approx.
	)

	D <- mutual_diff(vals_hyper_g_n)
	return(D)
}

####################################################################################################

# Robust
check_robust <- function(n,p_gamma,R2)
{
	vals_robust <- c(
		robust_bayarri1(n,p_gamma,R2), # 1 Robust via eq. 17 in package blma
		robust_bayarri2(n,p_gamma,R2), # 2 Robust via eq. 18 in package blma
		log.BF.robust.naive1(R2, n, p_gamma), # 3 Robust via eq. 17 in BayesFactors.R with gsl implementation of 2F1
		log.BF.robust.naive2(R2, n, p_gamma), # 4 Robust via eq. 17 in BayesFactors.R with BAS implementation of 2F1
		log.BF.robust.naive3(R2, n, p_gamma), # 5 Robust via eq. 17 in BayesFactors.R with BAS implementation of 2F1 using Laplace method
		log.BF.robust.quad1(R2, n, p_gamma),  # 6 Robust via eq. 17 in BayesFactors.R with trapezoidal integration
		log.BF.robust.quad2(R2, n, p_gamma),  # 7 Robust via eq. 17 in BayesFactors.R with trapezoidal integration after transform
		log.BF.robust.naive4(R2, n, p_gamma), # 8 Robust via eq. 18 in BayesFactors.R with gsl implementation of 2F1
		log.BF.robust.safe1(R2, n, p_gamma),  # 9 Robust via eq. 18 in BayesFactors.R with log.hyperg_2F1
		log.BF.robust.safe2(R2, n, p_gamma)   # 10 Robust via eq. 18 in BayesFactors.R with gsl implementation of 2F1
	)

	D <- mutual_diff(vals_robust)
	return(D)
}


####################################################################################################


# Cake

check_cake <- function(n,p_gamma,R2)
{
	vals_cake <- c(
		log.BF.cake(R2, n, p_gamma), # cake via eq. 22 in BayesFactors.R
		BIC(n,p_gamma,R2)    # cake via eq. 22 in blma
	)

	D <- mutual_diff(vals_cake)
	return(D)
}

####################################################################################################

# Cake

check_ZS <- function(n,p_gamma,R2)
{
	quadrule = gaussLaguerre(1000, a = 0)

	vals_ZS <- c(
	    log_BF_Zellner_Siow_quad(n,p_gamma,R2),
 		log.BF.ZS.Laplace(R2,n,p_gamma),
 		log.BF.ZellnerSiow.MC(R2,n,p_gamma, N = 1.0E6),
 		log.BF.ZellnerSiow.quad(R2,n,p_gamma, quadrule),
 		log.BF.ZellnerSiow.integer(R2, n, p_gamma, method=c("Exact"),M=1000)
	)

	D <- mutual_diff(vals_ZS)
	return(D)
}


####################################################################################################

print("Test case 1")

n <- 100
R2 <- 0.25
p_gamma <- 10

print(check_hyper_g(n,p_gamma,R2))
print(check_hyper_g_n(n,p_gamma,R2))
print(check_robust(n,p_gamma,R2))
print(check_cake(n,p_gamma,R2))
print(check_ZS(n,p_gamma,R2))

####################################################################################################

print("Test case 2")

n <- 1000
R2 <- 0.25
p_gamma <- 10
p <- -100

print(check_hyper_g(n,p_gamma,R2))
print(check_hyper_g_n(n,p_gamma,R2))
print(check_robust(n,p_gamma,R2))
print(check_cake(n,p_gamma,R2))
print(check_ZS(n,p_gamma,R2))

####################################################################################################

print("Test case 3")

n <- 10000
R2 <- 0.25
p_gamma <- 10
p <- -100

print(check_hyper_g(n,p_gamma,R2))
print(check_hyper_g_n(n,p_gamma,R2))
print(check_robust(n,p_gamma,R2))
print(check_cake(n,p_gamma,R2))
print(check_ZS(n,p_gamma,R2))

####################################################################################################



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


normalize <- function(y, X) {
	n <- length(y)
	p <- ncol(X)

	mu.y <- mean(y)
	sigma2.y <- (n - 1) * var(y) / n
	vy <- (y - mu.y) / sqrt(sigma2.y)

	# Normalise covariates
	mX <- matrix(0, n, p)
	mu.x <- c()
	sigma2.x <- c()
	for (j in 1:p)
	{
	mu.x[j] <- mean(X[, j])
	sigma2.x[j] <- (n - 1) * var(X[, j]) / n
	mX[, j] <- (X[, j] - mu.x[j]) / sqrt(sigma2.x[j])
	}

	return(list(vy = vy, mX = mX, mu.y = mu.y, sigma2.y = sigma2.y, mu.x = mu.x, sigma2.x = sigma2.x))
}

if (TRUE)
{
	datname <- "UScrime"
	#datname <- "Kakadu"
	#datname <- "VietNamI"

	print(datname)

	if (datname=="UScrime") {

		library(blma); library(MASS)
		dat <- UScrime
		dat[,-c(2,ncol(UScrime))] <- log(dat[,-c(2,ncol(UScrime))])
		vy <- dat$y
		mX <- data.matrix(cbind(dat[1:15]))
		colnames(mX) <- c("log(AGE)","S","log(ED)","log(Ex0)","log(Ex1)",
			"log(LF)","log(M)","log(N)","log(NW)","log(U1)","log(U2)","log(W)",
			"log(X)","log(prison)","log(time)")
	}

	if (datname=="Kakadu") {
		dat = read.csv("Kakadu.csv")
		vy <- as.vector(dat$income)
		mX <- dat[,c(2:21,23)]
		mX <- model.matrix(~.,data=mX)[,-1]
		varnames <- colnames(mX)
	}


	if (datname=="VietNamI") {
		library(Ecdat)
		dat <- VietNamI
		# Get y vector and X matrix
		vy <- as.vector(dat$lnhhexp)
		mX <- dat[,-2]
		mX <- model.matrix(~.,data=mX)[,-1]
		varnames <- colnames(mX)
	}


	res <- normalize(vy,mX)
	vy <- res$vy
	mX <- res$mX

	###################################################################################################################

	# BLMA

	time1 <- system.time({blma_result_BIC              <- blma(vy, mX, prior="BIC"); })[3]
	time2 <- system.time({blma_result_ZE               <- blma(vy, mX, prior="ZE"); })[3]
	time3 <- system.time({blma_result_liang_g1         <- blma(vy, mX, prior="liang_g1"); })[3]
	time4 <- system.time({blma_result_liang_g2         <- blma(vy, mX, prior="liang_g2"); })[3]
	
	timeX <- system.time({ blma_result_liang_g_n_appell <- blma(vy, mX, prior="liang_g_n_appell"); })
	
	blma_result_liang_g_n_appell$vlogp[1] <- 0
	blma_result_liang_g_n_appell$vlogp[is.na(blma_result_liang_g_n_appell$vlogp)] <- -10000
	blma_result_liang_g_n_appell$vlogp[is.nan(blma_result_liang_g_n_appell$vlogp)] <- -10000
	
	library(blma)
	mGamma <- graycode(p)
	
	res <- CalculateProbabilities(blma_result_liang_g_n_appell$vlogp, mGamma)
	
	blma_result_liang_g_n_appell$vinclusion_prob <- res$var.prob
	
 
	
	time5 <- system.time({blma_result_liang_g_n_quad   <- blma(vy, mX, prior="liang_g_n_quad"); })[3]
	time6 <- system.time({blma_result_robust_bayarri1  <- blma(vy, mX, prior="robust_bayarri1"); })[3]
	time7 <- system.time({blma_result_robust_bayarri2  <- blma(vy, mX, prior="robust_bayarri2"); })[3]
	time8 <- system.time({blma_result_ZS_gauss_legendre  <- blma(vy, mX, prior="zellner_siow_gauss_legendre"); })[3]

	###################################################################################################################

	print("BayesVarSel")
	library( BayesVarSel )
	
	source("Bvs.r")

	a1 <- proc.time()[3]
	bvs.rob <- my.Bvs(formula="vy~.",fixed.cov=c("Intercept"),data=data.frame(vy=vy,mX=mX),
		prior.betas="Robust",prior.models="Constant",time.test= FALSE, priorprobs=NULL,n.keep=50000)
	b1 <- proc.time()[3]
	t1 <- b1-a1
	print(t1)

	a2 <- proc.time()[3]
	bvs.liang <- my.Bvs(formula="vy~.",fixed.cov=c("Intercept"),data=data.frame(vy=vy,mX=mX),
		prior.betas="Liangetal",prior.models="Constant",time.test= FALSE, priorprobs=NULL,n.keep=50000)
	b2 <- proc.time()[3]
	t2 <- b2-a2
	print(t2)

	a3 <- proc.time()[3]
	bvs.zs <- my.Bvs(formula="vy~.",fixed.cov=c("Intercept"),data=data.frame(vy=vy,mX=mX),prior.betas="ZellnerSiow",prior.models="Constant", priorprobs=NULL,n.keep=50000)
	b3 <- proc.time()[3]
	t3 <- b3-a3
	print(t3)

	bvs.rob.prob   <- bvs.rob$inclprob
	bvs.liang.prob <- bvs.liang$inclprob
	bvs.zs.prob    <- bvs.zs$inclprob

	###################################################################################################################

	print("BAS")
	library(BAS)

	a4 <- proc.time()[3]
	bas.bic =  bas.lm(vy~mX,prior="BIC",modelprior=uniform(), initprobs="uniform")
	b4 <- proc.time()[3]
	t4 <- b4-a4
	cat(t4,"\n")

	a5 <- proc.time()[3]
	bas.hg =  bas.lm(vy~mX,prior="hyper-g",modelprior=uniform(), initprobs="uniform")
	b5 <- proc.time()[3]
	t5 <- b5-a5
	cat(t5,"\n")

	a6 <- proc.time()[3]
	bas.hgl =  bas.lm(vy~mX,prior="hyper-g-laplace",modelprior=uniform(), initprobs="uniform")
	b6 <- proc.time()[3]
	t6 <- b6-a6
	cat(t6,"\n")

	a7 <- proc.time()[3]
	bas.hgn =  bas.lm(vy~mX,prior="hyper-g-n",modelprior=uniform(), initprobs="uniform")
	b7 <- proc.time()[3]
	t7 <- b7-a7
	cat(t7,"\n")

	a8 <- proc.time()[3]
	bas.ZSN =  bas.lm(vy~mX,prior="ZS-null",modelprior=uniform(), initprobs="uniform")
	b8 <- proc.time()[3]
	t8 <- b8-a8
	cat(t8,"\n")

	a9 <- proc.time()[3]
	bas.ZSF =  bas.lm(vy~mX,prior="ZS-full",modelprior=uniform(), initprobs="uniform")
	b9 <- proc.time()[3]
	t9 <- b9-a9
	cat(t9,"\n")


	bas.bic.prob <- bas.bic$probne0
	bas.hg.prob  <- bas.hg$probne0
	bas.hgl.prob <- bas.hgl$probne0
	bas.hgn.prob <- bas.hgn$probne0
	bas.ZSN.prob <- bas.ZSN$probne0
	bas.ZSF.prob <- bas.ZSF$probne0

	###################################################################################################################

	print("BMS")
	library(BMS)

	a10 <- proc.time()[3]
	res.bms <- bms(cbind(vy,mX), nmodel = 10000, mcmc = "enumerate", g="hyper=3", g.stats = TRUE, mprior="uniform")
	b10 <- proc.time()[3]
	t10 <- b10-a10
	cat(t10,"\n")

	fit        <- predict(res.bms, exact=TRUE)
	coef.bms   <- coef(res.bms,order.by.pip=FALSE )
	prob.hyper <- coef.bms[,1]

	###################################################################################################################

	p <- ncol(mX)
 	if (!exists("bvs.liang.prob"))
 	{
 		bvs.liang.prob <- rep(NA,p)
 	} 
 	if (!exists("bvs.rob.prob"))
	{
		bvs.rob.prob <- rep(NA,p)
	}
	 	
	tab <- round(100*cbind(
		bas.hg.prob[-1],                                   # HG  BAS  not safe
		bas.hgl.prob[-1],                                  # HG  BAS  Laplace
		bvs.liang.prob,                                # HG  BVS  not safe 
		prob.hyper,                                    # HG  BMS  not safe
		blma_result_liang_g1$vinclusion_prob,          # HG  BLMA not safe
		blma_result_liang_g2$vinclusion_prob,          # HG  BLMA safe
		
		bas.hgn.prob[-1],                                  # HGN BAS Laplace
		blma_result_liang_g_n_appell$vinclusion_prob,  # HGN BLMA Exact
		blma_result_liang_g_n_quad$vinclusion_prob,    # HGN BLMA Quad
		
		bvs.rob.prob,                                  # ROB BVS Rob
#		blma_result_robust_bayarri1$vinclusion_prob,   # ROB BLMA exact not-safe
		blma_result_robust_bayarri2$vinclusion_prob,   # ROB BLMA exact safe
		
		bas.ZSN.prob[-1],                                  # ZS BAS Laplace
#		bas.ZSF.prob[-1],                                  # ZS BAS Laplace
		blma_result_ZS_gauss_legendre$vinclusion_prob  # ZS BLMA Quad
	),2)
	
	print(tab)
	
	times <- c(t5,t6,t2,t10,time3,time4,t7,0,time5,t1,time7,t8,time8)
	
	save.image(paste("Results_",datname,".Rdata",sep=""))
	
	
		
	source("BAS_functions.R")
	vR2 <- blma_result_BIC$vR2
	vp <- blma_result_BIC$vp_gamma
	n  <- length(vy)
	p <- ncol(mX)
	
	library(blma)
	mGamma <- graycode(p,0)
	
	logBF <- LogBF_GonN_Laplace(vR2, n, vp, alpha=3)$logmarg
	logBF <- LogBF_ZS_Laplace(vR2, n, vp)$logmarg
	logBF <- logBF - max(logBF)
	mpostprod <- exp(logBF)/sum(exp(logBF))
	
	vpostprod <- t(mGamma)%*%mpostprod
	
	
	
}


if (FALSE)
{
	datname <- "Communities"

	datname <- "building"
	
	datname <- "eye" 

	print(datname)

	if (datname=="Communities")
	{
		load("comData.Rdata")
		sum.na <- function(x) {
			sum(is.na(x))
		}

		inds <- which(apply(X, 2, sum.na) == 0)
		mX <- X[, inds]

		i <- 1
		vy <- Y[, i]
		vy <- vy - mean(vy)

		inds <- which(colnames(mX) %in% c("ownHousQrange", "rentUpperQ"))
		mX <- mX[, -inds]
	}

	if (datname=="building") {
		data(building)
		vy <- building$OUTPUTS
		mX <- as.matrix(building[,1:107])
	}
	
	if (FALSE) {
		library(tidyverse)
		library(mlbench)
		data(BostonHousing2)
		
		tib <- as.tibble(BostonHousing2)
		vy <- tib$medv
		dat <- tib[,-c(1,5,10)]
		
		
		form <- (~ -1+.^2)
		mX <- model.matrix(form, data = dat)
	}
	
	if (datname=="eye")
	{
		load("eyedata.Rdata")
		vy = y
		mX = x
	}

 
	res <- normalize(vy, mX)
	vy <- res$vy
	mX <- res$mX
	
	n <- length(vy)
	p <- ncol(mX)

	###################################################################################################################

 	N <- 10000
 
    if (TRUE) {
	    #blma_result_liang_g_n_appell <- blma(vy, mX, prior="liang_g_n_appell")
	    time5 <- system.time({blma_result_liang_g_n_quad   <- sampler(N, vy, mX, prior="liang_g_n_quad",modelprior="uniform"); })[3]
	    save(time5,file="time5.Rdata")
	    time6 <- system.time({blma_result_robust_bayarri1  <- sampler(N, vy, mX, prior="robust_bayarri1",modelprior="uniform"); })[3]
	    save(time6,file="time6.Rdata")
	    time7 <- system.time({blma_result_robust_bayarri2  <- sampler(N, vy, mX, prior="robust_bayarri2",modelprior="uniform"); })[3]
	    save(time7,file="time7.Rdata")
	    time8 <- system.time({blma_result_ZS_gauss_legendre  <- sampler(N, vy, mX, prior="zellner_siow_gauss_legendre",modelprior="uniform"); })[3]
		save(time8,file="time8.Rdata")
		
		
		time1 <- system.time({blma_result_BIC              <- sampler(N, vy, mX, prior="BIC",modelprior="uniform"); })[3]
		save(time1,file="time1.Rdata")
		time2 <- system.time({blma_result_ZE               <- sampler(N, vy, mX, prior="ZE",modelprior="uniform"); })[3]
		save(time2,file="time2.Rdata")
		time3 <- system.time({blma_result_liang_g1         <- sampler(N, vy, mX, prior="liang_g1",modelprior="uniform"); })[3]
		save(time3,file="time3.Rdata")
		time4 <- system.time({blma_result_liang_g2         <- sampler(N, vy, mX, prior="liang_g2",modelprior="uniform"); })[3]
		save(time4,file="time4.Rdata")
	} else {
		time5 <- system.time({blma_result_liang_g_n_quad   <- sampler(N, vy, mX, prior="liang_g_n_quad",modelprior="beta-binomial", modelpriorvec=c(1,p)); })[3]
		save(time5,file="time5.Rdata")
		time6 <- system.time({blma_result_robust_bayarri1  <- sampler(N, vy, mX, prior="robust_bayarri1",modelprior="beta-binomial", modelpriorvec=c(1,p)); })[3]
		save(time6,file="time6.Rdata") 
		time7 <- system.time({blma_result_robust_bayarri2  <- sampler(N, vy, mX, prior="robust_bayarri2",modelprior="beta-binomial", modelpriorvec=c(1,p)); })[3]
		save(time7,file="time7.Rdata")
		time8 <- system.time({blma_result_ZS_gauss_legendre  <- sampler(N, vy, mX, prior="zellner_siow_gauss_legendre",modelprior="beta-binomial", modelpriorvec=c(1,p)); })[3]
		save(time8,file="time8.Rdata")
		
		
		time1 <- system.time({blma_result_BIC              <- sampler(N, vy, mX, prior="BIC",modelprior="beta-binomial", modelpriorvec=c(1,p)); })[3]
		save(time1,file="time1.Rdata")
		time2 <- system.time({blma_result_ZE               <- sampler(N, vy, mX, prior="ZE",modelprior="beta-binomial", modelpriorvec=c(1,p)); })[3]
		save(time2,file="time2.Rdata")
		time3 <- system.time({blma_result_liang_g1         <- sampler(N, vy, mX, prior="liang_g1",modelprior="beta-binomial", modelpriorvec=c(1,p)); })[3]
		save(time3,file="time3.Rdata")
		time4 <- system.time({blma_result_liang_g2         <- sampler(N, vy, mX, prior="liang_g2",modelprior="beta-binomial", modelpriorvec=c(1,p)); })[3]
		save(time4,file="time4.Rdata")
	}
	
	##########################################################################################################

	if (FALSE) {
		# Note: All of these seem to return infinite Bayes factors.

		print("BayesVarSel")
		library( BayesVarSel )

		source("GibbsBvs.R")

		a1 <- proc.time()[3]
		bvs.rob <- my.GibbsBvs(formula="vy~.",data=data.frame(vy=vy,mX=mX),
			prior.betas="Robust",prior.models="Constant",
			n.iter = 10000, init.model = "Null", n.burnin = 500, n.thin = 1,
			time.test = TRUE, priorprobs = NULL, seed = 1)

		b1 <- proc.time()[3]
		t1 <- b1-a1
		print(t1)

		a2 <- proc.time()[3]
		bvs.liang <- my.GibbsBvs(formula="vy~.",data=data.frame(vy=vy,mX=mX),
			prior.betas="Liangetal",prior.models="Constant",
			n.iter = 10000, init.model = "Null", n.burnin = 500, n.thin = 1,
			time.test = TRUE, priorprobs = NULL, seed = 1)
		b2 <- proc.time()[3]
		t2 <- b2-a2
		print(t2)

		a3 <- proc.time()[3]
		bvs.zs <- my.GibbsBvs(formula="vy~.",data=data.frame(vy=vy,mX=mX),
			prior.betas="ZellnerSiow",prior.models="Constant",
			n.iter = 10000, init.model = "Null", n.burnin = 500, n.thin = 1,
			time.test = TRUE, priorprobs = NULL, seed = 1)
		b3 <- proc.time()[3]
		t3 <- b3-a3
		print(t3)

		bvs.rob.prob   <- bvs.rob$inclprob
		bvs.liang.prob <- bvs.liang$inclprob
		bvs.zs.prob    <- bvs.zs$inclprob
	}

	##########################################################################################################

	print("BAS")
	library(BAS)

	a4 <- proc.time()[3]
	bas.bic =  bas.lm(vy~-1+mX, n.models=10000, method="MCMC",
		prior="BIC", modelprior=uniform(), MCMC.iterations=N )

	b4 <- proc.time()[3]
	t4 <- b4-a4
	cat(t4,"\n")

	a5 <- proc.time()[3]
	bas.hg =  bas.lm(vy~-1+mX, n.models=10000, method="MCMC",
		prior="hyper-g", modelprior=uniform(), MCMC.iterations=N )
	b5 <- proc.time()[3]
	t5 <- b5-a5
	cat(t5,"\n")

	a6 <- proc.time()[3]
	bas.hgl =  bas.lm(vy~-1+mX, n.models=10000, method="MCMC",
		prior="hyper-g-laplace", modelprior=uniform(),  MCMC.iterations=N )
	b6 <- proc.time()[3]
	t6 <- b6-a6
	cat(t6,"\n")

	a7 <- proc.time()[3]
	bas.hgn =  bas.lm(vy~-1+mX, n.models=10000, method="MCMC",
		prior="hyper-g-n", modelprior=uniform(),  MCMC.iterations=N )
	b7 <- proc.time()[3]
	t7 <- b7-a7
	cat(t7,"\n")

	a8 <- proc.time()[3]
	bas.ZSN = bas.lm(vy~-1+mX, n.models=10000, method="MCMC",
		prior="ZS-null", modelprior=uniform(), MCMC.iterations=N )
	b8 <- proc.time()[3]
	t8 <- b8-a8
	cat(t8,"\n")

	a9 <- proc.time()[3]
	bas.ZSF = bas.lm(vy~-1+mX, n.models=10000, method="MCMC",
		prior="JZS", modelprior=uniform(),  MCMC.iterations=N )
	b9 <- proc.time()[3]
	t9 <- b9-a9
	cat(t9,"\n")

	###################################################################################################################

	# I tried BMS but for the hyper-g prior, but it returned NAs.

	print("BIC results")
	print(cbind( apply( blma_result_BIC$mGamma, 2, mean),  bas.bic$probne0.MCMC))

	print("hyper-g results")
	print(cbind(
				apply( blma_result_liang_g2$mGamma, 2, mean),
				bas.hg$probne0.MCMC,
				bas.hgl$probne0.MCMC))

	print("hyper g/n results")
	print(cbind( apply( blma_result_liang_g_n_quad$mGamma, 2, mean), bas.hgn$probne0.MCMC))

	print("ZS results")
	print(cbind( apply( blma_result_ZS_gauss_legendre$mGamma, 2, mean), bas.ZSN$probne0.MCMC, bas.ZSF$probne0.MCMC))
	
	save.image(paste("Results_",datname,".Rdata",sep=""))
}



