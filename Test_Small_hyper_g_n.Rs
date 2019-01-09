

library(blma)

source("./Common/BayesFactorsGonN.R",chdir=TRUE)
 
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
	
	meths <- c(
		"hyper-g/n via eq. 14 in package blma using gsl for appell",
		"hyper-g/n via eq. 13 in package blma with Gauss-Legendre",
		"hyper-g/n via old appell fn approx.",
		"hyper-g/n via eq. 14 in BayesFactors.R using gsl for appel",
		"hyper-g/n via eq. 13 in BayesFactors.R with Gauss-Legendre",
		"hyper-g/n via eq. 13 in BayesFactors.R with Laplace approx.")
	
	names(vals_hyper_g_n) <- meths

	D <- mutual_diff(vals_hyper_g_n)
	
	return(list(D=D,vals_hyper_g_n=vals_hyper_g_n,meths=meths))
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
	
	meths <- c(
		"1 Robust via eq. 17 in package blma",
		"2 Robust via eq. 18 in package blma",
		"3 Robust via eq. 17 in BayesFactors.R with gsl implementation of 2F1",
		"4 Robust via eq. 17 in BayesFactors.R with BAS implementation of 2F1",
		"5 Robust via eq. 17 in BayesFactors.R with BAS implementation of 2F1 using Laplace method",
		"6 Robust via eq. 17 in BayesFactors.R with trapezoidal integration",
		"7 Robust via eq. 17 in BayesFactors.R with trapezoidal integration after transform",
		"8 Robust via eq. 18 in BayesFactors.R with gsl implementation of 2F1",
		"9 Robust via eq. 18 in BayesFactors.R with log.hyperg_2F1",
		"10 Robust via eq. 18 in BayesFactors.R with gsl implementation of 2F1")
		
	names(vals_robust) <- meths
	
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
 		log.BF.ZellnerSiow.MC(R2,n,p_gamma, N = 1.0E7),
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
print(check_hyper_g_n(n,p_gamma,R2))


print("Test case 2")

n <- 1000
R2 <- 0.25
p_gamma <- 10
print(check_hyper_g_n(n,p_gamma,R2))

print("Test case 3")

n <- 10000
R2 <- 0.25
p_gamma <- 10
print(check_hyper_g_n(n,p_gamma,R2))

####################################################################################################
