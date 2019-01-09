
library(blma)
source("./Common/BayesFactors.R")

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

print("Test case 1")

n <- 100
R2 <- 0.25
p_gamma <- 10
print(check_robust(n,p_gamma,R2))

print("Test case 2")
n <- 1000
R2 <- 0.25
p_gamma <- 10
print(check_robust(n,p_gamma,R2))

print("Test case 3")
n <- 10000
R2 <- 0.25
p_gamma <- 10
print(check_robust(n,p_gamma,R2))

