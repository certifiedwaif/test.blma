
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

print("Test case 1")

n <- 100
R2 <- 0.25
p_gamma <- 10
print(check_hyper_g(n,p_gamma,R2))


print("Test case 2")

n <- 1000
R2 <- 0.25
p_gamma <- 10
print(check_hyper_g(n,p_gamma,R2))

print("Test case 3")

n <- 10000
R2 <- 0.25
p_gamma <- 10
p <- -100

print(check_hyper_g(n,p_gamma,R2))