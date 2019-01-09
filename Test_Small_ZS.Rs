

library(blma)
source("./Common/BayesFactorsZS.R")

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
print(check_ZS(n,p_gamma,R2))

print("Test case 2")
n <- 1000
R2 <- 0.25
p_gamma <- 10
print(check_ZS(n,p_gamma,R2))

print("Test case 3")
n <- 10000
R2 <- 0.25
p_gamma <- 10
print(check_ZS(n,p_gamma,R2))

