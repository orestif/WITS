############################################
#		SALMONELLA WITS MODEL 
#	SOLUTION OF THE MOMENTS EQUATION
############################################

library (expm)

library(foreach)
library(doParallel)
registerDoParallel(cores=16)

library(Rcpp)

setwd("~/Documents/Work/server/Salmonella/ChrisWITS")

# Version notes:
# - 2016/03/05: the gillespie simulation C++ code is now in file "gillespie.cpp", and has been corrected to exclude the first event after t_end and remove the cap N_max. 

# =============================================================================================================================
# 					 			MAIN FUNCTIONS TO USE 
# =============================================================================================================================

# ----------------------------- COMPUTE THE PREDICTED MOMENTS ----------------------------

# Break down the calculation of moments into a series of smaller steps 
# - t: time point to cimpute (starting from t=0)
# - par: named vector of 8 parameter values (cL,cS,eL,eS,kL,kS,rL,rS)
# - M0: vector of initial conditions: value of the 9 moments at t=0
# - interval: legnth of intervals to break down the calculation over [0,t]
# - solver is one of the functions below
# - h: time step for integration
# - exp.met: method for expm()
# Default values based on benchmarking tests

WITS.moment.sol.steps <- function(t,par,M0,interval=6,solver=WITS.moment.sol.4,h=1,exp.met="AlMohy-Hi09"){
	t.steps <- unique(c(seq(0,t,interval),t))
	t0 <- 0
	for(t1 in t.steps[-1]){
		M0 <- solver(t1-t0,par,M0,h,exp.met)
		t0 <- t1
	}
	return(M0)
}


# ----------------------------- KL Divergence -----------------------------

# Calculate the KL divergence between two multivariate normal distributions 
# with mean vectors mu0 and mu1, and cov matrices cov0 and cov1
KL.div <- function(mu0,cov0,mu1,cov1){
	k <- length(mu0)
	inv.C1 <- solve(cov1)
	as.numeric(sum(diag(inv.C1 %*% cov0)) + (mu1-mu0) %*% inv.C1 %*% (mu1-mu0) - k + log(det(cov1)/det(cov0)))/2
}

# Calculate the KL divergence from a predicted distribution to an observed sample, both characterised by their vectors of 9 moments
# - pred: predicted moments
# - obs: observed moments
# - sub: vector specifying which organs to include, e.g. c(1,2) for B and L. Default 1:3
KL.div.M2M <- function(pred,obs,sub=1:3){
	# Predicted moments
	mu0 <- pred[sub]
	cov0 <- matrix(pred[c(4,7,8,7,5,9,8,9,6)],3)[sub,sub]
	# Observed moments
	mu1 <- obs[sub]
	cov1 <- matrix(obs[c(4,7,8,7,5,9,8,9,6)],3)[sub,sub]
	# KL divergence
	KL.div(mu0,cov0,mu1,cov1)
}

# ----------------------------- Hellinger distance -----------------------------

# Calculate the Hellinger distance between two [multinormal] distributions with mean vectors mu0 and mu1, and cov matrices cov0 and cov1
Hell.dist <- function(mu0,cov0,mu1,cov1){
	P <- (cov0+cov1)/2
	m <- mu0-mu1
	as.numeric(m %*% solve(P) %*% m + log(det(P)/sqrt(det(cov0)*det(cov1)))/2)/8
}

# Calculate the Hell distance from an observed sample to a predicted distribution, both characterised by their vectors of 9 moments
# - obs: observed moments
# - pred: predicted moments
# - sub: vector specifying which organs to include, e.g. c(1,2) for B and L. Default 1:3
Hell.dist.M2M <- function(obs,pred,sub=1:3){
	# Observed moments
	mu0 <- obs[sub]
	cov0 <- matrix(obs[c(4,7,8,7,5,9,8,9,6)],3)[sub,sub]
	# Predicted moments
	mu1 <- pred[sub]
	cov1 <- matrix(pred[c(4,7,8,7,5,9,8,9,6)],3)[sub,sub]
	# Hell distance
	Hell.dist(mu0,cov0,mu1,cov1)
}



# ---------------------------- GILLESPIE SIMULATIONS -----------------------------------------

# Use WITS Gillespie model written in C++ 
#sourceCpp('gimh_cpp/gimh_simple.cpp')
sourceCpp('Moments/gillespie.cpp')

# Run a series of simulations in parallel
# In order to allow random variations in initial conditions, init.fun has to be a function that generates a vector of initial conditions when called
# Return a matrix with n.sim rows and 6 columns: nB(0), nL(0), nS(0), nB(t), nL(t), nS(t)
WITS.gillespie.sim <- function(t,par,init.fun,n.sim){
	foreach(i=1:n.sim, .combine=rbind) %dopar% {
		init <- init.fun()
		c(init,single_wits_gillespie_r(0,t,init,par))
	}
}

# Calculate the means, variances and covariances in blood, liver and spleen from a matrix of simulated values.
WITS.gillespie.moments <- function(sim){
	x <- c(mean(sim[,4]),mean(sim[,5]),mean(sim[,6]),var(sim[,4]),var(sim[,5]),var(sim[,6]),cov(sim[,4],sim[,5]),cov(sim[,4],sim[,6]),cov(sim[,5],sim[,6]))
	names(x) <- mom.names
	x
}


# ------------------------------- OTHER USEFUL FUNCTIONS ------------------------------------


# Takes the vector of variance-covriance M2 = (NB, NL, NS, NB*NL, NB*NS, NL*NS) and returns the covariance matrix
M2cov <- function (M){
	r <- c(1,4,5,4,2,6,5,6,3)
	matrix(M[r],3,3)
}



# Substitute individual parameter values by name
# - x and rep are named vectors.
replace.par <- function(x,rep)
{
	if(length(rep)==0) return(x)
	pos <- sapply(names(rep),function(n) which(names(x)==n))
	replace(x,pos,rep)
}

# Labels for the 9 moments (expectation and variances)
mom.names <- c('E.NB','E.NL','E.NS','V.NB','V.NL','V.NS','V.NB.NL','V.NB.NS','V.NL.NS')




# =============================================================================================================================
# 					 			SUBORDINATE FUNCTIONS 
# =============================================================================================================================


# ======================================= Second-moment solution ===============================================================

# All the functions below calculate the 9 moments for a given set of parameter values:
# - t is the time (starting from t=0)
# - par must be a named vector (cL,cS,eL,eS,kL,kS,rL,rS)
# - M0 is a vector of the initial means,variances and covariances: (E NB, E NL, E NS,V NB, V NL, V NS, V NB*NL, V NB*NS, V NL*NS)
# - h is the step size for the numerical integration
# Returns a single vector of length 9: c(M.1(t),M.2(t))

# -------------------------- Version 1 -------------------------------
# Integrates the vector ecba(s) %*% M1(0) by linear approximation
WITS.moment.sol.1 <- function(t,par,M0,h=1e-3,met="Higham08.b"){with(as.list(par),{
	# First moment: M.1(t) = (E[NB],E[NL],E[NS]) is solution of M.1'(t) = A * M.1
	A <- matrix(c(-cL-cS,cL,cS,eL,rL-kL-eL,0,eS,0,rS-kS-eS), 3,3)
	# Solution
	M.1 <- expAtv(A,M0[1:3],t)$eAtv
	names(M.1) <- c('E.NB','E.NL','E.NS')
	# Second moment: M.2(t) = (V(NB),V(NL),V(NS),V(NB,NL),V(NB,NS),V(NL,NS)) is solution of M.2' = B*M.1 + C*M.2
	B <- matrix(c(cL+cS,cL,cS,-cL,-cS,0, eL,rL+kL+eL,0,-eL,0,0, eS,0,rS+kS+eS,0,-eS,0), 6,3)
	C <- matrix(c(
		-2*cL-2*cS, 0, 0, 2*eL, 2*eS, 0,
		0, 2*rL-2*kL-2*eL, 0, 2*cL, 0, 0,
		0, 0 ,2*rS-2*kS-2*eS,0,2*cS,0,
		cL,eL,0,-cL-cS-eL+rL-kL,0,eS,
		cS,0,eS,0,-cL-cS-eS+rS-kS,eL,
		0,0,0,cL,cS,-eL-eS-kL-kS+rL+rS
	),6,6,byrow=T)
	# ecba.M1(s) returns the vector exp(-s*C) %*% B %*% [exp(s*A) %*% M.1(0)], with length 6
	ecba.M1 <- function(s) {as.numeric(expm(-s*C,m=met) %*% (B %*% expm(s*A,m=met) %*% M0[1:3]))}
	# - Step 1 (long): Calculate ecba.M1 for s in [0,t] => 6-row matrix
	fs <- matrix(sapply(seq(0,t,h), ecba.M1),6)
	# - Step 2 (short): summation of each row
	n <- ncol(fs)
	int.ecba.M1 <- apply(fs,1,function(fx){ sum(fx[-n]*h/2)+sum(fx[-1]*h/2) })
	# Calculate M2(t)
	M.2 <- expAtv(C, M0[4:9] + int.ecba.M1, t)$eAtv
	names(M.2) <- c('V.NB','V.NL','V.NS','V.NB.NL','V.NB.NS','V.NL.NS')
	return(c(M.1,M.2))
})}


# -------------------------- Version 2 -------------------------------
# Integrates the vector ecba(s) %*% M1(0) using Simpson cubic approximation
WITS.moment.sol.2 <- function(t,par,M0,h=1e-3,met="Higham08.b"){with(as.list(par),{
	# First moment: M.1(t) = (E[NB],E[NL],E[NS]) is solution of M.1'(t) = A * M.1
	A <- matrix(c(-cL-cS,cL,cS,eL,rL-kL-eL,0,eS,0,rS-kS-eS), 3,3)
	# Solution
	M.1 <- expAtv(A,M0[1:3],t)$eAtv
	names(M.1) <- c('E.NB','E.NL','E.NS')
	# Second moment: M.2(t) = (V(NB),V(NL),V(NS),V(NB,NL),V(NB,NS),V(NL,NS)) is solution of M.2' = B*M.1 + C*M.2
	B <- matrix(c(cL+cS,cL,cS,-cL,-cS,0, eL,rL+kL+eL,0,-eL,0,0, eS,0,rS+kS+eS,0,-eS,0), 6,3)
	C <- matrix(c(
		-2*cL-2*cS, 0, 0, 2*eL, 2*eS, 0,
		0, 2*rL-2*kL-2*eL, 0, 2*cL, 0, 0,
		0, 0 ,2*rS-2*kS-2*eS,0,2*cS,0,
		cL,eL,0,-cL-cS-eL+rL-kL,0,eS,
		cS,0,eS,0,-cL-cS-eS+rS-kS,eL,
		0,0,0,cL,cS,-eL-eS-kL-kS+rL+rS
	),6,6,byrow=T)
	# ecba.M1(s) returns the vector exp(-s*C) %*% B %*% [exp(s*A) %*% M.1(0)], with length 6
	ecba.M1 <- function(s) {as.numeric(expm(-s*C,m=met) %*% (B %*% expm(s*A,m=met) %*% M0[1:3]))}
	# - Step 1 (long): Calculate ecba.M1 for s in [0,t] => 6-row matrix
	fs <- matrix(sapply(seq(0,t,h), ecba.M1),6)
	# - Step 2 (short): Simpson 3/8 integration
	int.ecba.M1 <- apply(fs,1,function(fx){
		coeff <- c(1,rep(c(3,3,2),length.out=length(fx)-2),1)
		3*h/8*sum(fx*coeff)
	})
	# Calculate M2(t)
	M.2 <- expAtv(C, M0[4:9] + int.ecba.M1, t)$eAtv
	names(M.2) <- c('V.NB','V.NL','V.NS','V.NB.NL','V.NB.NS','V.NL.NS')
	return(c(M.1,M.2))
})}


# -------------------------- Version 3 -------------------------------
# This version uses pracma::quadv to integrate the vector.
WITS.moment.sol.3 <- function(t,par,M0,met="Higham08.b"){with(as.list(par),{
	require(pracma)
	# First moment: M.1(t) = (E[NB],E[NL],E[NS]) is solution of M.1'(t) = A * M.1
	A <- matrix(c(-cL-cS,cL,cS,eL,rL-kL-eL,0,eS,0,rS-kS-eS), 3,3)
	# Solution
	M.1 <- expAtv(A,M0[1:3],t)$eAtv
	names(M.1) <- c('E.NB','E.NL','E.NS')
	# Second moment: M.2(t) = (V(NB),V(NL),V(NS),V(NB,NL),V(NB,NS),V(NL,NS)) is solution of M.2' = B*M.1 + C*M.2
	B <- matrix(c(cL+cS,cL,cS,-cL,-cS,0, eL,rL+kL+eL,0,-eL,0,0, eS,0,rS+kS+eS,0,-eS,0), 6,3)
	C <- matrix(c(
		-2*cL-2*cS, 0, 0, 2*eL, 2*eS, 0,
		0, 2*rL-2*kL-2*eL, 0, 2*cL, 0, 0,
		0, 0 ,2*rS-2*kS-2*eS,0,2*cS,0,
		cL,eL,0,-cL-cS-eL+rL-kL,0,eS,
		cS,0,eS,0,-cL-cS-eS+rS-kS,eL,
		0,0,0,cL,cS,-eL-eS-kL-kS+rL+rS
	),6,6,byrow=T)
	# ecba returns the vector exp(-s*C) B [exp(s*A) M.1(0)]
	ecba <- function(s) {as.numeric(expm::expm(-s*C,m=met) %*% (B %*% expm::expm(s*A,m=met) %*% M0[1:3]))}
	int.ecba <- quadv(ecba,0,t)$Q
	M.2 <- expAtv(C,int.ecba+M0[4:9],t)$eAtv
	names(M.2) <- c('V.NB','V.NL','V.NS','V.NB.NL','V.NB.NS','V.NL.NS')
	detach(package:pracma)
	return(c(M.1,M.2))
})}



# -------------------------- Version 4 -------------------------------
# This version integrates the matrix with a linear approximation
WITS.moment.sol.4 <- function(t,par,M0,h=1e-3,met="Higham08.b"){with(as.list(par),{
	# First moment: M.1(t) = (E[NB],E[NL],E[NS]) is solution of M.1'(t) = A * M.1
	A <- matrix(c(-cL-cS,cL,cS,eL,rL-kL-eL,0,eS,0,rS-kS-eS), 3,3)
	# Solution
	M.1 <- expAtv(A,M0[1:3],t)$eAtv
	names(M.1) <- c('E.NB','E.NL','E.NS')
	# Second moment: M.2(t) = (V(NB),V(NL),V(NS),V(NB,NL),V(NB,NS),V(NL,NS)) is solution of M.2' = B*M.1 + C*M.2
	B <- matrix(c(cL+cS,cL,cS,-cL,-cS,0, eL,rL+kL+eL,0,-eL,0,0, eS,0,rS+kS+eS,0,-eS,0), 6,3)
	C <- matrix(c(
		-2*cL-2*cS, 0, 0, 2*eL, 2*eS, 0,
		0, 2*rL-2*kL-2*eL, 0, 2*cL, 0, 0,
		0, 0 ,2*rS-2*kS-2*eS,0,2*cS,0,
		cL,eL,0,-cL-cS-eL+rL-kL,0,eS,
		cS,0,eS,0,-cL-cS-eS+rS-kS,eL,
		0,0,0,cL,cS,-eL-eS-kL-kS+rL+rS
	),6,6,byrow=T)
	# ecba returns the 6x3 matrix exp(-s*C) %*% B %*% exp(s*A) flattened into a vector of length 18
	ecba <- function(s) {as.double((expm(-s*C,m=met) %*% B) %*% expm(s*A,m=met))}
	# - Step 1 (time-consuming): Calculate ecba(s) for s in [0,t] => matrix with 18 rows
	fs <- sapply(seq(0,t,h), ecba)
	# - Step 2: summation
	n <- ncol(fs)
	int.ecba <- apply(fs,1,function(fx){ sum(fx[-n]*h/2)+sum(fx[-1]*h/2) })
	ecba.M.1 <- matrix(int.ecba,nrow=6) %*% M0[1:3]
	M.2 <- as.numeric(expm(t*C,m=met) %*% (M0[4:9] + ecba.M.1))
	names(M.2) <- c('V.NB','V.NL','V.NS','V.NB.NL','V.NB.NS','V.NL.NS')
	return(c(M.1,M.2))
})}


# -------------------------- Version 5 -------------------------------
# This version integrates the matrix, using Simpson's cubic approximation
WITS.moment.sol.5 <- function(t,par,M0,h=1e-3,met="Higham08.b"){with(as.list(par),{
	# First moment: M.1(t) = (E[NB],E[NL],E[NS]) is solution of M.1'(t) = A * M.1
	A <- matrix(c(-cL-cS,cL,cS,eL,rL-kL-eL,0,eS,0,rS-kS-eS), 3,3)
	# Solution
	M.1 <- expAtv(A,M0[1:3],t)$eAtv
	names(M.1) <- c('E.NB','E.NL','E.NS')
	# Second moment: M.2(t) = (V(NB),V(NL),V(NS),V(NB,NL),V(NB,NS),V(NL,NS)) is solution of M.2' = B*M.1 + C*M.2
	B <- matrix(c(cL+cS,cL,cS,-cL,-cS,0, eL,rL+kL+eL,0,-eL,0,0, eS,0,rS+kS+eS,0,-eS,0), 6,3)
	C <- matrix(c(
		-2*cL-2*cS, 0, 0, 2*eL, 2*eS, 0,
		0, 2*rL-2*kL-2*eL, 0, 2*cL, 0, 0,
		0, 0 ,2*rS-2*kS-2*eS,0,2*cS,0,
		cL,eL,0,-cL-cS-eL+rL-kL,0,eS,
		cS,0,eS,0,-cL-cS-eS+rS-kS,eL,
		0,0,0,cL,cS,-eL-eS-kL-kS+rL+rS
	),6,6,byrow=T)
	# ecba returns the 6x3 matrix exp(-s*C) %*% B %*% exp(s*A) flattened into a vector of length 18
	ecba <- function(s) {as.numeric(expm(-s*C,m=met) %*% (B %*% expm(s*A,m=met)))}
	# Simpson 3/8 integration
	# - Step 1 (time-consuming): Calculate ecba(s) for s in [0,t] => matrix with 18 rows
	fs <- sapply(seq(0,t,h), ecba)
	# - Step 2: summation using Simpon's "3/8" rule (cubic approximation)
	int.ecba <- apply(fs,1,function(fx){
		coeff <- c(1,rep(c(3,3,2),length.out=length(fx)-2),1)
		3*h/8*sum(fx*coeff)
	})
	ecba.M.1 <- matrix(int.ecba,nrow=6) %*% M0[1:3]
	M.2 <- as.numeric(expm(t*C,m=met) %*% (M0[4:9] + ecba.M.1))
	names(M.2) <- c('V.NB','V.NL','V.NS','V.NB.NL','V.NB.NS','V.NL.NS')
	return(c(M.1,M.2))
})}







# ---------------------------------------- For testing purposes only ----------------------------------------------
# Return a matrix with the integrand values (fs) calculated at every time step.
WITS.moment.sol.fs <- function(t,par,M0,h=1e-3,met="Higham08.b"){with(as.list(par),{
	# First moment: M.1(t) = (E[NB],E[NL],E[NS]) is solution of M.1'(t) = A * M.1
	A <- matrix(c(-cL-cS,cL,cS,eL,rL-kL-eL,0,eS,0,rS-kS-eS), 3,3)
	# Solution
	M.1 <- expAtv(A,M0[1:3],t)$eAtv
	names(M.1) <- c('E.NB','E.NL','E.NS')
	# Second moment: M.2(t) = (V(NB),V(NL),V(NS),V(NB,NL),V(NB,NS),V(NL,NS)) is solution of M.2' = B*M.1 + C*M.2
	B <- matrix(c(cL+cS,cL,cS,-cL,-cS,0, eL,rL+kL+eL,0,-eL,0,0, eS,0,rS+kS+eS,0,-eS,0), 6,3)
	C <- matrix(c(
		-2*cL-2*cS, 0, 0, 2*eL, 2*eS, 0,
		0, 2*rL-2*kL-2*eL, 0, 2*cL, 0, 0,
		0, 0 ,2*rS-2*kS-2*eS,0,2*cS,0,
		cL,eL,0,-cL-cS-eL+rL-kL,0,eS,
		cS,0,eS,0,-cL-cS-eS+rS-kS,eL,
		0,0,0,cL,cS,-eL-eS-kL-kS+rL+rS
	),6,6,byrow=T)
	# ecba returns the 6x3 matrix exp(-s*C) %*% B %*% exp(s*A) flattened into a vector of length 18
	ecba <- function(s) {as.numeric((expm(-s*C,m=met) %*% B) %*% expm(s*A,m=met))}
	# - Step 1 (time-consuming): Calculate ecba(s) for s in [0,t] => matrix with 18 rows
	fs <- sapply(seq(0,t,h), ecba)
	return(fs)
})}
