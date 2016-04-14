#############################
#   SALMONELLA WITS MODEL
#     MOMENTS EQUATION
#          TESTS
#############################

# PARAMETER INFERENCE USING Kullbeck-Leibler Divergence

source("Moments/WITS_moments.R")

library(dplyr)
library(ggplot2)
library(tidyr)
library(powell)
library(ellipse)

# ========================================= KL Divergence =============================================================

# Calculate the KL divergence between two distributions with mean vectors mu0 and mu1, and cov matrices cov0 and cov1
# Note: this 
KL.div <- function(mu0,cov0,mu1,cov1){
	k <- length(mu0)
	inv.C1 <- solve(cov1)
	as.numeric(sum(diag(inv.C1 %*% cov0)) + (mu1-mu0) %*% inv.C1 %*% (mu1-mu0) - k + log(det(cov1)/det(cov0)))/2
}

# Calculate the KL divergence from a predicted distribution to an observed sample, both characterised by their vectors of 9 moments
# - pred: predicted moments
# - obs: observed moments
# - sub: vector specifying which organs to include, e.g. c(1,2) for B and L. Default 1:3
KL.div.p2o <- function(pred,obs,sub=1:3){
	# Predicted moments
	mu0 <- pred[sub]
	cov0 <- matrix(pred[c(4,7,8,7,5,9,8,9,6)],3)[sub,sub]
	# Observed moments
	mu1 <- obs[sub]
	cov1 <- matrix(obs[c(4,7,8,7,5,9,8,9,6)],3)[sub,sub]
	# KL divergence
	KL.div(mu0,cov0,mu1,cov1)
}



# =====================================================================================================================
# ============================================ Tests with IBD model ===================================================
# =====================================================================================================================


# Blood -> Liver only

# ------------------------------ Reference parameter values -------------------------------
ibd.par.ref <- c(cL=0.4,cS=0,eL=0,eS=0,kL=0.6,kS=0,rL=0.4,rS=0)
ibd.M0.ref <- c(1000,0,0, 1000,0,0, 0,0,0)
ibd.init.rpois <- function() {c(rpois(1,ibd.M0.ref[1]),0,0)}
ibd.t.ref <- 12

# ------------------------------ IBD Gillespie simulations -------------------------------

# Generate 10,000 Gillespie simulations - Poisson inoculum
ibd.gillespie.sim <- WITS.gillespie.sim(ibd.t.ref,ibd.par.ref,ibd.init.rpois,10000)
ibd.gillespie.moments <- WITS.gillespie.moments(ibd.gillespie.sim)

# Predicted moments
ibd.pred.mom <- WITS.moment.sol.steps(ibd.t.ref,ibd.par.ref,ibd.M0.ref)

# Comparison of moments
round(cbind(ibd.gillespie.moments,ibd.pred.mom),3)[c(1,2,4,5,7),]

# KL divergence from pred to simulated
KL.div.p2o(ibd.pred.mom,ibd.gillespie.moments,1:2)


# --------------------- Gillespie noise -------------------------------
# Distribution of KL divergences across many samples of simulations
n.wits <- 10*c(1,5,20,50,100)
ibd.KL.sim.dist <- foreach(n = n.wits, .combine=cbind) %dopar% {
	sapply(1:1000,function(i){
		KL.div.p2o(ibd.pred.mom, WITS.gillespie.moments(ibd.gillespie.sim[sample(nrow(ibd.gillespie.sim),n,T),]),1:2)
	})
}

colnames(ibd.KL.sim.dist) <- n.wits 

boxplot(ibd.KL.sim.dist,log="y",xlab="Number of simulations (WITS)",ylab="KL divergence")


# --------------- Alternative observation times -----------------------

# Observation time points
t.obs <- c(0.5,1,2,4,8,16)

# Reference simulations for each observation time point
ibd.gillespie.sim.mt <- lapply(t.obs, function(t) WITS.gillespie.sim(t,ibd.par.ref,ibd.init.rpois,10000))
ibd.gillespie.moments.mt <- lapply(ibd.gillespie.sim.mt, WITS.gillespie.moments)

# Distribution KL divergences from correct predictions to sets of 80 simulations at each time point
n.sim <- 80
ibd.KL.sim.t.obs <- foreach(i=1:length(t.obs), .combine=cbind) %dopar% {
	ibd.pred.mom.t <- WITS.moment.sol.steps(t.obs[i],ibd.par.ref,ibd.M0.ref,h = 0.1)
	sapply(1:1000,function(x){
		KL.div.p2o(ibd.pred.mom.t, WITS.gillespie.moments(ibd.gillespie.sim.mt[[i]][sample(nrow(ibd.gillespie.sim),n.sim,T),]),1:2)
	})
}
colnames(ibd.KL.sim.t.obs) <- t.obs

boxplot(ibd.KL.sim.t.obs,log="y",xlab="Observation time",ylab="KL divergence to sample of 80 simulations")




# --------------------- Sensitivity to model parameters --------------------

# -------- Sensitivity to cL -----
cL.range <- (1:10)/10
# Based on 10 mice
n.sim <- 80

ibd.KL.sim.cL <- foreach(cL.i = cL.range, .combine=cbind) %dopar% {
	ibd.pred.mom.i <- WITS.moment.sol.steps(ibd.t.ref,replace.par(ibd.par.ref,c(cL=cL.i)),ibd.M0.ref)
	sapply(1:1000,function(x){
		KL.div.p2o(ibd.pred.mom.i, WITS.gillespie.moments(ibd.gillespie.sim[sample(nrow(ibd.gillespie.sim),n.sim,T),]),1:2)
	})
}
colnames(ibd.KL.sim.cL) <- cL.range

boxplot(ibd.KL.sim.cL,log="y",xlab="Value of cL for prediction",ylab="KL divergence to sample of 80 simulations")

# -------- Sensitivity to kL -----
kL.range <- (1:10)/10
# Based on 10 mice
n.sim <- 80

ibd.KL.sim.kL <- foreach(kL.i = kL.range, .combine=cbind) %dopar% {
	ibd.pred.mom.i <- WITS.moment.sol.steps(ibd.t.ref,replace.par(ibd.par.ref,c(kL=kL.i)),ibd.M0.ref)
	sapply(1:1000,function(x){
		KL.div.p2o(ibd.pred.mom.i, WITS.gillespie.moments(ibd.gillespie.sim[sample(nrow(ibd.gillespie.sim),n.sim,T),]),1:2)
	})
}
colnames(ibd.KL.sim.kL) <- kL.range

boxplot(ibd.KL.sim.kL,log="y",xlab="Value of kL for prediction",ylab="KL divergence to sample of 80 simulations")


# -------- Sensitivity to rL -----
rL.range <- (1:10)/10
# Based on 10 mice
n.sim <- 80

ibd.KL.sim.rL <- foreach(rL.i = rL.range, .combine=cbind) %dopar% {
	ibd.pred.mom.i <- WITS.moment.sol.steps(ibd.t.ref,replace.par(ibd.par.ref,c(rL=rL.i)),ibd.M0.ref)
	sapply(1:1000,function(x){
		KL.div.p2o(ibd.pred.mom.i, WITS.gillespie.moments(ibd.gillespie.sim[sample(nrow(ibd.gillespie.sim),n.sim,T),]),1:2)
	})
}
colnames(ibd.KL.sim.rL) <- rL.range

boxplot(ibd.KL.sim.rL,log="y",xlab="Value of rL for prediction",ylab="KL divergence to sample of 80 simulations")


# -------- Sensitivity to kL when keeping kL-rL=0.2 -----
kL.2.range <- (3:10)/10
# Based on 10 mice
n.sim <- 80

ibd.KL.sim.kL.2 <- foreach(kL.i = kL.2.range, .combine=cbind) %dopar% {
	ibd.pred.mom.i <- WITS.moment.sol.steps(ibd.t.ref,replace.par(ibd.par.ref,c(kL=kL.i,rL=kL.i-0.2)),ibd.M0.ref)
	sapply(1:1000,function(x){
		KL.div.p2o(ibd.pred.mom.i, WITS.gillespie.moments(ibd.gillespie.sim[sample(nrow(ibd.gillespie.sim),n.sim,T),]),1:2)
	})
}
colnames(ibd.KL.sim.kL.2) <- kL.2.range

boxplot(ibd.KL.sim.kL.2,log="y",xlab="Value of kL for prediction",ylab="KL divergence to sample of 80 simulations")


# -------- Sensitivity to kL when keeping kL-rL=0.2 -----
kL.3.range <- (3:10)/10
# Based on 20 mice
n.sim <- 1000

ibd.KL.sim.kL.3 <- foreach(kL.i = kL.3.range, .combine=cbind) %dopar% {
	ibd.pred.mom.i <- WITS.moment.sol.steps(ibd.t.ref,replace.par(ibd.par.ref,c(kL=kL.i,rL=kL.i-0.2)),ibd.M0.ref)
	sapply(1:1000,function(x){
		KL.div.p2o(ibd.pred.mom.i, WITS.gillespie.moments(ibd.gillespie.sim[sample(nrow(ibd.gillespie.sim),n.sim,T),]),1:2)
	})
}
colnames(ibd.KL.sim.kL.3) <- kL.3.range

boxplot(ibd.KL.sim.kL.3,log="y",xlab="Value of kL for prediction",ylab="KL divergence to sample of 1000 simulations")




# ===================================== Parameter inference by KL divergence minimisation ============================================

# Estimate the parameter values that minimise the KL divergence
ibd.KL.optim <- powell(c(cL=0.1,kL=0.1,rL=0.1), function(par){
	if(min(par)<0) return(1E100)
	names(par) <- c("cL","kL","rL")
	par.i <- replace.par(ibd.par.ref,par)
	# Calculate the moments
	mom.i <- WITS.moment.sol.steps(ibd.t.ref,par.i,ibd.M0.ref)
	# Calculate the KL divergence to the set of simulations from the predicted moments
	div <- KL.div.p2o(mom.i,ibd.gillespie.moments,1:2)
	print(c(par,div=div))
	div
},control=list(rhoend=1E-5))

c(ibd.KL.optim$par,log10.D=log10(ibd.KL.optim$value),iter=ibd.KL.optim$counts)
# cL            kL            rL       log10.D iter.function 
# 0.3996947     0.5958288     0.3958707    -4.2386044   167.0000000 


# ------------------- Bootstrapped KL estimates --------------------------

# Generate 1000 datasets of 80 WITS
ibd.KL.optim.bs <- foreach(i=1:1000, .combine=rbind) %dopar% {
	sub <- sample(1:nrow(ibd.gillespie.sim),80,replace = T)
	sub.mom <- WITS.gillespie.moments(ibd.gillespie.sim[sub,])
	sol <- powell(c(cL=0.1,kL=0.1,rL=0.1), function(par){
		if(min(par)<0) return(1E100)
		names(par) <- c("cL","kL","rL")
		par.i <- replace.par(ibd.par.ref,par)
		# Calculate the KL divergence  from the predicted moments to the subset of simulations
		KL.div.p2o(WITS.moment.sol.steps(ibd.t.ref,par.i,ibd.M0.ref),sub.mom,1:2)
	})
	c(sol$par,D=sol$value)
}

summary(ibd.KL.optim.bs)
boxplot(ibd.KL.optim.bs)
plot(ibd.KL.optim.bs[,'kL'],ibd.KL.optim.bs[,'rL'],pch=".")
lines(ellipse(cor(ibd.KL.optim.bs[,"kL"],ibd.KL.optim.bs[,"rL"]),centre=c(mean(ibd.KL.optim.bs[,"kL"]),mean(ibd.KL.optim.bs[,"rL"])),scale=c(sd(ibd.KL.optim.bs[,"kL"]),sd(ibd.KL.optim.bs[,"rL"]))),col="blue",lwd=2)



# ------------------- Inference across different parameter values ------------------------------

# Grid of parameter values (cL,kL,rL)
ibd.par.grid.KL <- expand.grid(cL=0.2*(1:3), kL=0.2*(1:3),rL=0.2*(1:3))

t.1 <- 2
t.2 <- 12

# Generate series of 1000 Gillespie simulations at t=2 and t=12h
ibd.sim.1.grid <- foreach(i=1:nrow(ibd.par.grid.KL)) %dopar% {
	WITS.gillespie.sim(t.1,replace.par(ibd.par.ref,unlist(ibd.par.grid.KL[i,])),ibd.init.rpois,1000)
}

ibd.sim.2.grid <- foreach(i=1:nrow(ibd.par.grid.KL)) %dopar% {
	WITS.gillespie.sim(t.2,replace.par(ibd.par.ref,unlist(ibd.par.grid.KL[i,])),ibd.init.rpois,1000)
}

# Plot simulations
par(mfrow=c(3,9),mar=c(2,2,1,0.5),oma=c(0,0,0,0))
for(i in 1:nrow(ibd.par.grid.KL)){
	plot(ibd.sim.1.grid[[i]][,4],ibd.sim.1.grid[[i]][,5],xlim=c(1,1000),ylim=c(1,100000),log="xy",pch=".",col="blue")
	points(ibd.sim.2.grid[[i]][,4],ibd.sim.2.grid[[i]][,5],pch=".",col="green2")
}




# Inference using t=2h only

# Estimate parameters that minimise the sum of the two KL divergences
ibd.KL.estimates.grid.1 <- foreach(k=1:nrow(ibd.par.grid.KL)) %do% {
	print(ibd.par.grid.KL[k,])
	# Estimate parameters for 128 subsets of 80 simulations
	foreach(i=1:128, .combine=rbind) %dopar% {
		sub.1 <- sample(1:nrow(ibd.sim.1.grid[[k]]),80,replace = T)
		sub.mom.1 <- WITS.gillespie.moments(ibd.sim.1.grid[[k]][sub.1,])
		op <- powell(c(cL=0.1,kL=0.1,rL=0.1), function(par){
			if(min(par)<0) return(1E100)
			names(par) <- c("cL","kL","rL")
			par.i <- replace.par(ibd.par.ref,par)
			KL.div.p2o(WITS.moment.sol.steps(t.1,par.i,ibd.M0.ref),sub.mom.1,1:2)
		})
		unlist(op$par)
	}
}

par(mfcol=c(3,9),mar=c(2,2,1,0.5),oma=c(0,0,0,0))
for(i in 1:nrow(ibd.par.grid.KL)){
	boxplot(ibd.KL.estimates.grid.1[[i]],ylim=c(0,1))
	points(1:3,ibd.par.grid.KL[i,],col="red",pch="+",cex=3)
}


# Inference using t=12h only

# Estimate parameters that minimise the sum of the two KL divergences
ibd.KL.estimates.grid.2 <- foreach(k=1:nrow(ibd.par.grid.KL)) %do% {
	print(ibd.par.grid.KL[k,])
	# Estimate parameters for 128 subsets of 80 simulations at each time point
	foreach(i=1:128, .combine=rbind) %dopar% {
		sub.2 <- sample(1:nrow(ibd.sim.2.grid[[k]]),80,replace = T)
		sub.mom.2 <- WITS.gillespie.moments(ibd.sim.2.grid[[k]][sub.2,])
		op <- powell(c(cL=0.1,kL=0.1,rL=0.1), function(par){
			if(min(par)<0) return(1E100)
			names(par) <- c("cL","kL","rL")
			par.i <- replace.par(ibd.par.ref,par)
			KL.div.p2o(WITS.moment.sol.steps(t.2,par.i,ibd.M0.ref),sub.mom.2,1:2)
		})
		unlist(op$par)
	}
}

par(mfrow=c(3,9),mar=c(2,2,1,0.5),oma=c(0,0,0,0))
for(i in 1:nrow(ibd.par.grid.KL)){
	boxplot(ibd.KL.estimates.grid.2[[i]],ylim=c(0,1))
	points(1:3,ibd.par.grid.KL[i,],col="red",pch="+",cex=3)
}


# ################################ Save objects for report ##################################

save(ibd.par.ref,ibd.M0.ref,ibd.pred.mom,ibd.gillespie.moments,ibd.KL.sim.dist,ibd.KL.sim.t.obs,ibd.KL.sim.cL,ibd.KL.sim.kL,ibd.KL.sim.rL,ibd.KL.sim.kL.2,ibd.KL.optim, ibd.KL.optim.bs,ibd.KL.estimates.grid.1,ibd.KL.estimates.grid.2,ibd.par.grid.KL,t.1,t.2, file="Moments/KL_tests_1.RData")





# =====================================================================================================================
# 								TESTS WITH B-E-S MODEL
# =====================================================================================================================


# ================================= Bottleneck-Spillover scenario =============================================

# -------------- Parameter values ------------------------
bes.M0 <- c(20,0,0, 20,0,0, 0,0,0)
bes.init.rpois <- function() {c(rpois(1,bes.M0[1]),0,0)}

bes.par.B <- c(cL=1,cS=1,eL=0,eS=0,kL=1,kS=1,rL=0.8,rS=0.8)
bes.t.B <- 6
bes.par.E <- c(cL=1,cS=1,eL=0,eS=0,kL=0.3,kS=0.3,rL=0.5,rS=0.5)
bes.t.E <- 24
bes.par.S <- c(cL=1,cS=1,eL=0.1,eS=0.1,kL=0.3,kS=0.3,rL=0.5,rS=0.5)
bes.t.S <- 48

bes.t.obs.0 <- 0.5


# ------------------ 10000 Gillespie Simulations ------------------------------------

bes.sim.0 <- WITS.gillespie.sim(bes.t.obs.0,bes.par.B,bes.init.rpois,10000)
bes.sim.B <- WITS.gillespie.sim(bes.t.B,bes.par.B,bes.init.rpois,10000)
bes.sim.E <- WITS.gillespie.sim(bes.t.E-bes.t.B, bes.par.E, function(){bes.sim.B[sample(1:nrow(bes.sim.B),1),4:6]},10000)
bes.sim.S <- WITS.gillespie.sim(bes.t.S-bes.t.E, bes.par.S, function(){bes.sim.E[sample(1:nrow(bes.sim.E),1),4:6]},10000)

# Moments
bes.sim.mom.0 <- WITS.gillespie.moments(bes.sim.0)
bes.sim.mom.B <- WITS.gillespie.moments(bes.sim.B)
bes.sim.mom.E <- WITS.gillespie.moments(bes.sim.E)
bes.sim.mom.S <- WITS.gillespie.moments(bes.sim.S)


# L-S Covariance ellipses form Gillespie simulations 
ell.bes.sim.0 <- ellipse(cor(bes.sim.0[,5],bes.sim.0[,6]),centre=c(mean(bes.sim.0[,5]),mean(bes.sim.0[,6])),scale=c(sd(bes.sim.0[,5]),sd(bes.sim.0[,6])))
ell.bes.sim.B <- ellipse(cor(bes.sim.B[,5],bes.sim.B[,6]),centre=c(mean(bes.sim.B[,5]),mean(bes.sim.B[,6])),scale=c(sd(bes.sim.B[,5]),sd(bes.sim.B[,6])))
ell.bes.sim.E <- ellipse(cor(bes.sim.E[,5],bes.sim.E[,6]),centre=c(mean(bes.sim.E[,5]),mean(bes.sim.E[,6])),scale=c(sd(bes.sim.E[,5]),sd(bes.sim.E[,6])))
ell.bes.sim.S <- ellipse(cor(bes.sim.S[,5],bes.sim.S[,6]),centre=c(mean(bes.sim.S[,5]),mean(bes.sim.S[,6])),scale=c(sd(bes.sim.S[,5]),sd(bes.sim.S[,6])))


# --------------- Plot simulations -----------------------------
par(mfrow=c(1,2),mar=c(5,5,2,1),cex.lab=1.4)

# cfu B-L
plot(bes.sim.0[,1],bes.sim.0[,3],col="black",xlim=c(0,60),ylim=c(0,60),main="Blood v Spleen",xlab="Blood", ylab="Spleen")
points(bes.sim.0[,4],bes.sim.0[,6],pch=2,col="blue")
points(bes.sim.B[,4],bes.sim.B[,6],pch=3,col="green3")
legend("topright",paste(c(0,0.5,6,24,48),'h'),pch=1:5,col=c("black","blue","green3","red","purple"))

# cfu L-S
plot(1,1,type="n",xlim=c(1,1E5),ylim=c(1,1E5),log="xy",xaxt="n",yaxt="n", xlab="Liver", ylab="Spleen", main="Liver v Spleen")
axis(1,c(1,1+10^(0:4)),c(0,10^(0:4)))
axis(2,c(1,1+10^(0:4)),c(0,10^(0:4)))
points(bes.sim.0[,5]+1,bes.sim.0[,6]+1,pch=2,col="blue")
points(bes.sim.B[,5]+1,bes.sim.B[,6]+1,pch=3,col="green3")
points(bes.sim.E[,5]+1,bes.sim.E[,6]+1,pch=4,col="red")
points(bes.sim.S[,5]+1,bes.sim.S[,6]+1,pch=5,col="purple")
abline(h=1.5,lty=3)
abline(v=1.5,lty=3)



# lines(ell.bes.sim.0+1,col="blue")
# lines(ell.bes.sim.B+1,col="green3")
# lines(ell.bes.sim.E+1,col="red")
# lines(ell.bes.sim.S+1,col="purple")



# ------------------------- Compute Moments ---------------------------------------------------------

bes.mom.0 <- WITS.moment.sol.steps(bes.t.obs.0,bes.par.B,bes.M0,h = 0.01)
bes.mom.B <- WITS.moment.sol.steps(bes.t.B,bes.par.B,bes.M0)
bes.mom.E <- WITS.moment.sol.steps(bes.t.E-bes.t.B,bes.par.E,bes.mom.B)
bes.mom.S <- WITS.moment.sol.steps(bes.t.S-bes.t.E,bes.par.S,bes.mom.E)


# Compare with simulated moments
round(cbind(bes.sim.mom.0,bes.mom.0),2)
round(cbind(bes.sim.mom.B,bes.mom.B),2)
round(cbind(bes.sim.mom.E,bes.mom.E),2)
round(cbind(bes.sim.mom.S,bes.mom.S))

# KL divergence
KL.div.p2o(bes.mom.0,bes.sim.mom.0,1:3)
KL.div.p2o(bes.mom.B,bes.sim.mom.B,2:3)
KL.div.p2o(bes.mom.E,bes.sim.mom.E,2:3)
KL.div.p2o(bes.mom.S,bes.sim.mom.S,2:3)
KL.div.p2o(bes.mom.S,bes.sim.mom.S,1:3)



# =================================== Parameter estimation by KL optimisation - Bottleneck (0-6h) ============================

bes.par.def <- c(cL=0.1,cS=0.1,eL=0.1,eS=0.1,kL=0.1,kS=0.1,rL=0.1,rS=0.1)
bes.par.ckr <- c(cL=0.1,cS=0.1,kL=0.1,kS=0.1,rL=0.1,rS=0.1)


# Bottleneck: Estimate the parameter values that minimise the KL divergence over the first 6 hours only.
# Assume eL=eS=0
# Minimise the sums of KL at 05 and 6h
bes.KL.optim.ckr.B <- powell(bes.par.ckr, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.ckr)
	par.i <- replace.par(bes.par.B,par)
	# Calculate the moments at 0.5 and 6
	mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.01)
	mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.05)
	# All compartments at t=0.5h 
	div.0 <- KL.div.p2o(mom.i.0,bes.sim.mom.0,1:3) 
	# Liver and Spleen only at t=6
	div.B <- KL.div.p2o(mom.i.B,bes.sim.mom.B,2:3)
	
	print(c(par,div.0=div.0,div.B=div.B))
	div.0+div.B
},control=list(rhoend=1E-5))


round(c(bes.KL.optim.ckr.B$par,log10.D=log10(bes.KL.optim.ckr.B$value),iter=as.numeric(bes.KL.optim.ckr.B$counts)),3)
# cL      cS      kL      kS      rL      rS log10.D    iter 
# 1.002   0.994   0.990   0.968   0.790   0.771  -3.529 674.000 

round(cbind(target=bes.par.B[names(bes.KL.optim.ckr.B$par)],estimate=bes.KL.optim.ckr.B$par),3)


# Assume eL=eS=0
# Minimise the product of KL at 05 and 6h
bes.KL.optim.ckr.B.2 <- powell(bes.par.ckr*rnorm(6,5,1), function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.ckr)
	par.i <- replace.par(bes.par.B,par)
	# Calculate the moments at 0.5 and 6
	mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.01)
	mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.05)
	# All compartments at t=0.5h 
	div.0 <- KL.div.p2o(mom.i.0,bes.sim.mom.0,1:3) 
	# Liver and Spleen only at t=6
	div.B <- KL.div.p2o(mom.i.B,bes.sim.mom.B,2:3)
	
	print(c(par,div.0=div.0,div.B=div.B))
	div.0*div.B
},control=list(rhoend=1E-5))


round(c(bes.KL.optim.ckr.B.2$par,log10.D=log10(bes.KL.optim.ckr.B.2$value),iter=as.numeric(bes.KL.optim.ckr.B.2$counts)),3)
# cL      cS      kL      kS      rL      rS log10.D    iter 
# 1.003   0.994   0.972   0.976   0.773   0.777 -10.051 779.000 

round(cbind(target=bes.par.B[names(bes.KL.optim.ckr.B$par)],est.sum=bes.KL.optim.ckr.B$par,est.prod=bes.KL.optim.ckr.B.2$par),3)



# Assume eL=eS=0
# Minimise the max of KL at 05 and 6h
bes.KL.optim.ckr.B.3 <- powell(bes.par.ckr*rnorm(6,5,1), function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.ckr)
	par.i <- replace.par(bes.par.B,par)
	# Calculate the moments at 0.5 and 6
	mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.01)
	mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.05)
	# All compartments at t=0.5h 
	div.0 <- KL.div.p2o(mom.i.0,bes.sim.mom.0,1:3) 
	# Liver and Spleen only at t=6
	div.B <- KL.div.p2o(mom.i.B,bes.sim.mom.B,2:3)
	
	print(c(par,div.0=div.0,div.B=div.B))
	max(div.0,div.B)
},control=list(rhoend=1E-5))


round(c(bes.KL.optim.ckr.B.3$par,log10.D=log10(bes.KL.optim.ckr.B.3$value),iter=as.numeric(bes.KL.optim.ckr.B.3$counts)),3)

round(cbind(target=bes.par.B[names(bes.KL.optim.ckr.B$par)],est.sum=bes.KL.optim.ckr.B$par,est.prod=bes.KL.optim.ckr.B.2$par,est.max=bes.KL.optim.ckr.B.3$par),3)
#    target est.sum est.prod est.max
# cL    1.0   1.002    1.003   1.002
# cS    1.0   0.994    0.994   0.995
# kL    1.0   0.990    0.972   0.933
# kS    1.0   0.968    0.976   0.922
# rL    0.8   0.790    0.773   0.740
# rS    0.8   0.771    0.777   0.729



# Assume eL=eS=0, ignore blood data, Minimise the sum of KL at 05 and 6h
bes.KL.optim.b.ckr.LS.sum <- powell(bes.par.ckr, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.ckr)
	par.i <- replace.par(bes.par.B,par)
	# Calculate the moments at 0.5 and 6
	mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.01)
	mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.05)
	# All compartments at t=0.5h 
	div.0 <- KL.div.p2o(mom.i.0,bes.sim.mom.0,2:3) 
	# Liver and Spleen only at t=6
	div.B <- KL.div.p2o(mom.i.B,bes.sim.mom.B,2:3)
	
	print(c(par,div.0=div.0,div.B=div.B))
	div.0 + div.B
},control=list(rhoend=1E-4))

round(c(bes.KL.optim.b.ckr.LS.sum$par,log10.D=log10(bes.KL.optim.b.ckr.LS.sum$value),iter=as.numeric(bes.KL.optim.b.ckr.LS.sum$counts)),3)
# cL      cS      kL      kS      rL      rS log10.D    iter 
# 1.006   0.998   0.989   0.967   0.789   0.770  -3.680 658.000 


# Assume eL=eS=0, ignore blood data, minimise the max of KL at 05 and 6h
bes.KL.optim.b.ckr.LS.max <- powell(bes.par.ckr, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.ckr)
	par.i <- replace.par(bes.par.B,par)
	# Calculate the moments at 0.5 and 6
	mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.01)
	mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.05)
	# All compartments at t=0.5h 
	div.0 <- KL.div.p2o(mom.i.0,bes.sim.mom.0,2:3) 
	# Liver and Spleen only at t=6
	div.B <- KL.div.p2o(mom.i.B,bes.sim.mom.B,2:3)
	
	print(c(par,div.0=div.0,div.B=div.B))
	max(div.0,div.B)
},control=list(rhoend=1E-4))

round(c(bes.KL.optim.b.ckr.LS.max$par,log10.D=log10(bes.KL.optim.b.ckr.LS.max$value),iter=as.numeric(bes.KL.optim.b.ckr.LS.max$counts)),3)
# cL       cS       kL       kS       rL       rS  log10.D     iter 
# 1.086    1.103    1.014    1.062    0.798    0.850   -2.502 4335.000 

# Assume eL=eS=0, ignore blood data, minimise the product of KL at 05 and 6h
bes.KL.optim.b.ckr.LS.prod <- powell(bes.par.ckr*2, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.ckr)
	par.i <- replace.par(bes.par.B,par)
	# Calculate the moments at 0.5 and 6
	mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.01)
	mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.05)
	# All compartments at t=0.5h 
	div.0 <- KL.div.p2o(mom.i.0,bes.sim.mom.0,2:3) 
	# Liver and Spleen only at t=6
	div.B <- KL.div.p2o(mom.i.B,bes.sim.mom.B,2:3)
	
	print(c(par,div.0=div.0,div.B=div.B))
	div.0*div.B
},control=list(rhoend=1E-5))

round(c(bes.KL.optim.b.ckr.LS.prod$par,log10.D=log10(bes.KL.optim.b.ckr.LS.prod$value),iter=as.numeric(bes.KL.optim.b.ckr.LS.prod$counts)),3)
# cL      cS      kL      kS      rL      rS log10.D    iter 
# 1.089   1.106   0.885   0.937   0.680   0.735  -6.388 807.000 




# Estimate all 8 parameters, assuming known inoculum distribution
# Minimise the sums of KL at 05 and 6h
xx <-0
bes.KL.optim.b.cekr.BLS.sum <- powell(bes.par.def, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.def)
	par.i <- replace.par(bes.par.B,par)
	# Calculate the moments at 0.5 and 6
	mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.01)
	mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.05)
	# All compartments at t=0.5h 
	div.0 <- KL.div.p2o(mom.i.0,bes.sim.mom.0,1:3) 
	# Liver and Spleen only at t=6
	div.B <- KL.div.p2o(mom.i.B,bes.sim.mom.B,2:3)
	xx <<- xx+1
	if(xx%%10==1) print(round(c(par,div.0=div.0,div.B=div.B),3))
	if(!is.finite(div.0+div.B)) return(1E100)
	div.0+div.B
},control=list(rhoend=1E-5))

round(bes.KL.optim.b.cekr.BLS.sum$par,5)



# ------------------------------------ Inter-sample variability ----------------------------------

# Generate 64 datasets of 100 WITS and estimate c,k,r parameters from each
bes.KL.optim.B.bs <- foreach(i=1:64, .combine=rbind) %dopar% {
	sub.mom.0 <- WITS.gillespie.moments(bes.sim.0[sample(1:nrow(bes.sim.0),100,replace = T),])
	sub.mom.B <- WITS.gillespie.moments(bes.sim.B[sample(1:nrow(bes.sim.B),100,replace = T),])
	sol <- powell(bes.par.ckr*rnorm(6,5,1), function(par){
		if(min(par)<0) return(1E100)
		names(par) <- names(bes.par.ckr)
		par.i <- replace.par(bes.par.B,par)
		# Calculate the moments at 0.5 and 6
		mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.02)
		mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.1)
		# All compartments at t=0.5h 
		div.0 <- KL.div.p2o(mom.i.0,sub.mom.0,1:3) 
		# Liver and Spleen only at t=6
		div.B <- KL.div.p2o(mom.i.B,sub.mom.B,2:3)
		if(!is.finite(div.0+div.B)) return(1E100)
		div.0+div.B
	},control=list(rhoend=1E-4))
	c(sol$par,D=sol$value)
}

summary(bes.KL.optim.B.bs)
boxplot(bes.KL.optim.B.bs)
points(1:6,bes.par.B[names(bes.KL.optim.ckr.B$par)],pch="+",col="red",cex=2)

plot(ibd.KL.optim.bs[,'kL'],ibd.KL.optim.bs[,'rL'],pch=".")
lines(ellipse(cor(ibd.KL.optim.bs[,"kL"],ibd.KL.optim.bs[,"rL"]),centre=c(mean(ibd.KL.optim.bs[,"kL"]),mean(ibd.KL.optim.bs[,"rL"])),scale=c(sd(ibd.KL.optim.bs[,"kL"]),sd(ibd.KL.optim.bs[,"rL"]))),col="blue",lwd=2)





# =================================== Parameter estimation by KL optimisation - Expansion (6-24h) ============================

bes.par.kr <- c(kL=0.1,kS=0.1,rL=0.1,rS=0.1)


# Assume eL=eS=0
# Only estimate k and r
# Calculate the moments at 24h for L and S only, with initial conditions at t=6h given by data
bes.KL.optim.kr.E <- powell(bes.par.kr, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.kr)
	par.i <- replace.par(bes.par.E,par)
	mom.i.E <- WITS.moment.sol.steps(bes.t.E-bes.t.B,par.i,bes.sim.mom.B,h=0.2)
	div.E <- KL.div.p2o(mom.i.E,bes.sim.mom.E,2:3)
#	print(c(par,div=div.E))
	if(is.finite(div.E)) div.E else 1E100
},control=list(rhoend=1E-5))

cbind(bes.par.E[names(bes.KL.optim.kr.E$par)],bes.KL.optim.kr.E$par)


# ------------------------------------ Inter-sample variability ----------------------------------

# Generate 64 datasets of 100 WITS and estimate k,r parameters from each
bes.KL.optim.E.bs <- foreach(i=1:64, .combine=rbind) %dopar% {
	sub.mom.B <- WITS.gillespie.moments(bes.sim.B[sample(1:nrow(bes.sim.B),100,replace = T),])
	sub.mom.E <- WITS.gillespie.moments(bes.sim.E[sample(1:nrow(bes.sim.E),100,replace = T),])
	sol <- powell(bes.par.kr, function(par){
		if(min(par)<0) return(1E100)
		names(par) <- names(bes.par.kr)
		par.i <- replace.par(bes.par.E,par)
		mom.i.E <- WITS.moment.sol.steps(bes.t.E-bes.t.B,par.i,sub.mom.B,h=0.05)
		div.E <- KL.div.p2o(mom.i.E,sub.mom.E,2:3)
		if(is.finite(div.E)) div.E else 1E100
	},control=list(rhoend=1E-5))
	c(sol$par,D=sol$value)
}

summary(bes.KL.optim.E.bs)
boxplot(bes.KL.optim.E.bs,ylim=c(0,1.2))
points(1:4,bes.par.E[names(bes.par.kr)],pch="+",col="red",cex=2)



# =================================== Parameter estimation by KL optimisation - Spillover (24-48h) ============================

# Calculate the moments at 48h for all 3 organs, with initial conditions at t=24h given by data
bes.KL.optim.S <- powell(bes.par.def*5, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.def)
	par.i <- replace.par(bes.par.S,par)
	mom.i.S <- WITS.moment.sol.steps(bes.t.S-bes.t.E,par.i,bes.sim.mom.E,h=0.2)
	div.S <- KL.div.p2o(mom.i.S,bes.sim.mom.S,1:3)
	xx <<- xx+1
	if(xx%%20==1) print(c(par,div=div.S),digits=4) 
	if(is.finite(div.S)) div.S else 1E100
},control=list(rhobeg=0.5,rhoend=1E-5))

round(cbind(bes.par.S,bes.KL.optim.S$par),3)

bes.KL.optim.S.2 <- powell(bes.KL.optim.S$par, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.def)
	par.i <- replace.par(bes.par.S,par)
	mom.i.S <- WITS.moment.sol.steps(bes.t.S-bes.t.E,par.i,bes.sim.mom.E,h=0.2)
	div.S <- KL.div.p2o(mom.i.S,bes.sim.mom.S,1:3)
	xx <<- xx+1
	if(xx%%20==1) print(c(par,div=div.S),digits=4) 
	if(is.finite(div.S)) div.S else 1E100
},control=list(rhoend=1E-7))

round(cbind(bes.par.S,bes.KL.optim.S$par,bes.KL.optim.S.2$par),3)

# ------------------------------------ Inter-sample variability ----------------------------------

# Generate 32 datasets of 100 WITS and estimate all parameters from each
bes.KL.optim.S.bs <- foreach(i=1:32, .combine=rbind) %dopar% {
	sub.mom.E <- WITS.gillespie.moments(bes.sim.E[sample(1:nrow(bes.sim.E),100,replace = T),])
	sub.mom.S <- WITS.gillespie.moments(bes.sim.S[sample(1:nrow(bes.sim.S),100,replace = T),])
	sol <- powell(bes.par.def*5, function(par){
		if(min(par)<0) return(1E100)
		names(par) <- names(bes.par.def)
		par.i <- replace.par(bes.par.S,par)
		mom.i.S <- WITS.moment.sol.steps(bes.t.S-bes.t.E,par.i,sub.mom.E,h=0.2)
		div.S <- KL.div.p2o(mom.i.S,sub.mom.S,1:3)
		if(is.finite(div.S)) div.S else 1E100
	},control=list(rhoend=1E-5))
	c(sol$par,D=sol$value)
}

boxplot(bes.KL.optim.S.bs,ylim=c(0,1.5))
points(1:8,bes.par.S[names(bes.par.def)],pch="+",col="red",cex=2)




#########################################################################################################################
# Save objects for report

save(list=ls()[grep("bes.",ls())], file="Moments/KL_tests_BES.RData")







