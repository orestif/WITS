#############################
#   SALMONELLA WITS MODEL
#     MOMENTS EQUATION
#          TESTS
#############################

# PARAMETER INFERENCE USING the Hellinger Distance

source("Moments/WITS_moments.R")

library(dplyr)
library(ggplot2)
library(tidyr)
library(powell)
library(ellipse)



# ========================================= Hellinger distance =============================================================

#library(distrEx)

# Calculate the squared Hellinger distance between two [multinormal] distributions with mean vectors mu0 and mu1, and cov matrices cov0 and cov1
Hell.dist <- function(mu0,cov0,mu1,cov1){
	P <- (cov0+cov1)/2
	m <- mu0-mu1
	B <- as.numeric((m %*% solve(P) %*% m)/8 + log(det(P)/sqrt(det(cov0)*det(cov1)))/2)
	1-exp(-B)
}

# Calculate the Hell distance from an observed sample to a predicted distribution, both characterised by their vectors of 9 moments
# - obs: observed moments
# - pred: predicted moments
# - sub: vector specifying which organs to include, e.g. c(1,2) for B and L. Default 1:3
Hell.dist.o2p <- function(obs,pred,sub=1:3){
	# Observed moments
	mu0 <- obs[sub]
	cov0 <- matrix(obs[c(4,7,8,7,5,9,8,9,6)],3)[sub,sub]
	# Predicted moments
	mu1 <- pred[sub]
	cov1 <- matrix(pred[c(4,7,8,7,5,9,8,9,6)],3)[sub,sub]
	# Hell distance
	Hell.dist(mu0,cov0,mu1,cov1)
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

# Hell distance from pred to simulated
Hell.dist.o2p(ibd.gillespie.moments,ibd.pred.mom,1:2)


# --------------------- Gillespie noise -------------------------------
# Distribution of Hell distances across many samples of simulations
n.wits <- 8*c(1,5,10,20)
ibd.Hell.sim.dist <- foreach(n = n.wits, .combine=cbind) %dopar% {
	sapply(1:1000,function(i){
		Hell.dist.o2p( WITS.gillespie.moments(ibd.gillespie.sim[sample(nrow(ibd.gillespie.sim),n,T),]),ibd.pred.mom,1:2)
	})
}

colnames(ibd.Hell.sim.dist) <- n.wits 

boxplot(ibd.Hell.sim.dist,log="y",xlab="Number of simulations (WITS)",ylab="Hell distance")


# --------------- Alternative observation times -----------------------

# Observation time points
t.obs <- c(seq(2,12,2),c(16,20))

# Reference simulations for each observation time point
ibd.gillespie.sim.mt <- lapply(t.obs, function(t) WITS.gillespie.sim(t,ibd.par.ref,ibd.init.rpois,10000))
ibd.gillespie.moments.mt <- lapply(ibd.gillespie.sim.mt, WITS.gillespie.moments)

# Distribution Hell distances from sets of 80 simulations to correct predictions at each time point
n.sim <- 80
ibd.Hell.sim.t.obs <- foreach(i=1:length(t.obs), .combine=cbind) %dopar% {
	ibd.pred.mom.t <- WITS.moment.sol.steps(t.obs[i],ibd.par.ref,ibd.M0.ref)
	sapply(1:1000,function(x){
		Hell.dist.o2p( WITS.gillespie.moments(ibd.gillespie.sim.mt[[i]][sample(nrow(ibd.gillespie.sim),n.sim,T),]),ibd.pred.mom.t,1:2)
	})
}
colnames(ibd.Hell.sim.t.obs) <- t.obs

boxplot(ibd.Hell.sim.t.obs,log="y",xlab="Observation time",ylab="Hell distance to sample of 80 simulations")




# --------------------- Sensitivity to model parameters --------------------

# -------- Sensitivity to cL -----
cL.range <- (1:10)/10
# Based on 10 mice
n.sim <- 80

ibd.Hell.sim.cL <- foreach(cL.i = cL.range, .combine=cbind) %dopar% {
	ibd.pred.mom.i <- WITS.moment.sol.steps(ibd.t.ref,replace.par(ibd.par.ref,c(cL=cL.i)),ibd.M0.ref)
	sapply(1:1000,function(x){
		Hell.dist.o2p(WITS.gillespie.moments(ibd.gillespie.sim[sample(nrow(ibd.gillespie.sim),n.sim,T),]),ibd.pred.mom.i, 1:2)
	})
}
colnames(ibd.Hell.sim.cL) <- cL.range

boxplot(ibd.Hell.sim.cL,log="y",xlab="Value of cL for prediction",ylab="Hell distance to sample of 80 simulations")

# -------- Sensitivity to kL -----
kL.range <- (1:10)/10
# Based on 10 mice
n.sim <- 80

ibd.Hell.sim.kL <- foreach(kL.i = kL.range, .combine=cbind) %dopar% {
	ibd.pred.mom.i <- WITS.moment.sol.steps(ibd.t.ref,replace.par(ibd.par.ref,c(kL=kL.i)),ibd.M0.ref)
	sapply(1:1000,function(x){
		Hell.dist.o2p( WITS.gillespie.moments(ibd.gillespie.sim[sample(nrow(ibd.gillespie.sim),n.sim,T),]),ibd.pred.mom.i,1:2)
	})
}
colnames(ibd.Hell.sim.kL) <- kL.range

boxplot(ibd.Hell.sim.kL,log="y",xlab="Value of kL for prediction",ylab="Hell distance to sample of 80 simulations")


# -------- Sensitivity to rL -----
rL.range <- (1:10)/10
# Based on 10 mice
n.sim <- 80

ibd.Hell.sim.rL <- foreach(rL.i = rL.range, .combine=cbind) %dopar% {
	ibd.pred.mom.i <- WITS.moment.sol.steps(ibd.t.ref,replace.par(ibd.par.ref,c(rL=rL.i)),ibd.M0.ref)
	sapply(1:1000,function(x){
		Hell.dist.o2p(WITS.gillespie.moments(ibd.gillespie.sim[sample(nrow(ibd.gillespie.sim),n.sim,T),]),ibd.pred.mom.i, 1:2)
	})
}
colnames(ibd.Hell.sim.rL) <- rL.range

boxplot(ibd.Hell.sim.rL,log="y",xlab="Value of rL for prediction",ylab="Hell distance to sample of 80 simulations")


# -------- Sensitivity to kL when keeping kL-rL=0.2 -----
kL.2.range <- (3:10)/10
# Based on 10 mice
n.sim <- 80

ibd.Hell.sim.kL.2 <- foreach(kL.i = kL.2.range, .combine=cbind) %dopar% {
	ibd.pred.mom.i <- WITS.moment.sol.steps(ibd.t.ref,replace.par(ibd.par.ref,c(kL=kL.i,rL=kL.i-0.2)),ibd.M0.ref)
	sapply(1:1000,function(x){
		Hell.dist.o2p(WITS.gillespie.moments(ibd.gillespie.sim[sample(nrow(ibd.gillespie.sim),n.sim,T),]),ibd.pred.mom.i, 1:2)
	})
}
colnames(ibd.Hell.sim.kL.2) <- kL.2.range

boxplot(ibd.Hell.sim.kL.2,log="y",xlab="Value of kL for prediction",ylab="Hell distance to sample of 80 simulations")





# ===================================== Parameter inference by Hell distance minimisation ============================================

# Estimate the parameter values that minimise the Hell distance
ibd.Hell.optim <- powell(c(cL=0.1,kL=0.1,rL=0.1), function(par){
	if(min(par)<0) return(1E100)
	names(par) <- c("cL","kL","rL")
	par.i <- replace.par(ibd.par.ref,par)
	# Calculate the moments
	mom.i <- WITS.moment.sol.steps(ibd.t.ref,par.i,ibd.M0.ref)
	# Calculate the Hell distance to the set of simulations from the predicted moments
	div <- Hell.dist.o2p(ibd.gillespie.moments,mom.i,1:2)
	print(c(par,div=div))
	div
},control=list(rhoend=1E-4))


c(ibd.Hell.optim$par,log10.D=log10(ibd.Hell.optim$value),iter=ibd.Hell.optim$counts)
# cL            kL            rL       log10.D iter.function 
# 0.4002704     0.6045877     0.4045978    -5.4424558   157.0000000 


# Values from Hell.o2p
# cL            kL            rL       log10.D iter.function 
# 0.3998718     0.5892282     0.3892500    -4.9647098   129.0000000 

# For comparison, values from Hell.p2o:
# cL            kL            rL       log10.D iter.function 
# 0.3998044     0.5989939     0.3990370    -4.0277796   148.0000000 


# ------------------- Bootstrapped Hell estimates --------------------------

# Generate 1000 datasets of 80 WITS
ibd.Hell.optim.bs <- foreach(i=1:1000, .combine=rbind) %dopar% {
	sub <- sample(1:nrow(ibd.gillespie.sim),80,replace = T)
	sub.mom <- WITS.gillespie.moments(ibd.gillespie.sim[sub,])
	sol <- powell(c(cL=0.1,kL=0.1,rL=0.1), function(par){
		if(min(par)<0) return(1E100)
		names(par) <- c("cL","kL","rL")
		par.i <- replace.par(ibd.par.ref,par)
		# Calculate the Hell distance  from the predicted moments to the subset of simulations
		Hell.dist.o2p(sub.mom,WITS.moment.sol.steps(ibd.t.ref,par.i,ibd.M0.ref),1:2)
	})
	c(sol$par,D=sol$value)
}

summary(ibd.Hell.optim.bs)
boxplot(ibd.Hell.optim.bs)
points(1:3, ibd.par.ref[c("cL","kL","rL")], col="red", pch="+", cex=3)

plot(ibd.Hell.optim.bs[,'kL'],ibd.Hell.optim.bs[,'rL'],pch=".")
lines(ellipse(cor(ibd.Hell.optim.bs[,"kL"],ibd.Hell.optim.bs[,"rL"]),centre=c(mean(ibd.Hell.optim.bs[,"kL"]),mean(ibd.Hell.optim.bs[,"rL"])),scale=c(sd(ibd.Hell.optim.bs[,"kL"]),sd(ibd.Hell.optim.bs[,"rL"]))),col="blue",lwd=2)



# ------------------- Inference across different parameter values ------------------------------

# Grid of parameter values (cL,kL,rL)
ibd.par.grid.Hell <- expand.grid(cL=0.2*(1:3), kL=0.2*(1:3),rL=0.2*(1:3))

t.1 <- 2
t.2 <- 12

# Generate series of 1000 Gillespie simulations at t=2 and t=12h
ibd.sim.1.grid <- foreach(i=1:nrow(ibd.par.grid.Hell)) %dopar% {
	WITS.gillespie.sim(t.1,replace.par(ibd.par.ref,unlist(ibd.par.grid.Hell[i,])),ibd.init.rpois,1000)
}

ibd.sim.2.grid <- foreach(i=1:nrow(ibd.par.grid.Hell)) %dopar% {
	WITS.gillespie.sim(t.2,replace.par(ibd.par.ref,unlist(ibd.par.grid.Hell[i,])),ibd.init.rpois,1000)
}

# Plot simulations
par(mfrow=c(3,9),mar=c(2,2,1,0.5),oma=c(0,0,0,0))
for(i in 1:nrow(ibd.par.grid.Hell)){
	plot(ibd.sim.1.grid[[i]][,4],ibd.sim.1.grid[[i]][,5],xlim=c(1,1000),ylim=c(1,100000),log="xy",pch=".",col="blue")
	points(ibd.sim.2.grid[[i]][,4],ibd.sim.2.grid[[i]][,5],pch=".",col="green2")
}




# Inference using t=2h only

# Estimate parameters that minimise the sum of the two Hell distances
ibd.Hell.estimates.grid.1 <- foreach(k=1:nrow(ibd.par.grid.Hell)) %do% {
	print(ibd.par.grid.Hell[k,])
	# Estimate parameters for 64 subsets of 80 simulations
	foreach(i=1:64, .combine=rbind) %dopar% {
		sub.1 <- sample(1:nrow(ibd.sim.1.grid[[k]]),80,replace = T)
		sub.mom.1 <- WITS.gillespie.moments(ibd.sim.1.grid[[k]][sub.1,])
		op <- powell(c(cL=0.1,kL=0.1,rL=0.1), function(par){
			if(min(par)<0) return(1E100)
			names(par) <- c("cL","kL","rL")
			par.i <- replace.par(ibd.par.ref,par)
			Hell.dist.o2p(sub.mom.1,WITS.moment.sol.steps(t.1,par.i,ibd.M0.ref),1:2)
		})
		unlist(op$par)
	}
}

par(mfrow=c(3,9),mar=c(2,2,1,0.5),oma=c(0,0,0,0))
for(i in 1:nrow(ibd.par.grid.Hell)){
	boxplot(ibd.Hell.estimates.grid.1[[i]],ylim=c(0,1))
	points(1:3,ibd.par.grid.Hell[i,],col="red",pch="+",cex=3)
}


# Inference using t=12h only

# Estimate parameters that minimise the sum of the two Hell distances
ibd.Hell.estimates.grid.2 <- foreach(k=1:nrow(ibd.par.grid.Hell)) %do% {
	print(ibd.par.grid.Hell[k,])
	# Estimate parameters for 128 subsets of 80 simulations at each time point
	foreach(i=1:128, .combine=rbind) %dopar% {
		sub.2 <- sample(1:nrow(ibd.sim.2.grid[[k]]),80,replace = T)
		sub.mom.2 <- WITS.gillespie.moments(ibd.sim.2.grid[[k]][sub.2,])
		op <- powell(c(cL=0.1,kL=0.1,rL=0.1), function(par){
			if(min(par)<0) return(1E100)
			names(par) <- c("cL","kL","rL")
			par.i <- replace.par(ibd.par.ref,par)
			Hell.dist.o2p(sub.mom.2,WITS.moment.sol.steps(t.2,par.i,ibd.M0.ref),1:2)
		})
		unlist(op$par)
	}
}

par(mfrow=c(3,9),mar=c(2,2,1,0.5),oma=c(0,0,0,0))
for(i in 1:nrow(ibd.par.grid.Hell)){
	boxplot(ibd.Hell.estimates.grid.2[[i]],ylim=c(0,1))
	points(1:3,ibd.par.grid.Hell[i,],col="red",pch="+",cex=3)
}


# ################################ Save objects for report ##################################

save(list=ls()[grep("ibd.",ls())], file="Moments/Hell_tests_2.RData")



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


# ------------------ 1000 Gillespie Simulations ------------------------------------

bes.sim.0 <- WITS.gillespie.sim(bes.t.obs.0,bes.par.B,bes.init.rpois,1000)
bes.sim.B <- WITS.gillespie.sim(bes.t.B,bes.par.B,bes.init.rpois,1000)
bes.sim.E <- WITS.gillespie.sim(bes.t.E-bes.t.B, bes.par.E, function(){bes.sim.B[sample(1:nrow(bes.sim.B),1),4:6]},1000)
bes.sim.S <- WITS.gillespie.sim(bes.t.S-bes.t.E, bes.par.S, function(){bes.sim.E[sample(1:nrow(bes.sim.E),1),4:6]},1000)

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




# =================================== Parameter estimation by Hell optimisation - Bottleneck (0-6h) ============================

bes.par.def <- c(cL=0.1,cS=0.1,eL=0.1,eS=0.1,kL=0.1,kS=0.1,rL=0.1,rS=0.1)
bes.par.ckr <- c(cL=0.1,cS=0.1,kL=0.1,kS=0.1,rL=0.1,rS=0.1)


# Bottleneck: Estimate the parameter values that minimise the Hell distance over the first 6 hours only.
# Assume eL=eS=0
# Minimise the sums of Hell at 05 and 6h
bes.Hell.optim.ckr.B <- powell(bes.par.ckr, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.ckr)
	par.i <- replace.par(bes.par.B,par)
	# Calculate the moments at 0.5 and 6
	mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.01)
	mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.05)
	# All compartments at t=0.5h 
	dist.0 <- Hell.dist.o2p(bes.sim.mom.0,mom.i.0,1:3) 
	# Liver and Spleen only at t=6
	dist.B <- Hell.dist.o2p(bes.sim.mom.B,mom.i.B,2:3)
	
	print(c(par,dist.0=dist.0,dist.B=dist.B))
	dist.0+dist.B
},control=list(rhoend=1E-4))


round(c(bes.Hell.optim.ckr.B$par,log10.D=log10(bes.Hell.optim.ckr.B$value),iter=as.numeric(bes.Hell.optim.ckr.B$counts)),3)

round(cbind(target=bes.par.B[names(bes.Hell.optim.ckr.B$par)],estimate=bes.Hell.optim.ckr.B$par),3)


# Assume eL=eS=0
# Minimise the product of Hell at 05 and 6h
bes.Hell.optim.ckr.B.2 <- powell(bes.par.ckr, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.ckr)
	par.i <- replace.par(bes.par.B,par)
	# Calculate the moments at 0.5 and 6
	mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.01)
	mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.05)
	# All compartments at t=0.5h 
	dist.0 <- Hell.dist.o2p(bes.sim.mom.0,mom.i.0,1:3) 
	# Liver and Spleen only at t=6
	dist.B <- Hell.dist.o2p(bes.sim.mom.B,mom.i.B,2:3)
	
	print(round(c(par,dist.0=dist.0,dist.B=dist.B),3))
	dist.0*dist.B
},control=list(rhoend=1E-4))


round(c(bes.Hell.optim.ckr.B.2$par,log10.D=log10(bes.Hell.optim.ckr.B.2$value),iter=as.numeric(bes.Hell.optim.ckr.B.2$counts)),3)
# cL      cS      kL      kS      rL      rS log10.D    iter 
# 1.052   1.050   1.028   1.008   0.818   0.807  -8.361 980.000 


# KL o2p
# cL      cS      kL      kS      rL      rS log10.D    iter 
# 1.078   1.088   1.001   1.045   0.808   0.857  -5.202 720.000 

# Estimates from KL p2o
# cL       cS       kL       kS       rL       rS  log10.D     iter 
# 1.075    1.092    0.886    0.938    0.681    0.736   -6.244 1203.000 

round(cbind(target=bes.par.B[names(bes.Hell.optim.ckr.B$par)],est.sum=bes.Hell.optim.ckr.B$par,est.prod=bes.Hell.optim.ckr.B.2$par),3)



# Assume eL=eS=0
# Minimise the max of Hell at 05 and 6h
bes.Hell.optim.ckr.B.3 <- powell(bes.par.ckr, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.ckr)
	par.i <- replace.par(bes.par.B,par)
	# Calculate the moments at 0.5 and 6
	mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.01)
	mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.05)
	# All compartments at t=0.5h 
	dist.0 <- Hell.dist.o2p(bes.sim.mom.0,mom.i.0,1:3) 
	# Liver and Spleen only at t=6
	dist.B <- Hell.dist.o2p(bes.sim.mom.B,mom.i.B,2:3)
	
	print(round(c(par,dist.0=dist.0,dist.B=dist.B),3))
	max(dist.0,dist.B)
},control=list(rhoend=1E-4))


round(c(bes.Hell.optim.ckr.B.3$par,log10.D=log10(bes.Hell.optim.ckr.B.3$value),iter=as.numeric(bes.Hell.optim.ckr.B.3$counts)),3)
# cL      cS      kL      kS      rL      rS log10.D    iter 
# 1.052   1.049   0.884   0.931   0.675   0.729  -3.597 776.000 


# KL o2p
# cL      cS      kL      kS      rL      rS log10.D    iter 
# 1.080   1.088   0.929   0.899   0.740   0.721  -2.372 847.000 


# Estimates from KL p2o
# cL       cS       kL       kS       rL       rS  log10.D     iter 
# 1.073    1.092    0.765    0.795    0.587    0.615   -1.876 2308.000 

round(cbind(target=bes.par.B[names(bes.Hell.optim.ckr.B$par)],est.sum=bes.Hell.optim.ckr.B$par,est.prod=bes.Hell.optim.ckr.B.2$par,est.max=bes.Hell.optim.ckr.B.3$par),3)



# Assume eL=eS=0, ignore blood data, Minimise the sum of Hell at 05 and 6h
bes.Hell.optim.b.ckr.LS.sum <- powell(bes.par.ckr, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.ckr)
	par.i <- replace.par(bes.par.B,par)
	# Calculate the moments at 0.5 and 6
	mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.01)
	mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.05)
	# All compartments at t=0.5h 
	dist.0 <- Hell.dist.o2p(bes.sim.mom.0,mom.i.0,2:3) 
	# Liver and Spleen only at t=6
	dist.B <- Hell.dist.o2p(bes.sim.mom.B,mom.i.B,2:3)
	
	print(round(c(par,dist.0=dist.0,dist.B=dist.B),3))
	dist.0 + dist.B
},control=list(rhoend=1E-4))

round(c(bes.Hell.optim.b.ckr.LS.sum$par,log10.D=log10(bes.Hell.optim.b.ckr.LS.sum$value),iter=as.numeric(bes.Hell.optim.b.ckr.LS.sum$counts)),3)
# cL      cS      kL      kS      rL      rS log10.D    iter 
# 1.034   1.031   1.001   1.015   0.791   0.812  -3.963 988.000 


# KL o2p
# cL      cS      kL      kS      rL      rS log10.D    iter 
# 1.033   1.042   0.982   1.016   0.791   0.829  -2.659 753.000 


# Estimates from KL p2o
# cL      cS      kL      kS      rL      rS log10.D    iter 
# 1.088   1.105   0.945   1.004   0.735   0.797  -2.266 940.000 





# Estimate all 8 parameters, assuming known inoculum distribution
# Minimise the sums of Hell at 05 and 6h
xx <-0
bes.Hell.optim.b.cekr.BLS.sum <- powell(bes.par.def*rnorm(8,5,0.5), function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.def)
	par.i <- replace.par(bes.par.B,par)
	# Calculate the moments at 0.5 and 6
	mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.01)
	mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.05)
	# All compartments at t=0.5h 
	dist.0 <- Hell.dist.o2p(bes.sim.mom.0,mom.i.0,1:3) 
	# Liver and Spleen only at t=6
	dist.B <- Hell.dist.o2p(bes.sim.mom.B,mom.i.B,2:3)
	xx <<- xx+1
	if(xx%%10==1) print(round(c(par,dist.0=dist.0,dist.B=dist.B),3))
	if(!is.finite(dist.0+dist.B)) return(1E100)
	dist.0+dist.B
},control=list(rhoend=1E-6))

round(bes.Hell.optim.b.cekr.BLS.sum$par,5)
# cL      cS      eL      eS      kL      kS      rL      rS 
# 1.08415 1.08423 0.17420 0.00000 0.63851 1.12478 0.52439 0.84326 


# KL p2o
# cL      cS      eL      eS      kL      kS      rL      rS 
# 1.07449 1.09206 0.00000 0.00000 0.49937 0.57298 0.35741 0.42879 


# Estimate all 8 parameters, assuming known inoculum distribution
# Minimise the product of Hell at 05 and 6h
bes.Hell.optim.b.cekr.BLS.prod <- powell(bes.par.def, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.def)
	par.i <- replace.par(bes.par.B,par)
	# Calculate the moments at 0.5 and 6
	mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.01)
	mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.05)
	# All compartments at t=0.5h 
	dist.0 <- Hell.dist.o2p(bes.sim.mom.0,mom.i.0,1:3) 
	# Liver and Spleen only at t=6
	dist.B <- Hell.dist.o2p(bes.sim.mom.B,mom.i.B,2:3)
	xx <<- xx+1
	if(xx%%10==1) print(round(c(par,dist.0=dist.0,dist.B=dist.B),4))
	if(!is.finite(dist.0*dist.B)) return(1E100)
	dist.0*dist.B
},control=list(rhobeg=0.5,rhoend=1E-5))

round(cbind(target=bes.par.B,est=bes.Hell.optim.b.cekr.BLS.prod$par),4)
#    target    est
# cL    1.0 1.2906
# cS    1.0 1.2190
# eL    0.0 0.2473
# eS    0.0 0.7416
# kL    1.0 1.2425
# kS    1.0 0.1565
# rL    0.8 0.6670
# rS    0.8 0.2888


# ------------------------------------ Inter-sample variability ----------------------------------

# Generate 64 datasets of 100 WITS and estimate c,k,r parameters from each
bes.Hell.optim.B.bs <- foreach(i=1:64, .combine=rbind) %dopar% {
	sub.mom.0 <- WITS.gillespie.moments(bes.sim.0[sample(1:nrow(bes.sim.0),100,replace = T),])
	sub.mom.B <- WITS.gillespie.moments(bes.sim.B[sample(1:nrow(bes.sim.B),100,replace = T),])
	sol <- powell(bes.par.ckr, function(par){
		if(min(par)<0) return(1E100)
		names(par) <- names(bes.par.ckr)
		par.i <- replace.par(bes.par.B,par)
		# Calculate the moments at 0.5 and 6
		mom.i.0 <- WITS.moment.sol.steps(bes.t.obs.0,par.i,bes.M0,h=0.02)
		mom.i.B <- WITS.moment.sol.steps(bes.t.B,par.i,bes.M0,h=0.1)
		# All compartments at t=0.5h 
		dist.0 <- Hell.dist.o2p(sub.mom.0,mom.i.0,1:3) 
		# Liver and Spleen only at t=6
		dist.B <- Hell.dist.o2p(sub.mom.B,mom.i.B,2:3)
		if(!is.finite(dist.0+dist.B)) return(1E100)
		dist.0+dist.B
	},control=list(rhoend=1E-5))
	c(sol$par,D=sol$value)
}

summary(bes.Hell.optim.B.bs)
boxplot(bes.Hell.optim.B.bs)
points(1:6,bes.par.B[names(bes.Hell.optim.ckr.B$par)],pch="+",col="red",cex=2)

plot(ibd.Hell.optim.bs[,'kL'],ibd.Hell.optim.bs[,'rL'],pch=".")
lines(ellipse(cor(ibd.Hell.optim.bs[,"kL"],ibd.Hell.optim.bs[,"rL"]),centre=c(mean(ibd.Hell.optim.bs[,"kL"]),mean(ibd.Hell.optim.bs[,"rL"])),scale=c(sd(ibd.Hell.optim.bs[,"kL"]),sd(ibd.Hell.optim.bs[,"rL"]))),col="blue",lwd=2)



# =================================== Parameter estimation by Hell optimisation - Expansion (6-24h) ============================

bes.par.kr <- c(kL=0.5,kS=0.5,rL=0.5,rS=0.5)


# Assume eL=eS=0
# Only estimate k and r
# Calculate the moments at 24h for L and S only, with initial conditions at t=6h given by data
bes.Hell.optim.kr.E <- powell(bes.par.kr, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.kr)
	par.i <- replace.par(bes.par.E,par)
	mom.i.E <- WITS.moment.sol.steps(bes.t.E-bes.t.B,par.i,bes.sim.mom.B,h=0.1)
	dist.E <- Hell.dist.o2p(bes.sim.mom.E,mom.i.E,2:3)
	#	print(c(par,div=dist.E))
	if(is.finite(dist.E)) dist.E else 1E100
},control=list(rhoend=1E-5))

cbind(bes.par.E[names(bes.Hell.optim.kr.E$par)],bes.Hell.optim.kr.E$par)


# ------------------------------------ Inter-sample variability ----------------------------------

bes.par.kr <- c(kL=0.5,kS=0.5,rL=0.6,rS=0.6)

# Generate 64 datasets of 100 WITS and estimate k,r parameters from each
bes.Hell.optim.E.bs <- foreach(i=1:64, .combine=rbind) %dopar% {
	sub.mom.B <- WITS.gillespie.moments(bes.sim.B[sample(1:nrow(bes.sim.B),100,replace = T),])
	sub.mom.E <- WITS.gillespie.moments(bes.sim.E[sample(1:nrow(bes.sim.E),100,replace = T),])
	sol <- powell(bes.par.kr, function(par){
		if(min(par)<0) return(1E100)
		names(par) <- names(bes.par.kr)
		par.i <- replace.par(bes.par.E,par)
		mom.i.E <- WITS.moment.sol.steps(bes.t.E-bes.t.B,par.i,sub.mom.B,h=0.2)
		dist.E <- Hell.dist.o2p(sub.mom.E,mom.i.E,2:3)
		if(is.finite(dist.E)) dist.E else 1E100
	},control=list(rhoend=1E-6))
	c(sol$par,D=sol$value)
}

summary(bes.Hell.optim.E.bs)
boxplot(bes.Hell.optim.E.bs,ylim=c(0,1.2))
points(1:4,bes.par.E[names(bes.par.kr)],pch="+",col="red",cex=2)



# =================================== Parameter estimation by Hell optimisation - Spillover (24-48h) ============================

xx <- 0

# Calculate the moments at 48h for all 3 organs, with initial conditions at t=24h given by data
bes.Hell.optim.S <- powell(bes.par.def*rnorm(8,4,0.4), function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.def)
	par.i <- replace.par(bes.par.S,par)
	mom.i.S <- WITS.moment.sol.steps(bes.t.S-bes.t.E,par.i,bes.sim.mom.E,h=1)
	dist.S <- Hell.dist.o2p(bes.sim.mom.S,mom.i.S,1:3)
	xx <<- xx+1
	if(xx%%20==1) print(c(par,div=dist.S),digits=4) 
	if(is.finite(dist.S)) dist.S else 1E100
},control=list(rhobeg=1,rhoend=1E-5))

round(cbind(bes.par.S,bes.Hell.optim.S$par),3)

bes.Hell.optim.S.2 <- powell(bes.Hell.optim.S$par, function(par){
	if(min(par)<0) return(1E100)
	names(par) <- names(bes.par.def)
	par.i <- replace.par(bes.par.S,par)
	mom.i.S <- WITS.moment.sol.steps(bes.t.S-bes.t.E,par.i,bes.sim.mom.E,h=0.5)
	dist.S <- Hell.dist.o2p(bes.sim.mom.S,mom.i.S,1:3)
	xx <<- xx+1
	if(xx%%20==1) print(c(par,div=dist.S),digits=4) 
	if(is.finite(dist.S)) dist.S else 1E100
},control=list(rhobeg=1,rhoend=1E-7))

round(cbind(bes.par.S,bes.Hell.optim.S$par,bes.Hell.optim.S.2$par),3)


# ------------------------------------ Inter-sample variability ----------------------------------

# Generate 32 datasets of 100 WITS and estimate all parameters from each
bes.Hell.optim.S.bs <- foreach(i=1:32, .combine=rbind) %dopar% {
	sub.mom.E <- WITS.gillespie.moments(bes.sim.E[sample(1:nrow(bes.sim.E),100,replace = T),])
	sub.mom.S <- WITS.gillespie.moments(bes.sim.S[sample(1:nrow(bes.sim.S),100,replace = T),])
	sol <- powell(bes.par.def*5, function(par){
		if(min(par)<0) return(1E100)
		names(par) <- names(bes.par.def)
		par.i <- replace.par(bes.par.S,par)
		mom.i.S <- WITS.moment.sol.steps(bes.t.S-bes.t.E,par.i,sub.mom.E,h=0.5)
		dist.S <- Hell.dist.o2p(sub.mom.S,mom.i.S,1:3)
		if(is.finite(dist.S)) dist.S else 1E100
	},control=list(rhoend=1E-7))
	c(sol$par,D=sol$value)
}

boxplot(bes.Hell.optim.S.bs,ylim=c(0,1.5))
points(1:8,bes.par.S[names(bes.par.def)],pch="+",col="red",cex=2)





#########################################################################################################################
# Save objects for report

save(list=ls()[grep("bes.",ls())], file="Moments/Hell_tests_BES_2.RData")


