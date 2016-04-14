#############################
#   SALMONELLA WITS MODEL
#     MOMENTS EQUATION
#          TESTS
#############################

# MOMENTS COMPUTATION BENCHMARKS

source("Moments/WITS_moments.R")


library(foreach)
library(doParallel)
registerDoParallel(cores=16)

library(ggplot2)
library(dplyr)
library(tidyr)

library(Rcpp)
library(ellipse)



#====================== Model: Parameter values and simulations ==========================

par.test <- c(cL=1,cS=1,eL=0.1,eS=0.1,kL=0.5,kS=0.5,rL=0.6,rS=0.6)
# M.0
inoc.test <- 100
init.test.fun <- function(){c(rpois(1,inoc.test),0,0)}
t.try <- c(24,48,72)

# Simulate 10,000 WITS, drawing the inoculum size from a Poisson distribution: E(NB) = V(NB)
sim.test <- lapply(t.try, function(dur) WITS.gillespie.sim(dur,par.test,init.test.fun,1E4))
sim.test.mom <- t(sapply(sim.test,WITS.gillespie.moments))

save(sim.test,file="Moments/Benchmark_sim.RData")

#sim.test.mom <- cbind(data.frame(met=rep("Gillespie",4),h=rep(NA,4),t=t.try, elapsed=rep(NA,4)),t(sapply(sim.test,WITS.gillespie.moments)))

#system.time(foreach(i=1:1000, .combine=cbind) %dopar% single_wits_gillespie_r(0,72,c(rpois(1,init.test[1]),0,0),par.test,c(10000000,10000000,10000000)))

# ============================================ Benchmark for direct moment calculation =====================================

# Functions
fun.try <- list(WITS.moment.sol.1,WITS.moment.sol.2,WITS.moment.sol.4,WITS.moment.sol.5)
fun.lab <- c("vec.lin", "vec.cub", "mat.lin","mat.cub")
# Intervals for integral computation
dt.try <- c(1,6,12)
# expm methods
met.try <- c("Higham08.b","Ward77","AlMohy-Hi09","Pade","Taylor")
# Step size for integration
h.try <- c(5E-3,5E-2,5E-1)

par.try <- expand.grid(1:length(fun.try),met.try,h.try,t.try,dt.try)
colnames(par.try) <- c("version","met","h","t","dt")

# Calculate the moments with method V2
moments.benchmark <- foreach(i=1:nrow(par.try), .combine=rbind) %dopar% {
#	print(par.try[i,])
	elapsed <- system.time(test <- WITS.moment.sol.steps(par.try$t[i],par.test,init.test,interval=par.try$dt[i],h=par.try$h[i],solver=fun.try[[par.try$version[i]]],exp.met=as.character(par.try$met[i])))["elapsed"]
	sim.ref <- sim.test.mom[which(abs(t.try-par.try$t[i])<0.01),]
	c(elapsed,(test-sim.ref)/sim.ref)
}
expm.moments.benchmark <- cbind(par.try,integration=fun.lab[par.try$version],moments.benchmark)

expm.moments.benchmark$max.err <- apply(abs(expm.moments.benchmark[,8:16]), 1,max) 

save(mom.names, par.test, init.test,t.try,dt.try,expm.moments.benchmark, file="Moments/expm_benchmark.RData")



# ====================================== Plots =================================================================

# # Time elapsed
# expm.moments.benchmark %>% ggplot(aes(integration,elapsed)) + geom_point(aes(col=as.factor(h), shape=met),size=4) + facet_grid(t~dt) + scale_y_log10(breaks=as.numeric(c(1,2,5) %o% 10^(-2:1)))
# 
# # Relative error
# expm.moments.benchmark %>% ggplot(aes(integration,log10(max.err))) + geom_point(aes(col=as.factor(h), shape=met),size=4) + facet_grid(t~dt)
# 
# # Error vs. Time
# expm.moments.benchmark %>% filter(met=="AlMohy-Hi09") %>% ggplot(aes(max.err,elapsed)) + geom_point(aes(col=as.factor(h), shape=integration),size=4) + facet_grid(t~dt) + scale_y_log10(breaks=as.numeric(c(1,2,5) %o% 10^(-2:1))) + scale_x_log10()


