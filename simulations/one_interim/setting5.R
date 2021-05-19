# Direct to the folder where the code is saved
setwd("./simulations")


# The following packages are needed to call the function
library(splines)
library(MASS)
library(survival)  
library(Matrix)
library(coda)


# Load in the function
source("Sim_Shell_Exponential.R")

# Delta=0.4
fnc2 <- function(x) {
  x2 <- exp(25*(x-0.5))
  y <- log(0.4)*x2/(1+x2)
  return(y)
}


# n = 500, marker prevalence threshold for dichotomization = 0.2, one interim analysis, Delta = 0.4
start_time = proc.time()

res <- Sim_Shell_Exponential(Nrep = 1000, seed = 500, n = 500, n_events = 400, m.prev.lwr = 0.2,							     
                             last.on.days = 365.25*2.5, check = c(0.50, 1.00), md.t = 365.25/2, percent.exp = 0.5, 
                             marker.specs = list(dist="beta", param=c(1,1), range=c(0,1), ref="median", censor=FALSE), 		                                   
                             FUN.ctl = function(x) {0*x}, FUN.exp = fnc2, 
                             MCMC.specs = list(burnin=25000,B=2000,thin=25,multiplier=0.35), maxIntKnots = 7, alpha = 0.1, 
                             p_eff_all_Bayes=rep(0.975,2), p_eff_grp_Bayes=rep(0.975,2), p_fut_all_Bayes=0.05, p_fut_grp_Bayes=0.05, 
                             m.increment = 0.01, p_int_freq = rep(0.05,2), 
                             p_eff_all_freq=rep(0.025,2), p_eff_grp_freq=rep(0.025,2), 
                             HR_fut_all_freq=1, HR_fut_grp_freq=1, stop_eff = TRUE, verbose=TRUE)

runtime = proc.time() - start_time
res[["runtime"]] = runtime

saveRDS(res, "./one_interim/output/setting5/trial_results.Rds")


