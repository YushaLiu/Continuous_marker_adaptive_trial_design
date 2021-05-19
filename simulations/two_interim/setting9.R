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
  if(x<=0.5) {
    log(0.4)*exp(30*(x-0.3))/(1+exp(30*(x-0.3)))
  }
  else {
    -log(0.4)*exp(30*(x-0.7))/(1+exp(30*(x-0.7)))+log(0.4)
  }
}


# n = 500, marker prevalence threshold for dichotomization = 0.2, two interim analyses, Delta = 0.4
start_time = proc.time()

res <- Sim_Shell_Exponential(Nrep = 1000, seed = 500, n = 500, n_events = 400, m.prev.lwr = 0.2,							     
                             last.on.days = 365.25*2.5, check = c(0.50, 0.75, 1.00), md.t = 365.25/2, percent.exp = 0.5,	     				                                            	                          
                             marker.specs = list(dist="beta", param=c(1,1), range=c(0,1), ref="median", censor=FALSE), 		                                   
                             FUN.ctl = function(x) {0*x}, FUN.exp = fnc2, 
                             MCMC.specs = list(burnin=25000,B=2000,thin=25,multiplier=0.35), maxIntKnots = 7, alpha = 0.1, 
                             p_eff_all_Bayes=rep(0.975,3), p_eff_grp_Bayes=c(0.975,0.99,0.99),  
                             p_fut_all_Bayes=rep(0.05,2), p_fut_grp_Bayes=rep(0.05,2), 
                             m.increment = 0.01, p_int_freq = rep(0.05,3), 
                             p_eff_all_freq=rep(0.025,3), p_eff_grp_freq=c(0.025,0.01,0.01), 
                             HR_fut_all_freq=rep(1,2), HR_fut_grp_freq=rep(1,2), stop_eff = TRUE, verbose=TRUE)

runtime = proc.time() - start_time
res[["runtime"]] = runtime

saveRDS(res, "./two_interim/output/setting9/trial_results.Rds")


