###################################################################
# Code to run clinical trial simulations in the setting of two arms, one experimental arm and one control arm #
# Time-to-event endpoint T (exponential) with right censoring allowed  # 
# Assume one continuous biomarker x is available for each patient #
# The distribution of marker x, and the treatment effect as a function of marker x are specified by the user # 
# Implement the Bayesian as well as a comparator frequentist adaptive design #
# If decide to dichotomize, record the marker status for each patient from each simulated trial#
###################################################################
         
            
Sim_Shell_Exponential <- function(Nrep = 1, seed = 12342423,		            # Number of iterations and seed to reproduce results
                                  n = 500, 							   	 	              # total number of patients expected to enroll (both arms combined)
                                  n_events = 400,						  		          # the number of events expected to be observed at the end of the trial
                                  m.prev.lwr = 0.25,						   		      # the lower bound allowed for the marker prevalence estimated from the splines fit
                                  last.on.days = 365*3,      				   	  	# date that the last patient enters (length of accural period in days)                   
                                  check = c(0.25, 0.50, 0.75, 1.00), 		  	# interim analysis times based on % of observed events for T
                                  md.t = 365,	                              # median of endpoint T (exponential) for baseline hazard of the control group 
                                  percent.exp = 0.5,                        # randomization ratio to the experimental group
                                  marker.specs = list(dist="beta",param=c(1,1),range=c(0,1),ref="median", censor=FALSE), 		# a list giving the distribution, range, reference point of continuous marker, and whether to include right censored subjects in calculation of interior knots location                                   
                                  FUN.ctl = function(x) {0*x},              # mapping function of the log hazard as a function of marker for the control group
                                  FUN.exp = function(x) {0*x},              # mapping function of the log hazard as a function of marker for the experimental group
                                  MCMC.specs = list(burnin=20000,B=2000,thin=20,multiplier=0.5),         		# Bayes: a list giving the MCMC specs for spline modeling, including burn-in length, number of posterior samples, thinning parameter and scaling factor for the size of random walk proposal
				     	                    maxIntKnots = 9,								          # Bayes: the maximum number of interior knots allowed for penalized splines
                                  alpha = 0.05,                             # Bayes: the significance level at which an (1-alpha)*100% simultaneous/pointwise credible band is calculated
				     	                    p_eff_all_Bayes = rep(0.975,4),  			    # Bayes: a vector of overall cohort posterior thresholds to claim efficacy (length=number of checks)
				     	                    p_eff_grp_Bayes = rep(0.9,4),  					  # Bayes: a vector of marker subgroup posterior thresholds to claim efficacy (length=number of checks)
				    	                    p_fut_all_Bayes = rep(0.05,3),						# Bayes: a vector of overall cohort posterior thresholds for early stopping due to futility (length=number of checks - 1)
				    	                    p_fut_grp_Bayes = rep(0.05,3),						# Bayes: a vector of marker subgroup posterior thresholds for early stopping due to futility (length=number of checks - 1)
				  	                      m.increment = 0.01, 							        # Freq: the increment of the sequence of marker cutpoint searched
				  	                      p_int_freq = rep(0.01,4),							    # Freq: a vector of minimum p-value cutoffs for the treatment-by-marker interaction to decide whether to dichotomize (length=number of checks)
				  	                      p_eff_all_freq = rep(0.025,4),  					# Freq: a vector of overall cohort one-sided log-rank test p-value cutoffs to claim efficacy (length=number of checks)
				  	                      p_eff_grp_freq = rep(0.1,4),  						# Freq: a vector of marker subgroup one-sided log-rank test p-value cutoffs to claim efficacy (length=number of checks)
				  	                      HR_fut_all_freq = rep(1,3),						    # Freq: a vector of overall cohort thresholds of HR above which to claim futility (length=number of checks - 1)
				  	                      HR_fut_grp_freq = rep(1,3),						    # Freq: a vector of marker subgroup thresholds of HR above which to claim futility (length=number of checks - 1)
				     	                    stop_eff = TRUE,                          # a logical indictor of whether to stop early for efficacy	
				    	                    verbose = FALSE) { 								        # a logical indictor of whether to print progress		 	                                              		

		
   set.seed(seed)
      
 
      
   ### Code here to load in the functions that will be called later
   source("simdata.R")
   source("Bayes_adaptive_design.R")   
   source("Freq_adaptive_design.R")


   ### Number of analyses (including interim and final analysis)
   K <- length(check)


   ### Create matrices and vectors to save results

   ### Bayesian design
   # Initialize indicator vector for failure of MCMC convergence in fitting Bayesian penalized splines for each simulated trial
   fail.MCMC <- rep(1, Nrep)

   # Initialize indicator matrices for efficacy/futility/continue(rather than stopping early) at each analysis
   stop.eff.pos.Bayes <- matrix(NA, nrow=Nrep, ncol=K)
   stop.eff.neg.Bayes <- matrix(NA, nrow=Nrep, ncol=K)
   stop.eff.all.Bayes <- matrix(NA, nrow=Nrep, ncol=K)

   stop.fut.pos.Bayes <- matrix(NA, nrow=Nrep, ncol=K-1)
   stop.fut.neg.Bayes <- matrix(NA, nrow=Nrep, ncol=K-1)
   stop.fut.all.Bayes <- matrix(NA, nrow=Nrep, ncol=K-1)

   continue.pos.Bayes <- matrix(NA, nrow=Nrep, ncol=K-1)
   continue.neg.Bayes <- matrix(NA, nrow=Nrep, ncol=K-1)
   continue.all.Bayes <- matrix(NA, nrow=Nrep, ncol=K-1)

   # Initialize vector to save final result for each simulated trial
   final.Bayes <- rep(NA, Nrep)

   # Initialize vector to save final sample sizes for each simulated trial
   sampsize.Bayes <- rep(NA, Nrep)
   sampsize.pos.Bayes <- rep(NA, Nrep)
   sampsize.neg.Bayes <- rep(NA, Nrep)
   
   # Initialize vector to save marker prevalence if dichotomized for each simulated trial
   marker.prevalence.Bayes <- rep(NA, Nrep)

   # Initialize matrix to store marker value and status of each patient if dichotomized, for each simulated trial
   marker.value.Bayes <- matrix(NA, nrow=Nrep, ncol=n)
   marker.status.Bayes <- matrix(NA, nrow=Nrep, ncol=n)


   ### Freq design
   # Initialize indicator matrices for efficacy/futility/continue(rather than stopping early) at each analysis
   stop.eff.pos.freq <- matrix(NA, nrow=Nrep, ncol=K)
   stop.eff.neg.freq <- matrix(NA, nrow=Nrep, ncol=K)
   stop.eff.all.freq <- matrix(NA, nrow=Nrep, ncol=K)

   stop.fut.pos.freq <- matrix(NA, nrow=Nrep, ncol=K-1)
   stop.fut.neg.freq <- matrix(NA, nrow=Nrep, ncol=K-1)
   stop.fut.all.freq <- matrix(NA, nrow=Nrep, ncol=K-1)

   continue.pos.freq <- matrix(NA, nrow=Nrep, ncol=K-1)
   continue.neg.freq <- matrix(NA, nrow=Nrep, ncol=K-1)
   continue.all.freq <- matrix(NA, nrow=Nrep, ncol=K-1)

   # Initialize vector to save final result for each simulated trial
   final.freq <- rep(NA, Nrep)

   # Initialize vector to save final sample sizes for each simulated trial
   sampsize.freq <- rep(NA, Nrep)
   sampsize.pos.freq <- rep(NA, Nrep)
   sampsize.neg.freq <- rep(NA, Nrep)
   
   # Initialize vector to save marker prevalence if dichotomized for each simulated trial
   marker.prevalence.freq <- rep(NA, Nrep)

   # Initialize matrix to store marker value and status of each patient if dichotomized, for each simulated trial
   marker.value.freq <- matrix(NA, nrow=Nrep, ncol=n)
   marker.status.freq <- matrix(NA, nrow=Nrep, ncol=n)

      
   #############################################
   #############################################

   ### MAIN LOOP OVER ITERATIONS:
   
   #############################################
   #############################################
   
   for (i in 1:Nrep) { 
      
      ##### Generate Patient Data #####
      res.data <- simdata(n=n, last.on.days=last.on.days, md.t=md.t, percent.exp=percent.exp, marker.specs=marker.specs, FUN.ctl=FUN.ctl, FUN.exp=FUN.exp)

      ##### Run Bayesian design #####
	    res.Bayes <- Bayes_adaptive_design(surdata=res.data$surdata, ref=res.data$ref, n_events=n_events, m.prev.lwr=m.prev.lwr, check=check, marker.specs=marker.specs,
	                                       MCMC.specs=MCMC.specs, maxIntKnots=maxIntKnots, alpha=alpha, p_eff_all=p_eff_all_Bayes, p_eff_grp=p_eff_grp_Bayes,
	                                       p_fut_all=p_fut_all_Bayes, p_fut_grp=p_fut_grp_Bayes, stop_eff=stop_eff)
	    
	    if(!res.Bayes$fail.MCMC) {	
	      fail.MCMC[i] <- 0

	      stop.eff.pos.Bayes[i,] <- res.Bayes$stop.eff.pos
	      stop.eff.neg.Bayes[i,] <- res.Bayes$stop.eff.neg
	      stop.eff.all.Bayes[i,] <- res.Bayes$stop.eff.all

	      stop.fut.pos.Bayes[i,] <- res.Bayes$stop.fut.pos
	      stop.fut.neg.Bayes[i,] <- res.Bayes$stop.fut.neg
	      stop.fut.all.Bayes[i,] <- res.Bayes$stop.fut.all  
						
	      continue.pos.Bayes[i,] <- res.Bayes$continue.pos
	      continue.neg.Bayes[i,] <- res.Bayes$continue.neg
	      continue.all.Bayes[i,] <- res.Bayes$continue.all

	      final.Bayes[i] <- res.Bayes$res.final

	      sampsize.Bayes[i] <- res.Bayes$sampsize
	      sampsize.pos.Bayes[i] <- res.Bayes$sampsize.pos
	      sampsize.neg.Bayes[i] <- res.Bayes$sampsize.neg
	      
	      marker.prevalence.Bayes[i] <- res.Bayes$marker.prevalence
	      marker.value.Bayes[i,] <- res.Bayes$x
	      marker.status.Bayes[i,] <- res.Bayes$marker.status
	    }


	    ##### Run Freq design #####	
	    res.freq <- Freq_adaptive_design(surdata=res.data$surdata, n_events=n_events, m.prev.lwr=m.prev.lwr, check=check, m.increment=m.increment, p_int=p_int_freq,
	                                     p_eff_all=p_eff_all_freq, p_eff_grp=p_eff_grp_freq, HR_fut_all=HR_fut_all_freq, HR_fut_grp=HR_fut_grp_freq, stop_eff=stop_eff)
	    
	    stop.eff.pos.freq[i,] <- res.freq$stop.eff.pos
	    stop.eff.neg.freq[i,] <- res.freq$stop.eff.neg
	    stop.eff.all.freq[i,] <- res.freq$stop.eff.all

	    stop.fut.pos.freq[i,] <- res.freq$stop.fut.pos
	    stop.fut.neg.freq[i,] <- res.freq$stop.fut.neg
	    stop.fut.all.freq[i,] <- res.freq$stop.fut.all  
						
	    continue.pos.freq[i,] <- res.freq$continue.pos
	    continue.neg.freq[i,] <- res.freq$continue.neg
	    continue.all.freq[i,] <- res.freq$continue.all
	    
	    final.freq[i] <- res.freq$res.final
	    
	    sampsize.freq[i] <- res.freq$sampsize
	    sampsize.pos.freq[i] <- res.freq$sampsize.pos
	    sampsize.neg.freq[i] <- res.freq$sampsize.neg

	    marker.prevalence.freq[i] <- res.freq$marker.prevalence
	    marker.value.freq[i,] <- res.freq$x
	    marker.status.freq[i,] <- res.freq$marker.status
	
	    if(verbose){
	      print(sprintf("Iter %d Completed", i))
	    } 
      
   }   #  end loop over trial (iteration)


   
   return(list(fail.MCMC=fail.MCMC, 
               stop.eff.pos.Bayes=stop.eff.pos.Bayes, stop.eff.neg.Bayes=stop.eff.neg.Bayes, stop.eff.all.Bayes=stop.eff.all.Bayes, 
               stop.fut.pos.Bayes=stop.fut.pos.Bayes, stop.fut.neg.Bayes=stop.fut.neg.Bayes, stop.fut.all.Bayes=stop.fut.all.Bayes,
               continue.pos.Bayes=continue.pos.Bayes, continue.neg.Bayes=continue.neg.Bayes, continue.all.Bayes=continue.all.Bayes,
               sampsize.Bayes=sampsize.Bayes, sampsize.pos.Bayes=sampsize.pos.Bayes, sampsize.neg.Bayes=sampsize.neg.Bayes, 
               marker.prevalence.Bayes=marker.prevalence.Bayes, final.Bayes=final.Bayes,
               marker.value.Bayes=marker.value.Bayes, marker.status.Bayes=marker.status.Bayes,
               stop.eff.pos.freq=stop.eff.pos.freq, stop.eff.neg.freq=stop.eff.neg.freq, stop.eff.all.freq=stop.eff.all.freq, 
               stop.fut.pos.freq=stop.fut.pos.freq, stop.fut.neg.freq=stop.fut.neg.freq, stop.fut.all.freq=stop.fut.all.freq,
               continue.pos.freq=continue.pos.freq, continue.neg.freq=continue.neg.freq, continue.all.freq=continue.all.freq,
               sampsize.freq=sampsize.freq, sampsize.pos.freq=sampsize.pos.freq, sampsize.neg.freq=sampsize.neg.freq, 
               marker.prevalence.freq=marker.prevalence.freq, final.freq=final.freq,
               marker.value.freq=marker.value.freq, marker.status.freq=marker.status.freq))   

} #end function



