##########################################################################################
#### Implement Bayesian penalized splines to estimate the log hazard ratio for endpoint T as a function of continuous marker x.	####
#### Two arms (one control group and one experimental group) are assumed. The LHR is relative to the reference point in the control group.	####
#### This function is called by Sim_Shell_Exponential.			####
##########################################################################################


Bayes_penalized_splines <- function(data, 				
                                    #####################################################################################################	 	
                                    ### data = a data frame with the observations arranged by row, and including the columns below (in the following order):
                                    ###        1) futime: time to event or right censoring, whichever happens first
                                    ###        2) status: the censoring status indicator (=1 for event; =0 for right censored)
                                    ###        3) x: the value of the continuous marker
                                    ###        4) trt: the group indicator (=1 for experimental group; =0 for control group)
                                    ##########
                                    xrange, 							    # a vector of length 2 which gives the range of marker values
                                    numIntKnots=9,						# the number of interior knots (should scale with number of events)
                                    censor=FALSE,						  # whether marker values of right censored patients are included to decide the interior knots locations
                                    MCMC.specs=list(burnin=10000,B=2000,thin=10,multiplier=0.5), 
                                    ##########
                                    ### MCMC.specs = a list of MCMC specifications including 
                                    ###        	1) burnin: the number of burnin samples
                                    ###        	2) B: the number of posterior samples
                                    ###        	3) thin: the thinning parameter
                                    ###        	4) multiplier: a scaling factor for the size of random walk proposal at the beginning of a MCMC chain
                                    ### Users might need to tune MCMC.specs$thin and MCMC.specs$multiplier to make sure that the acceptance probabilities and trace plots are reasonable
                                    ##########
                                    ref, 								     # the reference value for the marker in the control group
                                    alpha=0.05){						 # the significance level at which an (1-alpha)*100% simultaneous credible band is calculated
  ######################################################################################################### 

	### Code here to load in the functions that will be called later
	source("Get_splines.R")
	source("MCMC.R")
	source("jointband.R")


	# Compute the O'Sullivan splines design matrix for the continuous marker x, and the interaction term for x and trt
	res_splines <- Get_splines(data,xrange,numIntKnots,censor)


	# A numerical stability check for the spline design matrix
	# If not passing the stability check, stop the function call now
	if(res_splines$stabCheck)
	{
		return(list(x=sort(data$x), success=FALSE))
	}

	
	# If passing the stability check, perform MCMC
	else
	{
		data <- res_splines$data
		formula <- res_splines$formula
		ncolz <- res_splines$ncolz 

		for(count in 1:15)
		{			
			res_MCMC <- MCMC(data, formula, ncolz, prior=list(alpha=0.01, beta=0.01),
					     multiplier1=MCMC.specs$multiplier, multiplier2=0.5, multiplier3=0.25,		# scaling factors to adaptively adjust the size of random walk proposal
					     burnin=MCMC.specs$burnin, B=MCMC.specs$B, thin=MCMC.specs$thin)
			
			if(res_MCMC$success) { 
				break
			}
		}

		if(count==15)
		{
			return(list(x=sort(data$x), success=FALSE))
		}

					
		# Using the posterior samples of the spline coefficients,
		# get the posterior samples of the log hazard ratio relative to the control group at the reference value,
		# as a function of marker value for the control group.
		# This approach results in pinching in the credible band at the reference value. 
		data_g <- data[order(data$x),3:(ncolz+3)]
		idx <- which.min(abs(data_g$x-ref))
		data_g <- sweep(data_g,2,as.numeric(data_g[idx,]))

		posterior_g <- res_MCMC$theta[,2:(ncolz+2)] %*% t(data_g)
		simCR <- jointband(posterior_g,alpha)
		upr_g <- simCR$upr_CI
		lwr_g <- simCR$lwr_CI
		upr_g_pt <- apply(posterior_g, 2, function(x){quantile(x,1-alpha/2)})
		lwr_g_pt <- apply(posterior_g, 2, function(x){quantile(x,alpha/2)})


		# Using the posterior samples of the spline coefficients,
		# get the posterior samples of the difference in log hazard between the two arms, as a function of marker value.
		# Note that the group difference does not depend on the choice of the reference value
		data_h <- data[order(data$x),3:(ncolz+3)]
		data_h <- cbind(1,data_h)

		posterior_h <- res_MCMC$theta[,-c(1:(ncolz+2))] %*% t(data_h)
		simCR <- jointband(posterior_h,alpha)
		upr_h <- simCR$upr_CI
		lwr_h <- simCR$lwr_CI
		upr_h_pt <- apply(posterior_h, 2, function(x){quantile(x,1-alpha/2)})
		lwr_h_pt <- apply(posterior_h, 2, function(x){quantile(x,alpha/2)})

  
		# Using the posterior samples of the spline coefficients,
		# get the posterior samples of the log hazard ratio relative to the control group at the reference value,
		# as a function of marker value for the experimental group.
		posterior_exp <- posterior_g + posterior_h
		simCR <- jointband(posterior_exp,alpha)
		upr_exp <- simCR$upr_CI
		lwr_exp <- simCR$lwr_CI    
		upr_exp_pt <- apply(posterior_exp, 2, function(x){quantile(x,1-alpha/2)})
		lwr_exp_pt <- apply(posterior_exp, 2, function(x){quantile(x,alpha/2)})


		# return a list with the MCMC samples, posterior mean estimate and credible band for the log hazard ratio as a function of marker value
		return(list(x=sort(data$x), success=TRUE,										# sorted marker values x and indicator of success for the function call
			 theta=res_MCMC$theta, sigma2=res_MCMC$sigma2,								# MCMC samples of the model parameters in Bayesian penalized splines
			 est_ctl=apply(posterior_g,2,mean), upr_ctl=upr_g, lwr_ctl=lwr_g,					# estimate and joint band for the LHR of the control group
			 upr_ctl_pt=upr_g_pt, lwr_ctl_pt=lwr_g_pt,								# pointwise band for the LHR of the control group
			 est_diff=apply(posterior_h,2,mean), upr_diff=upr_h, lwr_diff=lwr_h,				# estimate and joint band for the difference in LHR between two arms
			 upr_diff_pt=upr_h_pt, lwr_diff_pt=lwr_h_pt,								# pointwise band for the LHR difference between two arms
			 est_exp=apply(posterior_exp,2,mean), upr_exp=upr_exp, lwr_exp=lwr_exp,				# estimate and joint band for the LHR of the experimental group
			 upr_exp_pt=upr_exp_pt, lwr_exp_pt=lwr_exp_pt,								# pointwise band for the LHR of the experimental group
			 MCMC_ctl=posterior_g, MCMC_diff=posterior_h, MCMC_exp=posterior_exp,				# MCMC samples of the LHR as a function of marker
			 acpt_ctl=res_MCMC$acpt_g, acpt_diff=res_MCMC$acpt_h, 						# acceptance rates for two MCMC chains
			 geweke_theta=res_MCMC$geweke_theta, geweke_sigma2=res_MCMC$geweke_sigma2 ))			# Geweke test statistics 

	}
}


   
