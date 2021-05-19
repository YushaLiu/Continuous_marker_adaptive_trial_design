##########################################################################################
#### Implement Bayesian penalized splines to estimate the log hazard ratio for endpoint T as a function of continuous marker x=(x_1, ..., x_J)'.####
#### Two arms (one control group and one experimental group) are assumed.	####
#### The prognostic and predictive effects of multiple markers are assumed to be additive. ####
#### For prognostic marker effects, the LHR is relative to the reference point in the control group specified by the user. ####
#### For predictive marker effects, the LHR does not depend on the choice of the reference point. #### 
#### This function is called by Bayes_adaptive_design. ####
##########################################################################################


Bayes_penalized_splines <- function(data, 				
		                                #####################################################################################################	 	
		                                ### data = a data frame with the observations arranged by row, and including the columns (in the following order):
		                                ###        1) futime: time to event or right censoring, whichever happens first
		                                ###        2) status: the censoring status indicator (=1 for event; =0 for right censored)
		                                ###        3) trt: the group indicator (=1 for experimental group; =0 for control group)
		                                ##########
		                                X,                # a data frame giving continuous marker values, with rows as patients and columns as biomarkers 
						                        xrange=cbind(rep(0, ncol(X)), rep(1, ncol(X))), 					# a matrix with 2 columns which gives the range of marker values, one marker per row
						                        numIntKnots=5,		# the number of interior knots (should scale with number of events)
						                        censor=FALSE,			# whether marker values of right censored patients are included to decide the interior knots locations
						                        MCMC.specs=list(burnin=10000, B=2000, thin=10, multiplier=0.5), 
		                                ##########
		                                ### MCMC.specs = a list of MCMC specifications including 
		                                ###        	1) burnin: the number of burnin samples
		                                ###        	2) B: the number of posterior samples
		                                ###        	3) thin: the thinning parameter
		                                ###        	4) multiplier: a scaling factor for the size of random walk proposal at the beginning of a MCMC chain
		                                ### Users might need to tune thin and multiplier to make sure that the acceptance probabilities and trace plots are reasonable
		                                ##########
						                        ref=apply(as.matrix(X), 2, median), 					# the reference value for each marker in the control group, default to the column-wise median of X
						                        alpha=0.05) {				# the significance level at which an (1-alpha)*100% simultaneous credible band is calculated
		                                ######################################################################################################### 

	### Code here to load in the functions that will be called later
	source("Get_splines_general.R")
	source("MCMC_general.R")
	source("util.R")


	# Compute the O'Sullivan splines design matrix for each continuous marker, and the interaction term between each marker and trt
	res_splines <- Get_splines(data, X, xrange, numIntKnots, censor)


	# A numerical stability check for the spline design matrix
	# If not passing the stability check, stop the function call now
	if(any(res_splines$stabCheck)){
		return(list(X=X, success=FALSE))
	}

	
	# If passing the stability check, perform MCMC
	else{
	  data <- res_splines$data
		formula <- res_splines$formula
		ncolz <- res_splines$ncolz 

		for(count in 1:10){			
			res_MCMC <- MCMC(data=data, formula=formula, ncolz=ncolz, prior=list(alpha=0.01, beta=0.01),
			                 multiplier1=MCMC.specs$multiplier, multiplier2=0.5, multiplier3=0.25,		# scaling factors to adaptively adjust the size of random walk proposal
			                 burnin=MCMC.specs$burnin, B=MCMC.specs$B, thin=MCMC.specs$thin)
			
			if(res_MCMC$success){ 
				break
			}
		}

		if(count==10){
			return(list(X=X, success=FALSE))
		}

		
		### Get the posterior samples of the prognostic and predictive marker effects using posterior samples of the spline coefficients
		data.X <- as.matrix(data[, -c(1:2)])
		theta <- res_MCMC$theta[, -c(1)]
		
		# Get the respective column indices of design matrices related to prognostic and predictive effects 
		for(j in 1:ncol(X)){
		  if(j==1){
		    # prognostic marker effects
		    idx.jg.start <- 1
		    idx.jg.end <- idx.jg.start + ncolz
		    idx.g <- idx.jg.start:idx.jg.end
		    
		    # predictive marker effects
		    idx.jh.start <- idx.jg.end + 1
		    idx.jh.end <- idx.jh.start + ncolz + 1
		    idx.h <- idx.jh.start:idx.jh.end
		  }
		  
		  else{
		    # prognostic marker effects
		    idx.jg.start <- idx.jh.end + 1
		    idx.jg.end <- idx.jg.start + ncolz
		    idx.g <- c(idx.g, idx.jg.start:idx.jg.end)
		    
		    # predictive marker effects
		    idx.jh.start <- idx.jg.end + 1
		    idx.jh.end <- idx.jh.start + ncolz
		    idx.h <- c(idx.h, idx.jh.start:idx.jh.end)
		  }
		}
		
		
		# Using the posterior samples of the spline coefficients,
		# get the posterior samples of the log hazard ratio relative to the control group at the reference value 
		# as a function of marker value for the control group.
		# This approach results in pinching in the credible band at the reference value. 
		data_g <- data.X[, idx.g]
		idx.ref <- which.min(dist_marker(ref, X))
		data_g_cen <- sweep(data_g, 2, as.numeric(data_g[idx.ref,]))

		posterior_g <- theta[, idx.g] %*% t(data_g_cen)
		simCR <- jointband(posterior_g, alpha)
		upr_g <- simCR$upr_CI
		lwr_g <- simCR$lwr_CI
		upr_g_pt <- apply(posterior_g, 2, function(x){quantile(x,1-alpha/2)})
		lwr_g_pt <- apply(posterior_g, 2, function(x){quantile(x,alpha/2)})
		prob_pos_g <- apply(posterior_g, 2, function(x){mean(x>0)})
		prob_neg_g <- apply(posterior_g, 2, function(x){mean(x<0)})


		# Using the posterior samples of the spline coefficients,
		# get the posterior samples of the difference in log hazard between the two arms, as a function of marker value.
		# Note that the group difference does not depend on the choice of the reference value
		data_h <- cbind(1, data.X[, idx.g])
		
		posterior_h <- theta[, idx.h] %*% t(data_h)
		simCR <- jointband(posterior_h, alpha)
		upr_h <- simCR$upr_CI
		lwr_h <- simCR$lwr_CI
		upr_h_pt <- apply(posterior_h, 2, function(x){quantile(x,1-alpha/2)})
		lwr_h_pt <- apply(posterior_h, 2, function(x){quantile(x,alpha/2)})
		prob_pos_h <- apply(posterior_h, 2, function(x){mean(x>0)})
		prob_neg_h <- apply(posterior_h, 2, function(x){mean(x<0)})
		
  
		# Using the posterior samples of the spline coefficients,
		# get the posterior samples of the log hazard ratio relative to the control group at the reference value,
		# as a function of marker value for the experimental group.
		posterior_exp <- posterior_g + posterior_h
		simCR <- jointband(posterior_exp, alpha)
		upr_exp <- simCR$upr_CI
		lwr_exp <- simCR$lwr_CI    
		upr_exp_pt <- apply(posterior_exp, 2, function(x){quantile(x,1-alpha/2)})
		lwr_exp_pt <- apply(posterior_exp, 2, function(x){quantile(x,alpha/2)})
		prob_pos_exp <- apply(posterior_exp, 2, function(x){mean(x>0)})
		prob_neg_exp <- apply(posterior_exp, 2, function(x){mean(x<0)})


		# return a list with the MCMC samples, posterior mean estimate and credible band for the log hazard ratio as a function of marker value
		return(list(X=X,                            # matrix of marker values
		            success=TRUE,										# indicator of success for the function call
		            theta=res_MCMC$theta, sigma2_g=res_MCMC$sigma2_g,	sigma2_h=res_MCMC$sigma2_h,						# MCMC samples of the model parameters in Bayesian penalized splines
		            est_ctl=apply(posterior_g,2,mean), upr_ctl=upr_g, lwr_ctl=lwr_g,					# estimate and joint band for the LHR of the control group
		            prob_pos_ctl=prob_pos_g, prob_neg_ctl=prob_neg_g,                         # posterior probability of being positive/negative for the LHR
		            upr_ctl_pt=upr_g_pt, lwr_ctl_pt=lwr_g_pt,								  # pointwise band for the LHR of the control group
		            est_diff=apply(posterior_h,2,mean), upr_diff=upr_h, lwr_diff=lwr_h,				# estimate and joint band for the difference in LHR between two arms
		            upr_diff_pt=upr_h_pt, lwr_diff_pt=lwr_h_pt,								# pointwise band for the LHR difference between two arms
		            prob_pos_diff=prob_pos_h, prob_neg_diff=prob_neg_h,                       # posterior probability of being positive/negative for the LHR
		            est_exp=apply(posterior_exp,2,mean), upr_exp=upr_exp, lwr_exp=lwr_exp,				# estimate and joint band for the LHR of the experimental group
		            upr_exp_pt=upr_exp_pt, lwr_exp_pt=lwr_exp_pt,							# pointwise band for the LHR of the experimental group
		            prob_pos_exp=prob_pos_exp, prob_neg_exp=prob_neg_exp,                     # posterior probability of being positive/negative for the LHR
		            MCMC_ctl=posterior_g, MCMC_diff=posterior_h, MCMC_exp=posterior_exp,			# MCMC samples of the LHR as a function of marker
		            acpt_ctl=res_MCMC$acpt_g, acpt_diff=res_MCMC$acpt_h, 			# acceptance rates for two MCMC chains
		            geweke_theta=res_MCMC$geweke_theta, geweke_sigma2_g=res_MCMC$geweke_sigma2_g, geweke_sigma2_h=res_MCMC$geweke_sigma2_h))			# Geweke test statistics 
	 }
}


   
