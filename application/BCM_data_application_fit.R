###################################################################
# Code to do Bayesian splines estimate of biomarker effect in the setting of two arms, one experimental arm and one control arm #
# Time-to-event endpoint T with right censoring allowed  # 
# Assume L (>=1) continuous biomarkers x=(x_1, ..., x_L)' are available for each patient #
# Bayesian penalized splines is used to model the effect of biomarker x on endpoint T in terms of log HR in each arm, and the difference between two arms #
# The prognostic and predictive effects of multiple biomarkers are assumed to be additive # 
###################################################################
         
            
BCM_data_application_fit <- function(surdata,
                                     #####################################################################################################	 	
                                     ### surdata = a data frame with the observations arranged by row, and including the following columns:
                                     ### treatment arm indicator trt, entry date enter, event flag status (1=event, 0=censor), and event or censor date T_date. ###  
                                     ##########			              
                                     X,                                   # a data frame giving continuous marker values, with rows as patients and columns as biomarkers 
				     	                       n_events = 400,						  			  # the number of events expected to be observed at the end of the trial
                                     check = c(0.25, 0.50, 0.75, 1.00), 			   		# interim analysis times based on % of observed events for T
				     	                       maxIntKnots = 9,								      # the maximum number of interior knots allowed for penalized splines
				     	                       censor = FALSE,			                # whether marker values of right censored patients are included to decide the interior knots locations
				     	                       MCMC.specs = list(burnin=20000,B=2000,thin=20,multiplier=0.5),         		# a list giving the MCMC specs for spline modeling, including burn-in length, number of posterior samples, thinning parameter and initial scaling factor for the size of random walk proposal
				     	                       ref=apply(as.matrix(X), 2, median), 					  # the reference value for each marker in the control group, default to the column-wise median of X
				     	                       alpha = 0.1,                                           		    # the significance level at which an (1-alpha)*100% simultaneous/pointwise credible band is calculated
				     	                       marker.name,                                                   # the name of the marker used as part of the output file name
                                     path="./") {  	                                              	# the directory where the output should be saved
	 	                                              		
		
      
  ### Code here to load in the functions that will be called later
  source("util.R")
  source("Bayes_penalized_splines_general.R")  

  ### Number of analyses (including interim and final analysis)
	K <- length(check)


  ### Extract the columns from the data frame surdata
  trt <- surdata$trt
  enter <- surdata$enter
	status <- surdata$status
  T.date <- surdata$T_date 
      
  #############################################
  #############################################

  #############################################
  ###### Loop over interim analysis times #####
  #############################################
      
  for(j in 1:K){
    
    #### Censor data appropriately #### 
		event.date <- T.date[status==1]        
    date <- sort(event.date)[ceiling(n_events*check[j])]
         
    index <- which(enter < date)	       		# only consider those accrued by current date
    n.j <- length(index)			        	    # sample size at current date
    trt.j <- trt[index] 		       		      # trt group of those accrued by current date
    X.j <- X[index, , drop=FALSE]           # continuous marker values of those accrued by current date
 
         
    #### Time to Event Data T ####
    enter.j <- enter[index]			  		      # entry dates of those accured by current date
		status.j <- status[index]				        # event flag status (1=event, 0=censor) of those accrued by current date
    T.date.j <- T.date[index]           	  # event or censor dates of those accrued by current date 
    T.j <- T.date.j - enter.j		      	    # follow-up times of those accrued by current date        


	  #### Initialize time and censoring vectors ####
    T.cen <- rep(NA,n.j)					          # T.cen=1 for an observed event; =0 for right censoring.
    censorLimit.T <- rep(NA,n.j)				    # censorLimit.T is the time to event or right censoring (whichever happens first).

         
    for(k in 1:n.j){  					            # Loop over individuals on study at this interim time point (n.j set)
			# if event or censor occurs before interim j 
      if(T.date.j[k] <= date) {
        T.cen[k] <- status.j[k]
        censorLimit.T[k] <- T.j[k]
      }

			# if event or censor occurs after interim j
			else {
			  T.cen[k] <- 0
			  censorLimit.T[k] <- date - enter.j[k]
			}
    }
                  
    
    data.j = data.frame(enter=enter.j, futime=censorLimit.T, status=T.cen, trt=trt.j)


		##### Call the function to fit Bayesian penalized splines for the LHR as a function of marker based on data available at this time
    result.j = Bayes_penalized_splines(data=data.j[,c("futime","status","trt")], X=X.j, xrange=t(apply(as.matrix(X), 2, range)), 
                                       numIntKnots=pmin(round(sum(data.j$status)/20), maxIntKnots),  
                                       censor=censor, MCMC.specs=MCMC.specs, ref=ref, alpha=alpha)

	 	##### If the function call is a success, save the output at this analysis time #####
    if(result.j$success){
      
      model.fit = list(X=result.j$X,													                                            # matrix of marker values for patients accrued up to now
                       est_ctl=result.j$est_ctl, upr_ctl=result.j$upr_ctl, lwr_ctl=result.j$lwr_ctl,			# estimate and joint band for the LHR of the control group, relative to the reference value in the control group
				 	             upr_ctl_pt=result.j$upr_ctl_pt, lwr_ctl_pt=result.j$lwr_ctl_pt,					          # pointwise band for the LHR of the control group, relative to the reference value in the control group
				 	             prob_pos_ctl=result.j$prob_pos_ctl, prob_neg_ctl=result.j$prob_neg_ctl,            # posterior probability of being positive/negative for the LHR of the control group, relative to the reference value in the control group
			 	 	             est_exp=result.j$est_exp, upr_exp=result.j$upr_exp, lwr_exp=result.j$lwr_exp,			# estimate and joint band for the LHR of the experimental group, relative to the reference value in the control group
				 	             upr_exp_pt=result.j$upr_exp_pt, lwr_exp_pt=result.j$lwr_exp_pt,					          # pointwise band for the LHR of the experimental group, relative to the reference value in the control group
				 	             prob_pos_exp=result.j$prob_pos_exp, prob_neg_exp=result.j$prob_neg_exp,            # posterior probability of being positive/negative for the LHR of the experimental group, relative to the reference value in the control group
			   	 	           est_diff=result.j$est_diff, upr_diff=result.j$upr_diff, lwr_diff=result.j$lwr_diff,		# estimate and joint band for the LHR between two arms (predictive marker effect that is independent of reference value)
				 	             upr_diff_pt=result.j$upr_diff_pt, lwr_diff_pt=result.j$lwr_diff_pt,                # pointwise band for the LHR difference between two arms (predictive marker effect that is independent of reference value)
				 	             prob_pos_diff=result.j$prob_pos_diff, prob_neg_diff=result.j$prob_neg_diff)				# posterior probability of being positive/negative for the LHR between two arms (predictive marker effect that is independent of reference value)        

 
			acpt.rate = list(ctl=result.j$acpt_ctl, diff=result.j$acpt_diff)								          # acceptance rate of the random walk proposal for penalized spline coefficients in MCMC

			geweke.test = list(theta=result.j$geweke_theta, geweke_sigma2_g=result.j$geweke_sigma2_g, geweke_sigma2_h=result.j$geweke_sigma2_h)						# Geweke test statistics for the MCMC samples of each model parameter

      res = list(data=data.j, ref=ref, model.fit=model.fit, acpt.rate=acpt.rate, geweke.test=geweke.test) 	    # create a list with the components above

      saveRDS(res, file=paste0(path, marker.name, "_Bayes_fit_check_", j, ".Rds"))								              # write the output at this analysis time to the designated directory
    }

	 	else{

		 	print("WARNING: MCMC FAILS DUE TO NUMERICAL INSTABILITY")									                # print out a warning message to the screen if the function call fails
	 	
	 	}

      	        	   
    ######################################################                      

 }  # end loop over interim checks
      

} # end function



