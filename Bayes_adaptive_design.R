###################################################################
# Code to do Bayesian adaptive design in the setting of two arms, one experimental arm and one control arm #
# Time-to-event endpoint T with right censoring allowed  # 
# Assume one continuous biomarker x is available for each patient #
# Bayesian penalized splines is used to model the effect of biomarker x on endpoint T in terms of log HR in each arm, and the difference between two arms #
# The pointwise credible bands from Bayesian penalized splines fit are used to dichotomize the marker #
# Interim and final analysis for efficacy and futility are performed based on the posterior thresholds specified by the user #
# If decide to dichotomize, record the marker status for each patient #
###################################################################
         
            
Bayes_adaptive_design <- function(surdata,
				#####################################################################################################	 	
				### surdata = a data frame with the observations arranged by row, and including the columns below (in the following order) ###
				### treatment arm indicator trt, entry dates enter, continuous marker values x, event times T, and observed TTE dates T_date. ###  
				##########		                   
                                  ref, 							   			# the reference point of the continuous marker
				     	    n_events = 400,						  		# the number of events expected to be observed at the end of the trial
				     	    m.prev.lwr = 0.25,						   		# the lower bound allowed for the marker prevalence estimated from the splines fit
                                  check = c(0.25, 0.50, 0.75, 1.00), 			   		# interim analysis times based on % of observed events for T
                                  marker.specs = list(dist="beta",param=c(1,1),range=c(0,1),ref="median", censor=FALSE), 		# a list giving the distribution, range, reference point of continuous marker, and whether to include right censored subjects in calculation of interior knots location                                   
                                  MCMC.specs = list(burnin=20000,B=2000,thin=20,multiplier=0.5),         		# a list giving the MCMC specs for spline modeling, including burn-in length, number of posterior samples, thinning parameter and scaling factor for the size of random walk proposal
				     	    maxIntKnots = 9,								# the maximum number of interior knots allowed for penalized splines
                                  alpha = 0.05,                                              	# the significance level at which an (1-alpha)*100% simultaneous/pointwise credible band is calculated
				     	    p_eff_all = rep(0.99,4),  						# a vector of overall cohort posterior thresholds to claim efficacy (length=number of checks)
				     	    p_eff_grp = rep(0.9,4),  							# a vector of marker subgroup posterior thresholds to claim efficacy (length=number of checks)
				    	    p_fut_all = rep(0.05,3),							# a vector of overall cohort posterior thresholds for early stopping due to futility (length=number of checks - 1)
				    	    p_fut_grp = rep(0.05,3),							# a vector of marker subgroup posterior thresholds for early stopping due to futility (length=number of checks - 1)
				     	    stop_eff = TRUE) { 								# a logical indictor of whether to stop early for efficacy		 	                                              		
		
      
   	### Code here to load in the functions that will be called later
   	source("Bayes_penalized_splines.R")   


   	### Number of analyses (including interim and final analysis)
   	K <- length(check)


   	### Create vectors to save results
   	# Initialize indicator vectors for stopping (overall, cohort), endpoint, and reasons
   	stop.eff.pos <- rep(0, K)
   	stop.eff.neg <- rep(0, K)
   	stop.eff.all <- rep(0, K)

   	stop.fut.pos <- rep(0, K-1)
   	stop.fut.neg <- rep(0, K-1)
   	stop.fut.all <- rep(0, K-1)

   	# Initialize indicator vectors for continuing (rather than stopping early)
   	continue.pos <- rep(0, K-1)
   	continue.neg <- rep(0, K-1)
   	continue.all <- rep(0, K-1)


   	# Initialize scalar to store final Bayesian sample sizes 
   	sampsize <- NA
   	sampsize.pos <- NA
   	sampsize.neg <- NA
   
   	# Initialize scalar to store marker prevalence and marker status for each patient if dichotomized
	marker.prevalence <- NA
   	marker.status <- rep(NA, nrow(surdata))


      ### Extract the columns from the data frame surdata
      trt <- surdata$trt
      enter <- surdata$enter
      x <- surdata$x 
      T.date <- surdata$T_date 
      

      ### Initialize the trial stopping indicator (of either cohort if marker is dichotomized), which is updated at each interim check 
      stop.any <- FALSE
      stop.pos <- FALSE
      stop.neg <- FALSE

      
   	#############################################
   	#############################################

      #############################################
      ###### Loop over interim analysis times #####
      #############################################
      
      for(j in 1:K){
      	         
         #####################################################
         #####################################################
         ################# Interim Analyses ##################
         #####################################################
         #####################################################

            
         #####################################################
         ###### Add estimation / interim decision code here ##### 

	   if(!stop.any) {

         	#### Censor data appropriately #### 
         
         	date <- sort(T.date)[ceiling(n_events*check[j])]

         
         	index <- which(enter < date)	       		# only consider those accrued by current date
         	enter.j <- enter[index]			  		# accrual dates of those accured by current date
         	n.j <- length(index)			        	# sample size at current date
         	trt.j <- trt[index] 		       		# trt group of those accrued by current date
         	x.j <- x[index]                          		# continuous marker value of those accrued by current date
 
         
        	#### Time to Event Data T ####
         	T.date.j <- T.date[index]           	     	# event dates of those accrued by current date
         	T.j <- T.date.j - enter.j		      	# event times of those accrued by current date
         

	   	#### Initialize time and censoring vectors ####
         	T.cen <- rep(1,n.j)					# T.cen=1 for an observed event; =0 for right censoring.
         	censorLimit.T <- rep(NA,n.j)				# censorLimit.T is the time to event or right censoring (whichever happens first).

         
         	for(k in 1:n.j){  					# Loop over individuals on study at this interim time point (n.j set)
                  
			#  at interim j for those accrued
            	if(T.date.j[k] <= date) {
            		# T.cen[k] remains 1
            		censorLimit.T[k] <- T.j[k]
            	}

			else {
            		T.cen[k] <- 0
            		censorLimit.T[k] <- date - enter.j[k]
            	}
         	}
                  
	   
         	data.j = data.frame(enter=enter.j, futime=censorLimit.T, status=T.cen, x=x.j, trt=trt.j)


		##### Call the function to fit Bayesian penalized splines for the LHR as a function of marker based on data available at this time
         	result.j = Bayes_penalized_splines(data=data.j[,c("futime","status","x","trt")], xrange=marker.specs$range, numIntKnots=pmin(round(sum(data.j$status)/20), maxIntKnots),  
						           censor=marker.specs$censor, MCMC.specs=MCMC.specs, ref=ref, alpha=alpha)

		##### If the function call fails due to numerical instability, terminate the current trial (iteration)
		if(!result.j$success) {
			return(list(fail.MCMC=TRUE))
		} 	


	   	##### If the function call is a success, dichotomize the marker or not based on the splines fit, and do interim or final analysis #####
         	else {

			##### Check whether to dichotomize the marker based on the splines fit #####
			fit.diff <- data.frame(x=result.j$x, est=result.j$est_diff, upr=result.j$upr_diff, lwr=result.j$lwr_diff,
					           upr_pt=result.j$upr_diff_pt, lwr_pt=result.j$lwr_diff_pt)

				
			##### Loop over all subjects to determine marker status based on the dichotomization rule from the splines fit #####
			marker <- rep(NA, length(x))
			marker.j <- rep(NA, length(x.j))

			if(any(fit.diff$upr_pt<0)) {

				# define marker positivity
				marker_pos <- fit.diff[fit.diff$upr_pt<0,]$x						
				marker_neg <- fit.diff[fit.diff$upr_pt>=0,]$x

				if(length(marker_neg)==0) {
					marker <- rep(1, length(x))
					marker.j <- rep(1, length(x.j))
				}

				else {
					for(k in 1:length(x)) {  								
						marker[k] <- min(abs(x[k]-marker_pos)) <= min(abs(x[k]-marker_neg))
					}

					for(k in 1:length(x.j)) {
						marker.j[k] <- as.numeric(x.j[k] %in% marker_pos)
					}
           		 	}
			}

			else {
				marker <- rep(0, length(x))
				marker.j <- rep(0, length(x.j))
			}


			##### If the estimated marker prevalence based on currently enrolled patients is in (m.prev.lwr, 1-m.prev.lwr), dichotomize the marker #####
			if( mean(marker.j==1)>m.prev.lwr & mean(marker.j==1)<1-m.prev.lwr ) {

				# dichotomize the continuous marker for individuals enrolled at this time
				data.j$marker <- marker.j

				samp.pos <- sum(data.j$marker==1)
				samp.neg <- sum(data.j$marker==0)


				# subset the data available at this time based on marker positivity and arm
				data.pos.exp <- data.j[(data.j$marker==1 & data.j$trt==1), c("enter","futime","status")]
				data.pos.ctl <- data.j[(data.j$marker==1 & data.j$trt==0), c("enter","futime","status")]
				data.neg.exp <- data.j[(data.j$marker==0 & data.j$trt==1), c("enter","futime","status")]
				data.neg.ctl <- data.j[(data.j$marker==0 & data.j$trt==0), c("enter","futime","status")]

 
				# simulate hazard rates for each marker cohort and arm
				hazard.pos.exp <- rgamma(100000, shape=1e-2+sum(data.pos.exp$status), rate=1+sum(data.pos.exp$futime) )
				hazard.pos.ctl <- rgamma(100000, shape=1e-2+sum(data.pos.ctl$status), rate=1+sum(data.pos.ctl$futime) )
				hazard.neg.exp <- rgamma(100000, shape=1e-2+sum(data.neg.exp$status), rate=1+sum(data.neg.exp$futime) )
				hazard.neg.ctl <- rgamma(100000, shape=1e-2+sum(data.neg.ctl$status), rate=1+sum(data.neg.ctl$futime) )


				###### check efficacy independently within each marker cohort ######
				eff.pos <- mean(hazard.pos.exp < hazard.pos.ctl)		# the marker positive cohort
				eff.neg <- mean(hazard.neg.exp < hazard.neg.ctl)		# the marker negative cohort


            		###### If this is the final check, skip the futility analysis ######
				if(j==K) {
					stop.eff.pos[j] <- as.numeric(eff.pos >= p_eff_grp[j])
					stop.eff.neg[j] <- as.numeric(eff.neg >= p_eff_grp[j])
					sampsize <- nrow(data.j)
					sampsize.pos <- samp.pos
					sampsize.neg <- samp.neg
					marker.prevalence <- mean(marker==1)
					marker.status <- marker
				}

				##### If this is not the final check, perform efficacy and futility analysis #####	
				else {
					if(stop_eff & (eff.pos >= p_eff_grp[j])) {				# stop the marker positive cohort due to efficacy
						stop.pos <- TRUE 
						stop.eff.pos[j] <- 1
					}

					if(stop_eff & (eff.neg >= p_eff_grp[j])) {				# stop the marker negative cohort due to efficacy
						stop.neg <- TRUE 
						stop.eff.neg[j] <- 1
					}
		
					if(stop.eff.pos[j]==0 | stop.eff.neg[j]==0) {				# if either marker cohort fails to stop early for efficacy, check futility

						x.rem <- x[-index]
						enter.rem <- enter[-index]
						trt.rem <- trt[-index]
						marker.rem <- marker[-index]

						enter.rem.pos.exp <- enter.rem[marker.rem==1 & trt.rem==1]
						enter.rem.pos.ctl <- enter.rem[marker.rem==1 & trt.rem==0]
						enter.rem.neg.exp <- enter.rem[marker.rem==0 & trt.rem==1]
						enter.rem.neg.ctl <- enter.rem[marker.rem==0 & trt.rem==0]

						n.rem.pos.exp <- sum(marker.rem==1 & trt.rem==1)
						n.rem.pos.ctl <- sum(marker.rem==1 & trt.rem==0)
						n.rem.neg.exp <- sum(marker.rem==0 & trt.rem==1)
						n.rem.neg.ctl <- sum(marker.rem==0 & trt.rem==0)
		
						win.pos <- 0
						win.neg <- 0

						# futility analysis loops	
						for(l in 1:2000) {

							# Generate new time to event data 

							# pos + exp
							T.fut.pos.exp <- data.pos.exp$futime
							T.fut.pos.exp[data.pos.exp$status==0] <- data.pos.exp$futime[data.pos.exp$status==0] + 
									     				     rgamma(sum(data.pos.exp$status==0), shape=1, rate=hazard.pos.exp[l])
							T.fut.pos.exp <- c(T.fut.pos.exp, rgamma(n.rem.pos.exp, shape=1, rate=hazard.pos.exp[l]))

							# pos + ctl
							T.fut.pos.ctl <- data.pos.ctl$futime
							T.fut.pos.ctl[data.pos.ctl$status==0] <- data.pos.ctl$futime[data.pos.ctl$status==0] + 
									     				     rgamma(sum(data.pos.ctl$status==0), shape=1, rate=hazard.pos.ctl[l])
							T.fut.pos.ctl <- c(T.fut.pos.ctl, rgamma(n.rem.pos.ctl, shape=1, rate=hazard.pos.ctl[l]))

							# neg + exp
							T.fut.neg.exp <- data.neg.exp$futime
							T.fut.neg.exp[data.neg.exp$status==0] <- data.neg.exp$futime[data.neg.exp$status==0] + 
									     				     rgamma(sum(data.neg.exp$status==0), shape=1, rate=hazard.neg.exp[l])
							T.fut.neg.exp <- c(T.fut.neg.exp, rgamma(n.rem.neg.exp, shape=1, rate=hazard.neg.exp[l]))

							# neg + ctl
							T.fut.neg.ctl <- data.neg.ctl$futime
							T.fut.neg.ctl[data.neg.ctl$status==0] <- data.neg.ctl$futime[data.neg.ctl$status==0] + 
									     				     rgamma(sum(data.neg.ctl$status==0), shape=1, rate=hazard.neg.ctl[l])
							T.fut.neg.ctl <- c(T.fut.neg.ctl, rgamma(n.rem.neg.ctl, shape=1, rate=hazard.neg.ctl[l]))


							# Combine current and future data
							enter.fut.pos.exp <- c(data.pos.exp$enter, enter.rem.pos.exp)
							T.date.fut.pos.exp <- T.fut.pos.exp + enter.fut.pos.exp

							enter.fut.pos.ctl <- c(data.pos.ctl$enter, enter.rem.pos.ctl)
							T.date.fut.pos.ctl <- T.fut.pos.ctl + enter.fut.pos.ctl

							enter.fut.neg.exp <- c(data.neg.exp$enter, enter.rem.neg.exp)
							T.date.fut.neg.exp <- T.fut.neg.exp + enter.fut.neg.exp

							enter.fut.neg.ctl <- c(data.neg.ctl$enter, enter.rem.neg.ctl)
							T.date.fut.neg.ctl <- T.fut.neg.ctl + enter.fut.neg.ctl
 

							# Calculate the time of final analysis for this futility analysis loop
							final <- sort(c(T.date.fut.pos.exp, T.date.fut.pos.ctl, T.date.fut.neg.exp, T.date.fut.neg.ctl))[n_events]


							# Censor data at the final analysis time
							futime.fut.pos.exp <- pmin(T.date.fut.pos.exp,final)-enter.fut.pos.exp
							status.fut.pos.exp <- T.date.fut.pos.exp <= final

							futime.fut.pos.ctl <- pmin(T.date.fut.pos.ctl,final)-enter.fut.pos.ctl
							status.fut.pos.ctl <- T.date.fut.pos.ctl <= final

							futime.fut.neg.exp <- pmin(T.date.fut.neg.exp,final)-enter.fut.neg.exp
							status.fut.neg.exp <- T.date.fut.neg.exp <= final

							futime.fut.neg.ctl <- pmin(T.date.fut.neg.ctl,final)-enter.fut.neg.ctl
							status.fut.neg.ctl <- T.date.fut.neg.ctl <= final


							# Update the posterior probability of efficacy
							hazard.fut.pos.exp <- rgamma(100000, shape=1e-2+sum(status.fut.pos.exp), rate=1+sum(futime.fut.pos.exp)) 
							hazard.fut.pos.ctl <- rgamma(100000, shape=1e-2+sum(status.fut.pos.ctl), rate=1+sum(futime.fut.pos.ctl))
							hazard.fut.neg.exp <- rgamma(100000, shape=1e-2+sum(status.fut.neg.exp), rate=1+sum(futime.fut.neg.exp)) 
							hazard.fut.neg.ctl <- rgamma(100000, shape=1e-2+sum(status.fut.neg.ctl), rate=1+sum(futime.fut.neg.ctl)) 

							# Determine if the trial is a success for the positive cohort
							if(mean(hazard.fut.pos.exp < hazard.fut.pos.ctl) > p_eff_grp[K]) {
								win.pos <- win.pos + 1
							}

							# Determine if the trial is a success for the negative cohort
							if(mean(hazard.fut.neg.exp < hazard.fut.neg.ctl) > p_eff_grp[K]) {
								win.neg <- win.neg + 1
							}

						}

						
						# stop due to futility
						if (stop.eff.pos[j]==0) {	
							stop.pos <- win.pos/2000 < p_fut_grp[j]
							stop.fut.pos[j] <- as.numeric(win.pos/2000 < p_fut_grp[j])

							if(stop.fut.pos[j]==0) {
								continue.pos[j] <- 1
							}
						}

						if(stop.eff.neg[j]==0) {
							stop.neg <- win.neg/2000 < p_fut_grp[j]
							stop.fut.neg[j] <- as.numeric(win.neg/2000 < p_fut_grp[j])

							if(stop.fut.neg[j]==0) {
								continue.neg[j] <- 1
							}
						}
					}


					##### Determine if to stop (any or both) marker cohorts at this check
					stop.any <- (stop.pos | stop.neg)

					if(stop.any) {

						marker.prevalence <- mean(marker==1)
						marker.status <- marker

						if(stop.pos & stop.neg) {
							sampsize <- nrow(data.j)
							sampsize.pos <- samp.pos
							sampsize.neg <- samp.neg							
							break
						} 
					}

				}
			}


			##### If the estimated marker prevalence is outside (m.prev.lwr, 1-m.prev.lwr), not to dichotomize the marker #####
			else {

				# subset the data available at this time based on arm
				data.exp <- data.j[data.j$trt==1, c("enter","futime","status")]
				data.ctl <- data.j[data.j$trt==0, c("enter","futime","status")]

 
				# simulate hazard rates for each arm in the overall cohort
				hazard.exp <- rgamma(100000, shape=1e-2+sum(data.exp$status), rate=1+sum(data.exp$futime) )
				hazard.ctl <- rgamma(100000, shape=1e-2+sum(data.ctl$status), rate=1+sum(data.ctl$futime) )


				###### check efficacy for the overall cohort ######
				eff.all <- mean(hazard.exp < hazard.ctl)		


            		###### If this is the final check, skip the futility analysis ######
				if(j==K) {
					stop.eff.all[j] <- as.numeric(eff.all >= p_eff_all[j])
					sampsize <- nrow(data.j)
					sampsize.pos <- NA
					sampsize.neg <- NA	
					marker.prevalence <- NA
					marker.status <- rep(NA, nrow(surdata))
				}

				##### If this is not the final check, perform efficacy and futility analysis #####	
				else {
					if(stop_eff & (eff.all >= p_eff_all[j])) {				# stop the overall cohort due to efficacy
						stop.any <- TRUE 
						stop.eff.all[j] <- 1
						sampsize <- nrow(data.j)
						sampsize.pos <- NA
						sampsize.neg <- NA	
						marker.prevalence <- NA
						marker.status <- rep(NA, nrow(surdata))							
						break
					}
		
					else {										# if not, check futility for the overall cohort
						enter.rem <- enter[-index]
						trt.rem <- trt[-index]

						enter.rem.exp <- enter.rem[trt.rem==1]
						enter.rem.ctl <- enter.rem[trt.rem==0]
						n.rem.exp <- sum(trt.rem==1)
						n.rem.ctl <- sum(trt.rem==0)
		
						win <- 0

						# futility analysis loops	
						for(l in 1:2000) {

							# Generate new time to event data
							T.fut.exp <- data.exp$futime
							T.fut.exp[data.exp$status==0] <- data.exp$futime[data.exp$status==0] + rgamma(sum(data.exp$status==0), shape=1, rate=hazard.exp[l])
							T.fut.exp <- c(T.fut.exp, rgamma(n.rem.exp, shape=1, rate=hazard.exp[l]))

							T.fut.ctl <- data.ctl$futime
							T.fut.ctl[data.ctl$status==0] <- data.ctl$futime[data.ctl$status==0] + rgamma(sum(data.ctl$status==0), shape=1, rate=hazard.ctl[l])
							T.fut.ctl <- c(T.fut.ctl, rgamma(n.rem.ctl, shape=1, rate=hazard.ctl[l]))

							# Combine current and future data
							enter.fut.exp <- c(data.exp$enter, enter.rem.exp)
							T.date.fut.exp <- T.fut.exp + enter.fut.exp

							enter.fut.ctl <- c(data.ctl$enter, enter.rem.ctl)
							T.date.fut.ctl <- T.fut.ctl + enter.fut.ctl

							# Censor data at the final analysis time
							final <- sort(c(T.date.fut.exp, T.date.fut.ctl))[n_events]

							futime.fut.exp <- pmin(T.date.fut.exp,final)-enter.fut.exp
							status.fut.exp <- T.date.fut.exp <= final

							futime.fut.ctl <- pmin(T.date.fut.ctl,final)-enter.fut.ctl
							status.fut.ctl <- T.date.fut.ctl <= final


							# Update the posterior probability of efficacy
							hazard.fut.exp <- rgamma(100000, shape=1e-2+sum(status.fut.exp), rate=1+sum(futime.fut.exp)) 
							hazard.fut.ctl <- rgamma(100000, shape=1e-2+sum(status.fut.ctl), rate=1+sum(futime.fut.ctl)) 

							# Determine if the trial is a success 
							if(mean(hazard.fut.exp < hazard.fut.ctl) > p_eff_all[K]) {
								win <- win + 1
							}
						}
		
 						if(win/2000 < p_fut_all[j]) {				# stop the overall cohort due to futility
							stop.any <- TRUE 
							stop.fut.all[j] <- 1
							sampsize <- nrow(data.j)
							sampsize.pos <- NA
							sampsize.neg <- NA
							marker.prevalence <- NA
							marker.status <- rep(NA, nrow(surdata))								
							break
						}

						else {
							continue.all[j] <- 1
						}
					}
				}
			}
		}
	    }		


	    else if (stop.pos & (!stop.neg)) {							# only the marker negative cohort continues

		# extract the marker negative cohort based on the dichotomization rule
		trt.neg <- trt[marker==0] 
		enter.neg <- enter[marker==0]
		T.date.neg <- T.date[marker==0]


         	#### Censor data appropriately #### 
         
         	date <- sort(T.date.neg)[ceiling(n_events*mean(marker==0)*check[j])]

         
         	index <- which(enter.neg < date)	       		# only consider those accrued by current date
         	enter.j <- enter.neg[index]			  		# accrual dates of those accured by current date
         	n.j <- length(index)			            	# sample size at current date
         	trt.j <- trt.neg[index] 		       		# trt group of those accrued by current date
 
         
        	#### Time to Event Data T ####
         	T.date.j <- T.date.neg[index]           	     		# event dates of those accrued by current date
         	T.j <- T.date.j - enter.j		      		# event times of those accrued by current date
         

	   	#### Initialize time and censoring vectors ####
         	T.cen <- rep(1,n.j)						# T.cen=1 for an observed event; =0 for right censoring.
         	censorLimit.T <- rep(NA,n.j)					# censorLimit.T is the time to event or right censoring (whichever happens first).

         
         	for(k in 1:n.j){  						# Loop over individuals on study at this interim time point (n.j set)
            
            	#  at interim j for those accrued
            	if(T.date.j[k] <= date) {
            		# T.cen[k] remains 1
            		censorLimit.T[k] <- T.j[k]
            	}

			else {
            		T.cen[k] <- 0
            		censorLimit.T[k] <- date - enter.j[k]
            	}
         	}
                  
	   
         	data.j = data.frame(enter=enter.j, futime=censorLimit.T, status=T.cen, trt=trt.j)

		samp.neg <- nrow(data.j)


		# subset the data available at this time based on arm
		data.neg.exp <- data.j[data.j$trt==1, c("enter","futime","status")]
		data.neg.ctl <- data.j[data.j$trt==0, c("enter","futime","status")]

		# simulate hazard rates for each arm in the marker negative cohort
		hazard.neg.exp <- rgamma(100000, shape=1e-2+sum(data.neg.exp$status), rate=1+sum(data.neg.exp$futime) )
		hazard.neg.ctl <- rgamma(100000, shape=1e-2+sum(data.neg.ctl$status), rate=1+sum(data.neg.ctl$futime) )


		###### check efficacy for the marker negative cohort ######
		eff.neg <- mean(hazard.neg.exp < hazard.neg.ctl)


           	###### If this is the final check, skip the futility analysis ######
		if(j==K) {
			stop.eff.neg[j] <- as.numeric(eff.neg >= p_eff_grp[j])
			sampsize <- samp.pos + samp.neg
			sampsize.pos <- samp.pos
			sampsize.neg <- samp.neg
		}

		##### If this is not the final check, perform efficacy and futility analysis #####	
		else {
			if(stop_eff & (eff.neg >= p_eff_grp[j])) {					# stop the marker negative cohort due to efficacy
				stop.neg <- TRUE 
				stop.eff.neg[j] <- 1
				sampsize <- samp.pos + samp.neg
				sampsize.pos <- samp.pos
				sampsize.neg <- samp.neg
				break
			}
		
			else {											# if not, check futility for the marker negative cohort
				enter.rem <- enter.neg[-index]
				trt.rem <- trt.neg[-index]

				enter.rem.neg.exp <- enter.rem[trt.rem==1]
				enter.rem.neg.ctl <- enter.rem[trt.rem==0]
				n.rem.neg.exp <- sum(trt.rem==1)
				n.rem.neg.ctl <- sum(trt.rem==0)
		
				win <- 0	

				# futility analysis loops
				for(l in 1:2000) {

					# Generate new time to event data
					T.fut.neg.exp <- data.neg.exp$futime
					T.fut.neg.exp[data.neg.exp$status==0] <- data.neg.exp$futime[data.neg.exp$status==0] + 
									     		     rgamma(sum(data.neg.exp$status==0), shape=1, rate=hazard.neg.exp[l])
					T.fut.neg.exp <- c(T.fut.neg.exp, rgamma(n.rem.neg.exp, shape=1, rate=hazard.neg.exp[l]))

					T.fut.neg.ctl <- data.neg.ctl$futime
					T.fut.neg.ctl[data.neg.ctl$status==0] <- data.neg.ctl$futime[data.neg.ctl$status==0] + 
									     		     rgamma(sum(data.neg.ctl$status==0), shape=1, rate=hazard.neg.ctl[l])
					T.fut.neg.ctl <- c(T.fut.neg.ctl, rgamma(n.rem.neg.ctl, shape=1, rate=hazard.neg.ctl[l]))


					# Combine current and future data
					enter.fut.neg.exp <- c(data.neg.exp$enter, enter.rem.neg.exp)
					T.date.fut.neg.exp <- T.fut.neg.exp + enter.fut.neg.exp

					enter.fut.neg.ctl <- c(data.neg.ctl$enter, enter.rem.neg.ctl)
					T.date.fut.neg.ctl <- T.fut.neg.ctl + enter.fut.neg.ctl


					# Censor data at the final analysis time
					final <- sort(c(T.date.fut.neg.exp, T.date.fut.neg.ctl))[ceiling(n_events*mean(marker==0))]

					futime.fut.neg.exp <- pmin(T.date.fut.neg.exp,final)-enter.fut.neg.exp
					status.fut.neg.exp <- T.date.fut.neg.exp <= final
					futime.fut.neg.ctl <- pmin(T.date.fut.neg.ctl,final)-enter.fut.neg.ctl
					status.fut.neg.ctl <- T.date.fut.neg.ctl <= final


					# Update the posterior probability of efficacy
					hazard.fut.neg.exp <- rgamma(100000, shape=1e-2+sum(status.fut.neg.exp), rate=1+sum(futime.fut.neg.exp)) 
					hazard.fut.neg.ctl <- rgamma(100000, shape=1e-2+sum(status.fut.neg.ctl), rate=1+sum(futime.fut.neg.ctl)) 


					# Determine if the trial is a success 
					if(mean(hazard.fut.neg.exp < hazard.fut.neg.ctl) > p_eff_grp[K]) {
						win <- win + 1
					}
				}
		
 				if(win/2000 < p_fut_grp[j]) {							# stop the marker negative cohort due to futility
					stop.neg <- TRUE
					stop.fut.neg[j] <- 1
					sampsize <- samp.pos + samp.neg
					sampsize.pos <- samp.pos
					sampsize.neg <- samp.neg
					break
				}

				else {
					continue.neg[j] <- 1
				}					
			}
		}
	    }



	    else if (stop.neg & (!stop.pos)) {								# only the marker positive cohort continues

		# extract the marker positive cohort based on the dichotomization rule
		trt.pos <- trt[marker==1] 
		enter.pos <- enter[marker==1]
		T.date.pos <- T.date[marker==1]


         	#### Censor data appropriately #### 
         
         	date <- sort(T.date.pos)[ceiling(n_events*mean(marker==1)*check[j])]

         
         	index <- which(enter.pos < date)	       		# only consider those accrued by current date
         	enter.j <- enter.pos[index]			  		# accrual dates of those accured by current date
         	n.j <- length(index)			            	# sample size at current date
         	trt.j <- trt.pos[index] 		       		# trt group of those accrued by current date
 
         
        	#### Time to Event Data T ####
         	T.date.j <- T.date.pos[index]           	     		# event dates of those accrued by current date
         	T.j <- T.date.j - enter.j		      		# event times of those accrued by current date
         

	   	#### Initialize time and censoring vectors ####
         	T.cen <- rep(1,n.j)						# T.cen=1 for an observed event; =0 for right censoring.
         	censorLimit.T <- rep(NA,n.j)					# censorLimit.T is the time to event or right censoring (whichever happens first).

         
         	for(k in 1:n.j){  						# Loop over individuals on study at this interim time point (n.j set)
            
            	#  at interim j for those accrued
            	if(T.date.j[k] <= date) {
            		# T.cen[k] remains 1
            		censorLimit.T[k] <- T.j[k]
            	}

			else {
            		T.cen[k] <- 0
            		censorLimit.T[k] <- date - enter.j[k]
            	}
         	}
                  
	   
         	data.j = data.frame(enter=enter.j, futime=censorLimit.T, status=T.cen, trt=trt.j)

		samp.pos <- nrow(data.j)


		# subset the data available at this time based on arm
		data.pos.exp <- data.j[data.j$trt==1, c("enter","futime","status")]
		data.pos.ctl <- data.j[data.j$trt==0, c("enter","futime","status")]

		# simulate hazard rates for each arm in the marker positive cohort
		hazard.pos.exp <- rgamma(100000, shape=1e-2+sum(data.pos.exp$status), rate=1+sum(data.pos.exp$futime) )
		hazard.pos.ctl <- rgamma(100000, shape=1e-2+sum(data.pos.ctl$status), rate=1+sum(data.pos.ctl$futime) )


		###### check efficacy for the marker positive cohort ######
		eff.pos <- mean(hazard.pos.exp < hazard.pos.ctl)	


            ###### If this is the final check, skip the futility analysis ######
		if(j==K) {
			stop.eff.pos[j] <- as.numeric(eff.pos >= p_eff_grp[j])
			sampsize <- samp.pos + samp.neg
			sampsize.pos <- samp.pos
			sampsize.neg <- samp.neg
		}

		##### If this is not the final check, perform efficacy and futility analysis #####	
		else {
			if(stop_eff & (eff.pos >= p_eff_grp[j])) {					# stop the marker positive cohort due to efficacy
				stop.pos <- TRUE 
				stop.eff.pos[j] <- 1
				sampsize <- samp.pos + samp.neg
				sampsize.pos <- samp.pos
				sampsize.neg <- samp.neg
				break
			}
		
			else {											# if not, check futility for the marker positive cohort
				enter.rem <- enter.pos[-index]
				trt.rem <- trt.pos[-index]

				enter.rem.pos.exp <- enter.rem[trt.rem==1]
				enter.rem.pos.ctl <- enter.rem[trt.rem==0]
				n.rem.pos.exp <- sum(trt.rem==1)
				n.rem.pos.ctl <- sum(trt.rem==0)
		
				win <- 0	

				# futility analysis loops
				for(l in 1:2000) {

					# Generate new time to event data
					T.fut.pos.exp <- data.pos.exp$futime
					T.fut.pos.exp[data.pos.exp$status==0] <- data.pos.exp$futime[data.pos.exp$status==0] + 
									     		     rgamma(sum(data.pos.exp$status==0), shape=1, rate=hazard.pos.exp[l])
					T.fut.pos.exp <- c(T.fut.pos.exp, rgamma(n.rem.pos.exp, shape=1, rate=hazard.pos.exp[l]))

					T.fut.pos.ctl <- data.pos.ctl$futime
					T.fut.pos.ctl[data.pos.ctl$status==0] <- data.pos.ctl$futime[data.pos.ctl$status==0] + 
									     		     rgamma(sum(data.pos.ctl$status==0), shape=1, rate=hazard.pos.ctl[l])
					T.fut.pos.ctl <- c(T.fut.pos.ctl, rgamma(n.rem.pos.ctl, shape=1, rate=hazard.pos.ctl[l]))


					# Combine current and future data
					enter.fut.pos.exp <- c(data.pos.exp$enter, enter.rem.pos.exp)
					T.date.fut.pos.exp <- T.fut.pos.exp + enter.fut.pos.exp

					enter.fut.pos.ctl <- c(data.pos.ctl$enter, enter.rem.pos.ctl)
					T.date.fut.pos.ctl <- T.fut.pos.ctl + enter.fut.pos.ctl


					# Censor data at the final analysis time
					final <- sort(c(T.date.fut.pos.exp, T.date.fut.pos.ctl))[ceiling(n_events*mean(marker==1))]

					futime.fut.pos.exp <- pmin(T.date.fut.pos.exp,final)-enter.fut.pos.exp
					status.fut.pos.exp <- T.date.fut.pos.exp <= final
					futime.fut.pos.ctl <- pmin(T.date.fut.pos.ctl,final)-enter.fut.pos.ctl
					status.fut.pos.ctl <- T.date.fut.pos.ctl <= final


					# Update the posterior probability of efficacy
					hazard.fut.pos.exp <- rgamma(100000, shape=1e-2+sum(status.fut.pos.exp), rate=1+sum(futime.fut.pos.exp)) 
					hazard.fut.pos.ctl <- rgamma(100000, shape=1e-2+sum(status.fut.pos.ctl), rate=1+sum(futime.fut.pos.ctl)) 

					# Determine if the trial is a success 
					if(mean(hazard.fut.pos.exp < hazard.fut.pos.ctl) > p_eff_grp[K]) {
						win <- win + 1
					}
				}
		
 				if(win/2000 < p_fut_grp[j]) {								# stop the marker positive cohort due to futility
					stop.pos <- TRUE
					stop.fut.pos[j] <- 1
					sampsize <- samp.pos + samp.neg
					sampsize.pos <- samp.pos
					sampsize.neg <- samp.neg
					break
				}

				else {
					continue.pos[j] <- 1
				}					
			}
		}
	    }

	   
         ######################################################                      

      }   #  end loop over interim checks
      

	### Determine the final result for this trial (overall efficacy, subgroup efficacy, or overall futility)
	if(is.na(marker.prevalence)) {										# the marker is not dichotomized
		if(sum(stop.eff.all)>0) {										# overall efficacy
			res.final <- "eff_all"
		}
		else {
			res.final <- "eff_no"										
		}
	}

	else {														# the marker is dichotomized
		if(sum(stop.eff.pos)>0 & sum(stop.eff.neg)>0) {							# efficacy in both marker subgroups
			res.final <- "eff_all"
		}
		else if (sum(stop.eff.pos)==0 & sum(stop.eff.neg)==0) {					# no efficacy in either marker subgroup
			res.final <- "eff_no"
		}		
		else if (sum(stop.eff.pos)> 0 & sum(stop.eff.neg)==0) {
			res.final <- "eff_pos"										# efficacy in one marker subgroup but not the other
		}
		else {
			res.final <- "eff_neg"
		}
	}		


	return(list(stop.eff.pos=stop.eff.pos, stop.eff.neg=stop.eff.neg, stop.eff.all=stop.eff.all, stop.fut.pos=stop.fut.pos, stop.fut.neg=stop.fut.neg, stop.fut.all=stop.fut.all,
                  continue.pos=continue.pos, continue.neg=continue.neg, continue.all=continue.all, res.final=res.final,
	            sampsize=sampsize, sampsize.pos=sampsize.pos, sampsize.neg=sampsize.neg, 
			marker.prevalence=marker.prevalence, x=x, marker.status=marker.status,
			fail.MCMC=FALSE))   

} #end function



