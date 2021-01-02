###################################################################
# Code to do frequentist adaptive design in the setting of two arms, one experimental arm and one control arm #
# Time-to-event endpoint T with right censoring allowed  # 
# Assume one continuous biomarker x is available for each patient #
# The treatment-by-marker interaction effect is used to decide whether to dichotomize the continuous marker, and the minimum p-value cutoff is specified by the user #
# One-sided log-rank test is used to test for efficacy, and the p-value cutoff is specified by the user #
# Point estimate of HR is used to test for futility, and the threshold is specified by the user #
# If decide to dichotomize, record the marker status for each patient #
###################################################################
         
            
Freq_adaptive_design <- function(surdata,
				#####################################################################################################	 	
				### surdata = a data frame with the observations arranged by row, and including the columns below (in the following order) ###
				### treatment arm indicator trt, entry dates enter, continuous marker values x, event times T, and observed TTE dates T_date. ###  
				##########		     							   		   
				  	   n_events = 400,						         # the number of events expected to be observed at the end of the trial
				  	   m.prev.lwr = 0.25,						   	   # the lower bound of marker prevalence allowed for the grid search of cutpoint
                          	   check = c(0.25, 0.50, 0.75, 1.00), 			   	   # interim analysis times based on % of observed events for T
				  	   m.increment = 0.01, 							   # the increment of the sequence of marker cutpoint searched 
				  	   p_int = rep(0.01,4),							   # a vector of minimum p-value cutoffs for the treatment-by-marker interaction to decide whether to dichotomize (length=number of checks)
				  	   p_eff_all = rep(0.01,4),  						   # a vector of overall cohort one-sided log-rank test p-value cutoffs to claim efficacy (length=number of checks)
					   p_eff_grp = rep(0.1, 4),						   # a vector of marker subgroup one-sided log-rank test p-value cutoffs to claim efficacy (length=number of checks)
				  	   HR_fut_all = rep(1,3),						   # a vector of overall cohort thresholds of HR above which to claim futility (length=number of checks - 1)
					   HR_fut_grp = rep(1,3),						   # a vector of marker subgroup thresholds of HR above which to claim futility (length=number of checks - 1)
				  	   stop_eff = TRUE) { 							   # a logical indictor of whether to allow early stopping due to efficacy		 	                                              		
		
 
      
   	### Code here to load in the functions that will be called later
   	source("Freq_detect_marker.R")  


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

       	date <- sort(T.date)[ceiling(n_events*check[j])]


         	index <- which(enter < date)	       		# only consider those accrued by current date
         	enter.j <- enter[index]			  		# accrual dates of those accured by current date
         	n.j <- length(index)			        	# sample size at current date
         	trt.j <- trt[index] 		       		# trt group of those accrued by current date
         	x.j <- x[index]                      	    	# continuous marker value of those accrued by current date
 
         
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


		##### Call the function to decide if and how to dichotomize the marker based on data available at this time
		result.j = Freq_detect_marker(data=data.j[,c("futime","status","x","trt")], x.all=x, m.prev.lwr = m.prev.lwr, m.increment = m.increment, p_int = p_int[j])


		##### If dichotomize the marker #####
		if( result.j$dicho ) {

			# dichotomize the continuous marker for individuals enrolled at this time
			data.j$marker <- result.j$marker

			samp.pos <- sum(data.j$marker==1)
			samp.neg <- sum(data.j$marker==0)


			# subset the data available at this time based on marker positivity and arm
			data.pos <- data.j[data.j$marker==1, c("enter","futime","status","trt")]
			data.neg <- data.j[data.j$marker==0, c("enter","futime","status","trt")]


			###### do log-rank test and fit proportional hazard model within each marker cohort ######
			test.pos <- survdiff(Surv(futime, status) ~ trt, data=data.pos)			# the marker positive cohort
			idx.pos <- which(names(test.pos$n) == "trt=1")
			z.pos <- (test.pos$obs[idx.pos]-test.pos$exp[idx.pos])/sqrt(test.pos$var[idx.pos,idx.pos])
			pvalue.pos <- pnorm(z.pos, lower.tail = TRUE)
			ph.pos <- coxph(Surv(futime, status) ~ trt, data=data.pos)


			test.neg <- survdiff(Surv(futime, status) ~ trt, data=data.neg)			# the marker negative cohort
			idx.neg <- which(names(test.neg$n) == "trt=1")
			z.neg <- (test.neg$obs[idx.neg]-test.neg$exp[idx.neg])/sqrt(test.neg$var[idx.neg,idx.neg])
			pvalue.neg <- pnorm(z.neg, lower.tail = TRUE)
			ph.neg <- coxph(Surv(futime, status) ~ trt, data=data.neg)


            	###### If this is the final check, skip the futility analysis ######
			if(j==K) {
				stop.eff.pos[j] <- as.numeric(pvalue.pos < p_eff_grp[j])
				stop.eff.neg[j] <- as.numeric(pvalue.neg < p_eff_grp[j])
				sampsize <- nrow(data.j)
				sampsize.pos <- samp.pos
				sampsize.neg <- samp.neg
				marker.prevalence <- mean(result.j$marker.all==1)
				marker.status <- result.j$marker.all
			}

			##### If this is not the final check, perform efficacy and futility analysis #####	
			else {

				if(stop_eff & (pvalue.pos < p_eff_grp[j])) {					# stop the marker positive cohort due to efficacy
					stop.pos <- TRUE 
					stop.eff.pos[j] <- 1
				}

				else {											# if the marker positive cohort fails to stop early for efficacy, test futility
					stop.pos <- coef(ph.pos) > log(HR_fut_grp[j])
					stop.fut.pos[j] <- as.numeric( coef(ph.pos) > log(HR_fut_grp[j]) )

					if(stop.fut.pos[j]==0) {
						continue.pos[j] <- 1
					}
				}


				if(stop_eff & (pvalue.neg < p_eff_grp[j])) {					# stop the marker negative cohort due to efficacy
					stop.neg <- TRUE 
					stop.eff.neg[j] <- 1
				}

				else {											# if the marker negative cohort fails to stop early for efficacy, test futility
					stop.neg <- coef(ph.neg) > log(HR_fut_grp[j])
					stop.fut.neg[j] <- as.numeric( coef(ph.neg) > log(HR_fut_grp[j]) )

					if(stop.fut.neg[j]==0) {
						continue.neg[j] <- 1
					}
				}


				##### Determine if to stop (any or both) marker cohorts at this check
				stop.any <- (stop.pos | stop.neg)

				if(stop.any) {

					marker <- result.j$marker.all
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


		##### If not to dichotomize the marker #####
		else {

			###### do log-rank test and fit proportional hazard model in the overall cohort ######
			test.all <- survdiff(Surv(futime, status) ~ trt, data=data.j)		
			idx.all <- which(names(test.all$n) == "trt=1")
			z.all <- (test.all$obs[idx.all]-test.all$exp[idx.all])/sqrt(test.all$var[idx.all,idx.all])
			pvalue.all <- pnorm(z.all, lower.tail = TRUE)
			ph.all <- coxph(Surv(futime, status) ~ trt, data=data.j)
		

            	###### If this is the final check, skip the futility analysis ######
			if(j==K) {
				stop.eff.all[j] <- as.numeric(pvalue.all < p_eff_all[j])
				sampsize <- nrow(data.j)
				sampsize.pos <- NA
				sampsize.neg <- NA	
				marker.prevalence <- NA
				marker.status <- rep(NA, nrow(surdata))
			}

			##### If this is not the final check, perform efficacy and futility analysis #####	
			else {
				if(stop_eff & (pvalue.all < p_eff_all[j])) {					# stop the overall cohort due to efficacy
					stop.any <- TRUE 
					stop.eff.all[j] <- 1
					sampsize <- nrow(data.j)
					sampsize.pos <- NA
					sampsize.neg <- NA	
					marker.prevalence <- NA	
					marker.status <- rep(NA, nrow(surdata))						
					break
				}
		
				else if (coef(ph.all) > log(HR_fut_all[j])) {					# if not, test futility for the overall cohort
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


	    else if (stop.pos & (!stop.neg)) {							# only the marker negative cohort continues

		# extract the marker negative cohort based on the dichotomization rule
		trt.neg <- trt[marker==0] 
		enter.neg <- enter[marker==0]
		T.date.neg <- T.date[marker==0]


         	#### Censor data appropriately #### 
         
         	date <- sort(T.date.neg)[ceiling(n_events*mean(marker==0)*check[j])]
		
         
         	index <- which(enter.neg < date)	       					# only consider those accrued by current date
         	enter.j <- enter.neg[index]			  					# accrual dates of those accured by current date
         	n.j <- length(index)			            				# sample size at current date
         	trt.j <- trt.neg[index] 		       					# trt group of those accrued by current date
 
         
        	#### Time to Event Data T ####
         	T.date.j <- T.date.neg[index]           	     					# event dates of those accrued by current date
         	T.j <- T.date.j - enter.j		      					# event times of those accrued by current date
         

	   	#### Initialize time and censoring vectors ####
         	T.cen <- rep(1,n.j)									# T.cen=1 for an observed event; =0 for right censoring.
         	censorLimit.T <- rep(NA,n.j)								# censorLimit.T is the time to event or right censoring (whichever happens first).

         
         	for(k in 1:n.j){  									# Loop over individuals on study at this interim time point (n.j set)
            
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


		###### Do log-rank test and fit proportional hazard model in the marker-negative cohort ######
		test.neg <- survdiff(Surv(futime, status) ~ trt, data=data.j)		
		idx.neg <- which(names(test.neg$n) == "trt=1")
		z.neg <- (test.neg$obs[idx.neg]-test.neg$exp[idx.neg])/sqrt(test.neg$var[idx.neg,idx.neg])
		pvalue.neg <- pnorm(z.neg, lower.tail = TRUE)
		ph.neg <- coxph(Surv(futime, status) ~ trt, data=data.j)


           	###### If this is the final check, skip the futility analysis ######
		if(j==K) {
			stop.eff.neg[j] <- as.numeric(pvalue.neg < p_eff_grp[j])
			sampsize <- samp.pos + samp.neg
			sampsize.pos <- samp.pos
			sampsize.neg <- samp.neg
		}

		##### If this is not the final check, perform efficacy and futility analysis #####	
		else {
			if(stop_eff & (pvalue.neg < p_eff_grp[j])) {				# stop the marker negative cohort due to efficacy
				stop.neg <- TRUE 
				stop.eff.neg[j] <- 1
				sampsize <- samp.pos + samp.neg
				sampsize.pos <- samp.pos
				sampsize.neg <- samp.neg
				break
			}
		
			else if (coef(ph.neg) > log(HR_fut_grp[j])) {				# if not, test futility for the marker negative cohort
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



	    else if (stop.neg & (!stop.pos)) {							# only the marker positive cohort continues

		# extract the marker positive cohort based on the dichotomization rule
		trt.pos <- trt[marker==1] 
		enter.pos <- enter[marker==1]
		T.date.pos <- T.date[marker==1]


         	#### Censor data appropriately #### 
         
         	date <- sort(T.date.pos)[ceiling(n_events*mean(marker==1)*check[j])]

         
         	index <- which(enter.pos < date)	       					# only consider those accrued by current date
         	enter.j <- enter.pos[index]			  					# accrual dates of those accured by current date
         	n.j <- length(index)			            				# sample size at current date
         	trt.j <- trt.pos[index] 		       					# trt group of those accrued by current date
 
         
        	#### Time to Event Data T ####
         	T.date.j <- T.date.pos[index]           	     					# event dates of those accrued by current date
         	T.j <- T.date.j - enter.j		      					# event times of those accrued by current date
         

	   	#### Initialize time and censoring vectors ####
         	T.cen <- rep(1,n.j)									# T.cen=1 for an observed event; =0 for right censoring.
         	censorLimit.T <- rep(NA,n.j)								# censorLimit.T is the time to event or right censoring (whichever happens first).

         
         	for(k in 1:n.j){  									# Loop over individuals on study at this interim time point (n.j set)
            
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


		###### Do log-rank test and fit proportional hazard model in the marker-positive cohort ######
		test.pos <- survdiff(Surv(futime, status) ~ trt, data=data.j)		
		idx.pos <- which(names(test.pos$n) == "trt=1")
		z.pos <- (test.pos$obs[idx.pos]-test.pos$exp[idx.pos])/sqrt(test.pos$var[idx.pos,idx.pos])
		pvalue.pos <- pnorm(z.pos, lower.tail = TRUE)
		ph.pos <- coxph(Surv(futime, status) ~ trt, data=data.j)


            ###### If this is the final check, skip the futility analysis ######
		if(j==K) {
			stop.eff.pos[j] <- as.numeric(pvalue.pos < p_eff_grp[j])
			sampsize <- samp.pos + samp.neg
			sampsize.pos <- samp.pos
			sampsize.neg <- samp.neg
		}

		##### If this is not the final check, perform efficacy and futility analysis #####	
		else {
			if(stop_eff & (pvalue.pos < p_eff_grp[j])) {				# stop the marker positive cohort due to efficacy
				stop.pos <- TRUE 
				stop.eff.pos[j] <- 1
				sampsize <- samp.pos + samp.neg
				sampsize.pos <- samp.pos
				sampsize.neg <- samp.neg
				break
			}
		
			else if (coef(ph.pos) > log(HR_fut_grp[j])) {				# if not, test futility for the marker negative cohort
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
			marker.prevalence=marker.prevalence, x=x, marker.status=marker.status))   

} #end function



