##########################################################################################
#### Simulate clinical trial data including arm assignment, continuous marker, follow-up time and right censoring indicator.	####
#### Two arms (one control group and one experimental group) are assumed.	####
#### Exponential distribution is assumed for T. The median of T for baseline hazard should be specified.	####
#### The continuous marker x is assumed to affect T under a proportional hazards model.	####
#### The mapping functions of the log hazard for T, as a function of marker x, need to be specified for each arm.	 ####
#### This function is called by Sim_Shell_Exponential.	####
##########################################################################################


simdata <- function(n = 300, 								# total trial size (both arms combined)
                    last.on.days = 365*2,      				                # date that the last patient enters (length of accural period in days)                   
                    md.t = 365,	                                                   	# median of endpoint T (exponential dist) for baseline hazard of the control group 
                    percent.exp = 0.5,                                                	# randomization ratio to the experimental group
                    marker.specs = list(dist="beta", param=c(1,1),range=c(0,1), ref="median"),        	# a list giving the distribution, range, and referece point of the continuous marker x                                     
                    FUN.ctl = function(x) {0*x},                                       	# mapping function of the log hazard as a function of marker x for the control group
                    FUN.exp = function(x) {0*x} )                                      	# mapping function of the log hazard as a function of marker x for the experimental group
{
 

	# Compute number of patients in each arm
	nsize.exp <- round(n*percent.exp)
	nsize.ctl <- n-nsize.exp


	# Set up treatment arm indicators
	trt <- c(rep(0, nsize.ctl), rep(1, nsize.exp) )					# Treatment assignment 


	##### Entry Dates, Evenly Spaced within each arm  #####
	enter <- rep(NA,length=n)
	# 
	# # NOT RANDOM entry times, just evenly spaced times across indices within each arm and adding small bit of noise
	# 
	enter[trt==0] <- sample(jitter(seq(0,last.on.days,length=sum(trt==0))),sum(trt==0),replace=FALSE)
	enter[trt==1] <- sample(jitter(seq(0,last.on.days,length=sum(trt==1))),sum(trt==1),replace=FALSE)
	# 
	enter[enter<0] <- 0  # change negative dates to day zero (zero is first day of trial)
	


	##### Generate continuous marker x and Exponential endpoint T for each patient #####  

	if(marker.specs$dist=="beta")
	{
		# generate continuous marker x for the control group
		x.ctl <- rbeta(nsize.ctl, shape1=marker.specs$param[1], shape2=marker.specs$param[2])
		x.ctl <- marker.specs$range[1] + marker.specs$range[2]*x.ctl

            	# generate continuous marker x for the experimental group
		x.exp <- rbeta(nsize.exp, shape1=marker.specs$param[1], shape2=marker.specs$param[2])
		x.exp <- marker.specs$range[1] + marker.specs$range[2]*x.exp
	}

	# More options of distributions for marker x can be added here later.


	if(marker.specs$ref=="median")
	{
		ref <- median(x.ctl)									# determine the reference point (in the control group)
	}

	# More options for the choice of reference point can be added here later.


	# Guarantee that a value of x.ctl is exactly equal to the reference value
	# This is simply a trick to ease the calculation of log hazard ratio (as a function of marker) relative to the reference value
	x.ctl[which.min(abs(x.ctl-ref))] <- ref


	# Calculate the effect in terms of the log hazard ratio for T, as a function of continuous marker x, relative to the reference point in the control group
	LHR.ctl <- unlist(lapply(x.ctl,FUN.ctl))-FUN.ctl(ref)
	LHR.exp <- unlist(lapply(x.exp,FUN.exp))-FUN.ctl(ref)


	# Generate endpoint T based on log relative hazard as a function of marker value
	T.ctl <- -md.t*log(runif(nsize.ctl))/(exp(LHR.ctl)*log(2))					# the control group
	T.exp <- -md.t*log(runif(nsize.exp))/(exp(LHR.exp)*log(2))					# the experimental group



	##### Return a list with the following two components:	#####
	### 1.) A data frame surdata with treatment arm indicator trt, entry dates enter, continuous marker values x, event times T, and observed TTE dates T_date.  
	###     Each row represents a patient, and the patients are ordered by trt status (the control group first).	###   
	### 2.) A numerical scalar ref which gives the reference value of marker x.		###

	return(list(surdata = data.frame(trt=trt, enter=enter,  x=c(x.ctl, x.exp), T=c(T.ctl, T.exp), T_date=c(T.ctl, T.exp)+enter), ref = ref) )

}
      

   
    
   