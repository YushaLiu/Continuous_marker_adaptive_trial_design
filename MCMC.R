##########################################################################################
#### Perform MCMC to generate posterior samples of penalized spline coefficients and regularization parameters
#### This function is called by Bayes_penalized_splines, not directly called by the user
##########################################################################################

MCMC <- function(data,formula,ncolz,prior,multiplier1,multiplier2,multiplier3,burnin,B,thin)
{

   library(MASS)
   library(survival)  
   library(Matrix)
   library(coda)


   ####################################### MCMC initialization ##############################
   # determine the initial values for the PH regression coefficients except for the intercept 
   initfit <- coxph(formula=formula,data=data)
   theta <- initfit$coef

   theta_g <- theta[1:(ncolz+1)]
   u_g <- theta_g[-c(1)]

   theta_h <- theta[(ncolz+2):length(theta)]
   u_h <- theta_h[-c(1:2)]
  
   K <- ncolz-2
   sigma2_g <- 1/rgamma(1,shape=(K+2)/2+prior$alpha,rate=sum(u_g^2)/2+prior$beta)
   sigma2_h <- 1/rgamma(1,shape=(K+2)/2+prior$alpha,rate=sum(u_h^2)/2+prior$beta)

 
   # determine the initial proposal variance for Metropolis random walk 
   # of the PH regression coefficients except for the intercept term
   data_g <- data[,3:(ncolz+3)]
   data_h <- data[,(ncolz+4):ncol(data)]

   temp_h <- as.numeric(as.matrix(data_h) %*% theta_h)
   temp_g <- as.numeric(as.matrix(data_g) %*% theta_g)

   data_g$h <- temp_h
   data_h$g <- temp_g

   formula1 <- as.formula(paste("Surv(futime,status) ~ ", paste(colnames(data_g), collapse="+") ))
   initfit_g <- coxph(formula=formula1, data=cbind(data[,1:2],data_g))
   var_g <- initfit_g$var[-c(ncol(data_g)),-c(ncol(data_g))]

   formula2 <- as.formula(paste("Surv(futime,status) ~ ", paste(colnames(data_h), collapse="+") ))
   initfit_h <- coxph(formula=formula2, data=cbind(data[,1:2],data_h))
   var_h <- (multiplier1^2)*initfit_h$var[-c(ncol(data_h)),-c(ncol(data_h))]


   # determine the initial value and random walk proposal variance for the baseline hazard lambda0,
   # which has a conjugate Gamma full conditional
   alpha_lambda0 <- sum(data$status)
   beta_lambda0 <- sum(exp(temp_g+temp_h)*data$futime)
   lambda0 <- (alpha_lambda0-1)/beta_lambda0
   var_lambda0 <- alpha_lambda0/(beta_lambda0^2)
   
   
   # put together the initial values and proposal variance for all coefficients,
   # including the intercept term, i.e., log(lambda0)
   theta_g <- c(log(lambda0),theta_g)
   var_g <- (multiplier1^2)*as.matrix(bdiag(var_lambda0*diag(1)/(lambda0^2), var_g))
   theta <- c(theta_g,theta_h)

   
   
   ########################################### MCMC starts #####################################
   L <- burnin+B*thin
   theta.all <- matrix(NA,nrow=L,ncol=length(theta))
   theta.res <- matrix(NA,nrow=B,ncol=length(theta))

   sigma2.all <- matrix(NA,nrow=L,ncol=2)
   sigma2.res <- matrix(NA,nrow=B,ncol=2)

   ll <- 0
   acpt_g <- 0
   acpt_h <- 0


   # indices corresponding to the control group (i.e., without interaction terms)
   idx_g <- 1:(ncolz+2)


   # calculate initial log likelihood (only the part relevant to the acpt prob calculation)
   data_g <- as.matrix(cbind(1,data[,3:(ncolz+3)]))
   data_h <- as.matrix(data[,(ncolz+4):ncol(data)])
   coef_g <- as.numeric(data_g %*% theta_g)
   coef_h <- as.numeric(data_h %*% theta_h)



   for(l in 1:L)
   {
      ##########  update theta_g  ##########
      newtheta_g <- mvrnorm(1,theta_g,var_g)
      newu_g <- newtheta_g[-c(1:2)]

      newtheta <- theta
      newtheta[idx_g] <- newtheta_g

      newcoef_g <- as.numeric(data_g %*% newtheta_g)
      newloglik <- sum((newcoef_g+coef_h)*data$status - exp(newcoef_g+coef_h)*data$futime)

	coef_g <- as.numeric(data_g %*% theta_g)
	loglik <- sum((coef_g+coef_h)*data$status - exp(coef_g+coef_h)*data$futime)
      
      log_acpt <- newloglik-0.5*sum(newu_g^2)/sigma2_g-loglik+0.5*sum(u_g^2)/sigma2_g

      if(log(runif(1))<log_acpt)
      {
         theta <- newtheta
         theta_g <- newtheta_g
	   coef_g <- newcoef_g
         u_g <- newu_g
         acpt_g <- acpt_g+1
      }


      ##########  update sigma2_g  ##########
      sigma2_g <- 1/rgamma(1,shape=(K+2)/2+prior$alpha,rate=sum(u_g^2)/2+prior$beta)


      ##########  update theta_h  ##########
      newtheta_h <- mvrnorm(1,theta_h,var_h)
      newu_h <- newtheta_h[-c(1:2)]

	newtheta <- theta
	newtheta[-idx_g] <- newtheta_h

      newcoef_h <- as.numeric(data_h %*% newtheta_h)
      newloglik <- sum((coef_g+newcoef_h)*data$status - exp(coef_g+newcoef_h)*data$futime) 

	coef_h <- as.numeric(data_h %*% theta_h)
	loglik <- sum((coef_g+coef_h)*data$status - exp(coef_g+coef_h)*data$futime)

      log_acpt <- newloglik-0.5*sum(newu_h^2)/sigma2_h-loglik+0.5*sum(u_h^2)/sigma2_h

      if(log(runif(1))<log_acpt)
      {
         theta <- newtheta
         theta_h <- newtheta_h
	   coef_h <- newcoef_h
         u_h <- newu_h
         acpt_h <- acpt_h+1
      }


      ##########  update sigma2_h  ##########
      sigma2_h <- 1/rgamma(1,shape=(K+2)/2+prior$alpha,rate=sum(u_h^2)/2+prior$beta)


      ##########  save the values of the new iteration  ##########
      theta.all[l,] <- theta
      sigma2.all[l,] <- c(sigma2_g,sigma2_h)
      
      if( (l>burnin) && ((l-burnin)%%thin==0) )
      {
         ll <- ll+1
         theta.res[ll,] <- theta
         sigma2.res[ll,] <- c(sigma2_g,sigma2_h)
      }


      ##########  adaptively adjust the proposal variance based on the acceptance probability of the current MCMC chain  ##########
      if( (l>2000) && (l%%2000==1) )
      {
         var_g=2.4^2*(cov(theta.all[1:(l-1),idx_g])+1e-10*diag(rep(1,length(idx_g))))/length(idx_g)
         var_h=2.4^2*(cov(theta.all[1:(l-1),-idx_g])+1e-10*diag(rep(1,length(idx_g))))/length(idx_g)

         if(acpt_g/l < 0.1)
         {
            var_g <- (multiplier3^2)*var_g
         }

         else if(acpt_g/l < 0.2)
         {
            var_g <- (multiplier2^2)*var_g
         }

         if(acpt_h/l < 0.1)
         {
            var_h <- (multiplier3^2)*var_h
         }

         else if(acpt_h/l < 0.2)
         {
            var_h <- (multiplier2^2)*var_h
         }
      }

   }


   ########################################### MCMC diagnostics #####################################
   geweke.theta=rep(NA,ncol(theta.res))

   for(t in 1:ncol(theta.res))
   {
	geweke.theta[t] = geweke.diag(theta.res[,t],frac1=0.3,frac2=0.3)[[1]]
   }

   geweke.sigma2=rep(NA,2)

   for(t in 1:2)
   {
	geweke.sigma2[t] = geweke.diag(sigma2.res[,t],frac1=0.3,frac2=0.3)[[1]]
   }


   ########################################### return the result ##############################

   geweke.conv <- !( max(abs(geweke.theta))>4.25 | max(abs(geweke.sigma2))>4.25 )

   if(!is.na(geweke.conv) & geweke.conv)
   {
   	return(list(success=TRUE, theta=theta.res, sigma2=sigma2.res,				# MCMC samples 
		 	acpt_g=acpt_g/L, acpt_h=acpt_h/L,							# acceptance rates
		 	geweke_theta=geweke.theta, geweke_sigma2=geweke.sigma2 ))			# Geweke test statistics 
   }

   else
   {
	return(list(success=FALSE)) 
   }

}


