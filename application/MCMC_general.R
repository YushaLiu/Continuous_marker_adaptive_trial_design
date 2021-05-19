##########################################################################################
#### Perform MCMC to generate posterior samples of penalized spline coefficients and regularization parameters
#### This function is called by Bayes_penalized_splines, not directly called by the user
##########################################################################################

MCMC <- function(data, formula, ncolz, prior, multiplier1, multiplier2, multiplier3, burnin, B, thin){
  
  library(MASS)
  library(survival)  
  library(Matrix)
  library(coda)
  
  # calculate the number of markers from the input data
  J <- (ncol(data)-3)/(ncolz+1)/2

  # get the data frame with regressors only
  data.X <- data[, -c(1:2)]
  

  ####################################### MCMC initialization ##############################
  # determine the initial values for the PH regression coefficients
  initfit <- coxph(formula=formula, data=data)
  theta <- initfit$coef
  
  # determine the initial value and random walk proposal variance for the baseline hazard lambda0,
  # which has a conjugate Gamma full conditional
  alpha_lambda0 <- sum(data$status)
  beta_lambda0 <- sum(exp(as.matrix(data.X) %*% theta)*data$futime)
  lambda0 <- (alpha_lambda0-1)/beta_lambda0
  var_lambda0 <- alpha_lambda0/(beta_lambda0^2)
   
  # determine the initial proposal variance for Metropolis random walk of the PH regression coefficients for each marker
  # determine the initial value of sigma2_g and sigma2_h, i.e., regularization parameters for each marker
  var_g <- list(NULL)
  var_h <- list(NULL)
  sigma2_g <- rep(NA, J)
  sigma2_h <- rep(NA, J)
   
  # loop over each marker
  for(j in 1:J){
    if(j==1){
      # prognostic marker effects
      idx.jg.start <- 1
      idx.jg.end <- idx.jg.start + ncolz
      data_jg <- data.X[, idx.jg.start:idx.jg.end]
      theta_jg <- theta[idx.jg.start:idx.jg.end]
      u_jg <- theta_jg[-c(1)]
      data_jg$offset <- as.numeric(as.matrix(data.X[,-c(idx.jg.start:idx.jg.end)]) %*% theta[-c(idx.jg.start:idx.jg.end)])
      
      # predictive marker effects
      idx.jh.start <- idx.jg.end + 1
      idx.jh.end <- idx.jh.start + ncolz + 1
      data_jh <- data.X[, idx.jh.start:idx.jh.end]
      theta_jh <- theta[idx.jh.start:idx.jh.end]
      u_jh <- theta_jh[-c(1,2)]
      data_jh$offset <- as.numeric(as.matrix(data.X[,-c(idx.jh.start:idx.jh.end)]) %*% theta[-c(idx.jh.start:idx.jh.end)])
    }
    
    else{
      # prognostic marker effects
      idx.jg.start <- idx.jh.end + 1
      idx.jg.end <- idx.jg.start + ncolz
      data_jg <- data.X[, idx.jg.start:idx.jg.end]
      theta_jg <- theta[idx.jg.start:idx.jg.end]
      u_jg <- theta_jg[-c(1)]
      data_jg$offset <- as.numeric(as.matrix(data.X[,-c(idx.jg.start:idx.jg.end)]) %*% theta[-c(idx.jg.start:idx.jg.end)])   
      
      # predictive marker effects
      idx.jh.start <- idx.jg.end + 1
      idx.jh.end <- idx.jh.start + ncolz
      data_jh <- data.X[, idx.jh.start:idx.jh.end]
      theta_jh <- theta[idx.jh.start:idx.jh.end]
      u_jh <- theta_jh[-c(1)]
      data_jh$offset <- as.numeric(as.matrix(data.X[,-c(idx.jh.start:idx.jh.end)]) %*% theta[-c(idx.jh.start:idx.jh.end)])
    }
    
    # determine the initial proposal variance for Metropolis random walk for coefficients related to prognostic marker effects
    formula1 <- as.formula(paste("Surv(futime,status) ~ ", paste(colnames(data_jg), collapse="+")))
    initfit_jg <- coxph(formula=formula1, data=cbind(data[,1:2],data_jg))
    if(j==1){
      var_g[[j]] <- (multiplier1^2)*as.matrix(bdiag(var_lambda0*diag(1)/(lambda0^2), initfit_jg$var[-c(ncol(data_jg)),-c(ncol(data_jg))]))
    }
    else{
      var_g[[j]] <- (multiplier1^2)*initfit_jg$var[-c(ncol(data_jg)),-c(ncol(data_jg))]
    }
    
    # determine the initial value of sigma2_jg
    sigma2_g[j] <- 1/rgamma(1, shape=ncolz/2+prior$alpha, rate=sum(u_jg^2)/2+prior$beta)
    
    # determine the initial proposal variance for Metropolis random walk for coefficients related to predictive marker effects
    formula2 <- as.formula(paste("Surv(futime,status) ~ ", paste(colnames(data_jh), collapse="+") ))
    initfit_jh <- coxph(formula=formula2, data=cbind(data[,1:2],data_jh))
    var_h[[j]] <- (multiplier1^2)*initfit_jh$var[-c(ncol(data_jh)),-c(ncol(data_jh))]
    
    # determine the initial value of sigma2_jh
    sigma2_h[j] <- 1/rgamma(1, shape=ncolz/2+prior$alpha, rate=sum(u_jh^2)/2+prior$beta)
  } 
  
  
  # add the intercept term, i.e., log(lambda0)
  theta <- c(log(lambda0), theta)
  names(theta) <- c("intercept", colnames(data.X))
  data.X <- cbind(1, as.matrix(data.X))
  colnames(data.X) <- names(theta)

   
  ########################################### MCMC starts #####################################
  L <- burnin+B*thin
  
  theta.all <- matrix(NA, nrow=L, ncol=length(theta))
  colnames(theta.all) <- colnames(data.X)
  theta.res <- matrix(NA, nrow=B, ncol=length(theta))
  colnames(theta.res) <- colnames(data.X)
  sigma2.g.all <- matrix(NA, nrow=L, ncol=J)
  sigma2.g.res <- matrix(NA, nrow=B, ncol=J)
  sigma2.h.all <- matrix(NA, nrow=L, ncol=J)
  sigma2.h.res <- matrix(NA, nrow=B, ncol=J)

  ll <- 0
  acpt_g <- rep(0, J)
  acpt_h <- rep(0, J) 
  
  for(l in 1:L){
    ### update theta_jg, sigma2_jg, theta_jh, sigma2_jh for each marker j
    for(j in 1:J){
      # prognostic marker effects 
      if(j==1){
        # indices for theta_jg
        idx.jg.start <- 1
        idx.jg.end <- idx.jg.start + ncolz + 1
        
        # current value of theta_jg
        theta_jg <- theta[idx.jg.start:idx.jg.end]
        u_jg <- theta_jg[-c(1,2)]
      }
      
      else{
        # indices for theta_jg
        idx.jg.start <- idx.jh.end + 1
        idx.jg.end <- idx.jg.start + ncolz
        
        # current value of theta_jg
        theta_jg <- theta[idx.jg.start:idx.jg.end]
        u_jg <- theta_jg[-c(1)]
      }
      
      data_jg <- data.X[, idx.jg.start:idx.jg.end]
      coef_jg <- as.numeric(data_jg %*% theta_jg)
      coef_all <- as.numeric(data.X %*% theta)
      
      # propose new value for theta_jg
      newtheta_jg <- mvrnorm(1, theta_jg, var_g[[j]])
      if(j==1){
        newu_jg <- newtheta_jg[-c(1,2)]
      }
      else{
        newu_jg <- newtheta_jg[-c(1)]
      }
      
      newtheta <- theta
      newtheta[idx.jg.start:idx.jg.end] <- newtheta_jg
      newcoef_jg <- as.numeric(data_jg %*% newtheta_jg)
      newcoef_all <- as.numeric(data.X %*% newtheta)
      
      # calculate acceptance probability
      loglik <- sum(coef_jg*data$status - exp(coef_all)*data$futime)
      newloglik <- sum(newcoef_jg*data$status - exp(newcoef_all)*data$futime)
      log_acpt <- newloglik - 0.5*sum(newu_jg^2)/sigma2_g[j] - loglik + 0.5*sum(u_jg^2)/sigma2_g[j]
      
      # decide whether to accept the proposed value for theta_jg
      if(log(runif(1)) < log_acpt){
        theta <- newtheta
        u_jg <- newu_jg
        acpt_g[j] <- acpt_g[j] + 1
      }
      
      # update sigma2_jg
      sigma2_g[j] <- 1/rgamma(1, shape=ncolz/2+prior$alpha, rate=sum(u_jg^2)/2+prior$beta)
      
      
      ### predictive marker effects
      if(j==1){
        # indices for theta_jh
        idx.jh.start <- idx.jg.end + 1
        idx.jh.end <- idx.jh.start + ncolz + 1
        
        # current values of theta_jh
        theta_jh <- theta[idx.jh.start:idx.jh.end]
        u_jh <- theta_jh[-c(1,2)]
      }
      
      else{
        # indices for theta_jg
        idx.jh.start <- idx.jg.end + 1
        idx.jh.end <- idx.jh.start + ncolz
        
        # current values of theta_jh
        theta_jh <- theta[idx.jh.start:idx.jh.end]
        u_jh <- theta_jh[-c(1)]
      }
      
      data_jh <- data.X[, idx.jh.start:idx.jh.end]
      coef_jh <- as.numeric(data_jh %*% theta_jh)
      coef_all <- as.numeric(data.X %*% theta)
      
      # propose new value for theta_jh
      newtheta_jh <- mvrnorm(1, theta_jh, var_h[[j]])
      if(j==1){
        newu_jh <- newtheta_jh[-c(1,2)]
      }
      else{
        newu_jh <- newtheta_jh[-c(1)]
      }
      
      newtheta <- theta
      newtheta[idx.jh.start:idx.jh.end] <- newtheta_jh
      newcoef_jh <- as.numeric(data_jh %*% newtheta_jh)
      newcoef_all <- as.numeric(data.X %*% newtheta)
      
      # calculate acceptance probability
      loglik <- sum(coef_jh*data$status - exp(coef_all)*data$futime)
      newloglik <- sum(newcoef_jh*data$status - exp(newcoef_all)*data$futime)
      log_acpt <- newloglik - 0.5*sum(newu_jh^2)/sigma2_h[j] - loglik + 0.5*sum(u_jh^2)/sigma2_h[j]
    
      # decide whether to accept the proposed value for theta_jh
      if(log(runif(1)) < log_acpt){
        theta <- newtheta
        u_jh <- newu_jh
        acpt_h[j] <- acpt_h[j] + 1
      }
      
      # update sigma2_jh
      sigma2_h[j] <- 1/rgamma(1, shape=ncolz/2+prior$alpha, rate=sum(u_jh^2)/2+prior$beta)
    }
    
    
    ############################# save the values of the new iteration ###############################
    theta.all[l,] <- theta
    sigma2.g.all[l,] <- sigma2_g
    sigma2.h.all[l,] <- sigma2_h
      
    if( (l>burnin) && ((l-burnin)%%thin==0) ){
      ll <- ll+1
      theta.res[ll,] <- theta
      sigma2.g.res[ll,] <- sigma2_g
      sigma2.h.res[ll,] <- sigma2_h
    }
    
    
    ############# adaptively adjust the proposal variance based on the cumulative acceptance probability #############
    if( (l>2000) && (l%%2000==1) ){
      for(j in 1:J){
        ### prognostic marker effects
        if(j==1){
          idx.jg.start <- 1
          idx.jg.end <- idx.jg.start + ncolz + 1
        }
        else{
          idx.jg.start <- idx.jh.end + 1
          idx.jg.end <- idx.jg.start + ncolz
        }

        # calculate the covariance matrix of theta_jg
        var_jg <- 2.4^2*(cov(theta.all[1:(l-1), idx.jg.start:idx.jg.end]) + 
                           1e-10*diag(rep(1,length(idx.jg.start:idx.jg.end))))/length(idx.jg.start:idx.jg.end)
        
        # adjust adaptively
        if(acpt_g[j]/l < 0.1){
          var_jg <- (multiplier3^2)*var_jg
        }
        
        else if(acpt_g[j]/l < 0.2){
          var_jg <- (multiplier2^2)*var_jg
        }
        
        var_g[[j]] <- var_jg
        
        
        ### predictive marker effects
        if(j==1){
          idx.jh.start <- idx.jg.end + 1
          idx.jh.end <- idx.jh.start + ncolz + 1
        }
        else{
          idx.jh.start <- idx.jg.end + 1
          idx.jh.end <- idx.jh.start + ncolz
        }
        
        # calculate the covariance matrix of theta_jh
        var_jh <- 2.4^2*(cov(theta.all[1:(l-1), idx.jh.start:idx.jh.end]) + 
                           1e-10*diag(rep(1,length(idx.jh.start:idx.jh.end))))/length(idx.jh.start:idx.jh.end)
        
        # adjust adaptively
        if(acpt_h[j]/l < 0.1){
          var_jh <- (multiplier3^2)*var_jh
        }
        
        else if(acpt_h[j]/l < 0.2){
          var_jh <- (multiplier2^2)*var_jh
        }
        
        var_h[[j]] <- var_jh
      }
    }
  }


  ############################################ MCMC diagnostics ###############################################
  geweke.theta <- rep(NA, ncol(theta.res))
  names(geweke.theta) <- colnames(theta.res)
  
  for(t in 1:ncol(theta.res)){
    geweke.theta[t] <- geweke.diag(theta.res[,t], frac1=0.25, frac2=0.25)[[1]]
  }

  geweke.sigma2.g <- rep(NA, J)
  names(geweke.sigma2.g) <- paste0("x", 1:J)

  for(t in 1:J){
    geweke.sigma2.g[t] <- geweke.diag(sigma2.g.res[,t], frac1=0.25, frac2=0.25)[[1]]
  }
  
  geweke.sigma2.h <- rep(NA, J)
  names(geweke.sigma2.h) <- paste0("x", 1:J)
  
  for(t in 1:J){
    geweke.sigma2.h[t] <- geweke.diag(sigma2.h.res[,t], frac1=0.25, frac2=0.25)[[1]]
  }


  ########################################### return the result ################################################
  geweke.conv <- !(max(abs(geweke.theta))>4 | max(abs(geweke.sigma2.g))>4 | max(abs(geweke.sigma2.h))>4)

  if(!is.na(geweke.conv) & geweke.conv){
    return(list(success=TRUE, theta=theta.res, sigma2_g=sigma2.g.res, sigma2_h=sigma2.h.res,				# MCMC samples 
                acpt_g=acpt_g/L, acpt_h=acpt_h/L,							# acceptance rates
                geweke_theta=geweke.theta, geweke_sigma2_g=geweke.sigma2.g, geweke_sigma2_h=geweke.sigma2.h))			# Geweke test statistics 
  }

  else{
    return(list(success=FALSE)) 
  }

}


