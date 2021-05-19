##########################################################################################
#### Calculate (1-alpha)*100% simultaneous band based on the matrix of posterior samples MCMC_P
#### This function is called by Bayes_penalized_splines, not directly called by the user
##########################################################################################

jointband <- function(MCMC_P,alpha){
   B <- dim(MCMC_P)[1]

   mean_P <- apply(MCMC_P,2,mean)
   sd_P <- apply(MCMC_P,2,sd)

   z_P <- rep(NA,B)

   keep <- (!is.na(sd_P)) & (sd_P!=0) 
    
   for(j in 1:B)
   {
      z_P[j] <- max(abs((MCMC_P[j,keep]-mean_P[keep])/sd_P[keep]))
   }

   upr_CI <- mean_P + quantile(z_P,1-alpha)*sd_P
   lwr_CI <- mean_P - quantile(z_P,1-alpha)*sd_P

   return(list(upr_CI=upr_CI,lwr_CI=lwr_CI))
}



##########################################################################################
#### Calculate the distance between a given vector of marker values and each row of a matrix of marker values 
#### This function is called by Bayes_adaptive_design and Bayes_penalized_splines, not directly called by the user
##########################################################################################

dist_marker <- function(x0, X){
#   n <- nrow(X)
#   dist.X <- rep(NA, n)
   
#   for(i in 1:n){
#      dist.X[i] <- sqrt(sum((X[i,]-x0)^2))
#   }
 
   dist.X <- sweep(as.matrix(X), 2, as.numeric(x0))
   dist.X <- sqrt(rowSums(dist.X^2))
   return(dist.X)
}