#######################################################################################
#### This function constructs the O'Sullivan penalized splines design matrix (nonlinear term) 
#### based on the paper "On semiparametric regression with O'Sullivan penalized splines" (2008) by Wand and Ormerod
#### This function is called by Bayes_penalized_splines, not directly called by the user
### data = a data frame with the observations arranged by row, and including the columns below (in the following order):
###        1) futime: time to event or right censoring, whichever happens first
###        2) status: the censoring status indicator (=1 for event; =0 for right censored)
###        3) trt: the group indicator (=1 for experimental group; =0 for control group)
### X = a data frame giving continuous marker values, with rows as patients and columns as biomarkers 
### xrange is a matrix with 2 columns which gives the range of marker values, one marker per row
### numIntKnots is the number of interior knots (should scale with number of events)
### censor denotes whether marker values of right censored patients are included to decide the interior knots locations
#######################################################################################


Get_splines <- function(data, X, xrange=cbind(rep(0, ncol(X)), rep(1, ncol(X))), numIntKnots, censor=FALSE){
   
   library(splines)

   # Set up the design matrix and related quantities for each biomarker
   J <- ncol(X)
   stabCheck_idx <- rep(FALSE, J)
   
   for(j in 1:J){
      a <- xrange[j,1]
      b <- xrange[j,2]
      x <- as.numeric(X[,j])
      
      # Determine the location of interior knots
      if(censor)
      { 
         intKnots <- quantile(unique(x), seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))])
      }
      
      else
      {
         intKnots <- quantile(unique(x[data$status==1]), seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))])
      }
      
      names(intKnots) <- NULL
      B <- bs(x, knots=intKnots, degree=3, Boundary.knots=c(a,b), intercept=TRUE)
      
      # Create the Omega matrix
      formOmega <- function(a,b,intKnots)
      {
         allKnots <- c(rep(a,4), intKnots, rep(b,4))
         K <- length(intKnots) 
         L <- 3*(K+8)
         xtilde <- (rep(allKnots, each=3)[-c(1,(L-1),L)] + rep(allKnots, each=3)[-c(1,2,L)])/2
         wts <- rep(diff(allKnots), each=3)*rep(c(1,4,1)/6, K+7)
         Bdd <- spline.des(allKnots, xtilde, derivs=rep(2,length(xtilde)), outer.ok=TRUE)$design 
         Omega <- t(Bdd*wts)%*%Bdd  
         return(Omega)
      }
      
      Omega <- formOmega(a,b,intKnots)
      
      # Obtain the spectral decomposition of Omega
      eigOmega <- eigen(Omega)
      
      # Obtain the matrix for linear transformation of $\\bB$ to $\\bZ$
      indsZ <- 1:(numIntKnots+2)
      UZ <- eigOmega$vectors[,indsZ]
      LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))
      
      # Perform a stability check
      indsX <- (numIntKnots+3):(numIntKnots+4)
      UX <- eigOmega$vectors[,indsX]   
      L <- cbind(UX,LZ)
      stabCheck <- t(crossprod(L,t(crossprod(L,Omega))))  
      stabCheck_idx[j] <- sum(stabCheck^2) > 1.0001*(numIntKnots+2)
      
      # Form Z matrices
      Z <- B%*%LZ
      
      # Interaction between trt and Z matrices
      trtZ <- (data$trt)*Z     
      
      if(j==1){
         design.j <- cbind(x, Z, data$trt, (data$trt)*x, trtZ)
         colnames(design.j) <- c(paste0("x", j), paste0("x", j, "_z", 1:ncol(Z)), "trt", paste0("x", j, "trt"), paste0("x", j, "_z", 1:ncol(Z), "trt"))
         data.full <- cbind(data[,1:2], design.j)
      }
      else{
         design.j <- cbind(x, Z, (data$trt)*x, trtZ)
         colnames(design.j) <- c(paste0("x", j), paste0("x", j, "_z", 1:ncol(Z)), paste0("x", j, "trt"), paste0("x", j, "_z", 1:ncol(Z), "trt"))
         data.full <- cbind(data.full, design.j)
      }
   }
   
   # Return the formula for the PH regression model
   formula <- as.formula(paste("Surv(futime,status) ~ ", paste(colnames(data.full)[-c(1:2)], collapse="+"))) 

   return(list(data=data.full, formula=formula, stabCheck=stabCheck_idx, ncolz=numIntKnots+2))
}


 
   



