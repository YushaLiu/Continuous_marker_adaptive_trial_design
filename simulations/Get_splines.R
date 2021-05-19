#######################################################################################
#### This function constructs the O'Sullivan penalized splines design matrix (nonlinear term) 
#### based on the paper "On semiparametric regression with O'Sullivan penalized splines" (2008) by Wand and Ormerod
#### This function is called by Bayes_penalized_splines, not directly called by the user
### data = a data frame with the observations arranged by row, and including the columns below (in the following order):
###        1) futime: time to event or right censoring, whichever happens first
###        2) status: the censoring status indicator (=1 for event; =0 for right censored)
###        3) x: the value of the continuous marker
###        4) trt: the group indicator (=1 for experimental group; =0 for control group)
### xrange is a vector of length 2 which gives the range of marker values
### numIntKnots is the number of interior knots (should scale with number of events)
### censor denotes whether marker values of right censored patients are included to decide the interior knots locations
#######################################################################################


Get_splines <- function(data,xrange,numIntKnots,censor)
{

   library(splines)

   # Set up the design matrix and related quantities
   a <- xrange[1]
   b <- xrange[2]
   x <- data$x
 
   if(censor)
   { 
   	intKnots <- quantile(unique(x),seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))])
   }

   else
   {
   	intKnots <- quantile(unique(x[data$status==1]),seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))])
   }

   names(intKnots) <- NULL
   B <- bs(x,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE)


   # Create the Omega matrix
   formOmega <- function(a,b,intKnots)
   {
      allKnots <- c(rep(a,4),intKnots,rep(b,4))
      K <- length(intKnots) 
      L <- 3*(K+8)
      xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+rep(allKnots,each=3)[-c(1,2,L)])/2
      wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
      Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),outer.ok=TRUE)$design 
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
   stabCheck_idx <- sum(stabCheck^2) > 1.0001*(numIntKnots+2)
        

   # Form Z matrices:
   Z <- B%*%LZ


   # Interaction between trt and Z matrices:
   trtZ <- (data$trt)*Z

   
   # Add both matrices to the data frame:
   data2 <- cbind(data[,1:3], Z, data$trt, (data$x)*(data$trt), trtZ)
   colnames(data2) <- c(colnames(data)[1:3], as.vector(outer("z", 1:ncol(Z), paste0)), "trt", "xtrt", paste0(as.vector(outer("z", 1:ncol(Z), paste0)),"trt") )

   formula <- as.formula(paste("Surv(futime,status) ~ ", paste(colnames(data2)[-c(1:2)], collapse="+") )) 

   return(list(data=data2, formula=formula, stabCheck=stabCheck_idx, ncolz=ncol(Z)))
}


 
   



