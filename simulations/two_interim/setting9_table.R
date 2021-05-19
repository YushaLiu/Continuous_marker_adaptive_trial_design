# Direct to the folder where the output is saved
setwd("./simulations/two_interim/output/setting9")

library(xtable)

# n = 500, threshold for marker dichotomization = 0.2, two interim checks, p_eff_group = c(0.975, 0.99, 0.99)
res <- readRDS("trial_results.Rds")

fail.MCMC <- res$fail.MCMC
stop.eff.all.Bayes <- res$stop.eff.all.Bayes
stop.eff.pos.Bayes <- res$stop.eff.pos.Bayes
stop.eff.neg.Bayes <- res$stop.eff.neg.Bayes
stop.fut.all.Bayes <- res$stop.fut.all.Bayes
stop.fut.pos.Bayes <- res$stop.fut.pos.Bayes
stop.fut.neg.Bayes <- res$stop.fut.neg.Bayes
continue.all.Bayes <- res$continue.all.Bayes
continue.pos.Bayes <- res$continue.pos.Bayes
continue.neg.Bayes <- res$continue.neg.Bayes
sampsize.Bayes <- res$sampsize.Bayes
sampsize.pos.Bayes <- res$sampsize.pos.Bayes
sampsize.neg.Bayes <- res$sampsize.neg.Bayes
prev.Bayes <- res$marker.prevalence.Bayes
final.Bayes <- res$final.Bayes
marker.value.Bayes <- res$marker.value.Bayes
marker.status.Bayes <- res$marker.status.Bayes

stop.eff.all.freq <- res$stop.eff.all.freq
stop.eff.pos.freq <- res$stop.eff.pos.freq
stop.eff.neg.freq <- res$stop.eff.neg.freq
stop.fut.all.freq <- res$stop.fut.all.freq
stop.fut.pos.freq <- res$stop.fut.pos.freq
stop.fut.neg.freq <- res$stop.fut.neg.freq
continue.all.freq <- res$continue.all.freq
continue.pos.freq <- res$continue.pos.freq
continue.neg.freq <- res$continue.neg.freq
sampsize.freq <- res$sampsize.freq
sampsize.pos.freq <- res$sampsize.pos.freq
sampsize.neg.freq <- res$sampsize.neg.freq
prev.freq <- res$marker.prevalence.freq
final.freq <- res$final.freq
marker.value.freq <- res$marker.value.freq
marker.status.freq <- res$marker.status.freq


# check whether the Bayesian design fails for any iteration due to slow convergence of MCMC
sum(fail.MCMC)
# [1] 0



# logical indicator whether the trial still continues at the first/second interim analysis
continue.1.Bayes <- (continue.all.Bayes[,1]==1 | continue.pos.Bayes[,1]==1 | continue.neg.Bayes[,1]==1)
continue.2.Bayes <- (continue.all.Bayes[,2]==1 | continue.pos.Bayes[,2]==1 | continue.neg.Bayes[,2]==1)
continue.1.freq <- (continue.all.freq[,1]==1 | continue.pos.freq[,1]==1 | continue.neg.freq[,1]==1)
continue.2.freq <- (continue.all.freq[,2]==1 | continue.pos.freq[,2]==1 | continue.neg.freq[,2]==1)


# calculate true positive rate and true negative rate
fnc.trt <- function(x) {
  if(x<=0.5) {
    log(0.4)*exp(30*(x-0.3))/(1+exp(30*(x-0.3)))
  }
  else {
    -log(0.4)*exp(30*(x-0.7))/(1+exp(30*(x-0.7)))+log(0.4)
  }
}

TPR.all.Bayes <- rep(NA, length(final.Bayes))
TPR.sub.Bayes <- rep(NA, length(final.Bayes))

TPR.all.freq <- rep(NA, 1000)
TPR.sub.freq <- rep(NA, 1000)

TNR.all.Bayes <- rep(NA, length(final.Bayes))
TNR.sub.Bayes <- rep(NA, length(final.Bayes))

TNR.all.freq <- rep(NA, 1000)
TNR.sub.freq <- rep(NA, 1000)

for(i in 1:length(final.Bayes)) {
  if(final.Bayes[i]=="eff_all") {
    TPR.all.Bayes[i] <- 1
    TNR.all.Bayes[i] <- 0
  }
  else if(final.Bayes[i]=="eff_no") {
    TPR.all.Bayes[i] <- 0
    TNR.all.Bayes[i] <- 1
  }
  else {
    x.Bayes <- marker.value.Bayes[i,]
    idx1.Bayes <- which( unlist(lapply(x.Bayes,fnc.trt)) <= log(0.8) )
    idx2.Bayes <- which( unlist(lapply(x.Bayes,fnc.trt)) > log(0.95) )
    
    if(final.Bayes[i]=="eff_pos") {
      TPR.sub.Bayes[i] <- mean(marker.status.Bayes[i,idx1.Bayes]==1)
      TNR.sub.Bayes[i] <- mean(marker.status.Bayes[i,idx2.Bayes]==0)
    }
    else if (final.Bayes[i]=="eff_neg") {
      TPR.sub.Bayes[i] <- mean(marker.status.Bayes[i,idx1.Bayes]==0)
      TNR.sub.Bayes[i] <- mean(marker.status.Bayes[i,idx2.Bayes]==1)
    }				
  }	
}


for(i in 1:length(final.freq)) {
  if(final.freq[i]=="eff_all") {
    TPR.all.freq[i] <- 1
    TNR.all.freq[i] <- 0
  }
  else if(final.freq[i]=="eff_no") {
    TPR.all.freq[i] <- 0
    TNR.all.freq[i] <- 1
  }
  else {
    x.freq <- marker.value.freq[i,]
    idx1.freq <- which( unlist(lapply(x.freq,fnc.trt)) <= log(0.8) )
    idx2.freq <- which( unlist(lapply(x.freq,fnc.trt)) > log(0.95) )
    
    if(final.freq[i]=="eff_pos") {
      TPR.sub.freq[i] <- mean(marker.status.freq[i,idx1.freq]==1)
      TNR.sub.freq[i] <- mean(marker.status.freq[i,idx2.freq]==0)
    }
    else if (final.freq[i]=="eff_neg") {
      TPR.sub.freq[i] <- mean(marker.status.freq[i,idx1.freq]==0)
      TNR.sub.freq[i] <- mean(marker.status.freq[i,idx2.freq]==1)
    }				
  }
}



# create the result table
result <- matrix(NA, nrow=14, ncol=14)
rownames(result) <- c("No Marker Detected", "Overall Efficacy", "Overall Futility", "Overall Continue",
                      "Marker Detected", "M+ Eff/M- Eff", "M+ Eff/M- Cont", "M+ Eff/M- Fut", "M+ Cont/M- Eff",
                      "M+ Cont/M- Cont", "M+ Cont/M- Fut", "M+ Fut/M- Eff", "M+ Fut/M- Cont", "M+ Fut/M- Fut")
colnames(result) <- c("Interim1 Bayes", "Interim1 Freq", "Interim2 Bayes", "Interim2 Freq", "Final Bayes", "Final Freq", 
                      "Overall Bayes", "Overall Freq", "Sample Size Bayes", "Sample Size Freq", 
                      "TPR Bayes", "TPR Freq", "TNR Bayes", "TNR Freq")



# Interim and final analysis decision rate

# no marker detected: efficacy, futility and continue
result[2,c(1,3,5)] <- apply(stop.eff.all.Bayes,2,mean)
result[3,c(1,3)] <- apply(stop.fut.all.Bayes,2,mean)
result[3,5] <- mean(is.na(prev.Bayes) & apply(stop.eff.all.Bayes,1,sum)==0 & apply(stop.fut.all.Bayes,1,sum)==0)
result[4,c(1,3)] <- apply(continue.all.Bayes,2,mean)

result[2,c(2,4,6)] <- apply(stop.eff.all.freq,2,mean)
result[3,c(2,4)] <- apply(stop.fut.all.freq,2,mean)
result[3,6] <- mean(is.na(prev.freq) & apply(stop.eff.all.freq,1,sum)==0 & apply(stop.fut.all.freq,1,sum)==0)
result[4,c(2,4)] <- apply(continue.all.freq,2,mean)


# marker detected: efficacy, futility and continue
result[6,1]  <- mean(stop.eff.pos.Bayes[,1]>0 & stop.eff.neg.Bayes[,1]>0)
result[7,1]  <- mean(stop.eff.pos.Bayes[,1]>0 & continue.neg.Bayes[,1]>0)
result[8,1]  <- mean(stop.eff.pos.Bayes[,1]>0 & stop.fut.neg.Bayes[,1]>0)
result[9,1]  <- mean(continue.pos.Bayes[,1]>0 & stop.eff.neg.Bayes[,1]>0)
result[10,1] <- mean(continue.pos.Bayes[,1]>0 & continue.neg.Bayes[,1]>0)
result[11,1] <- mean(continue.pos.Bayes[,1]>0 & stop.fut.neg.Bayes[,1]>0)
result[12,1] <- mean(stop.fut.pos.Bayes[,1]>0 & stop.eff.neg.Bayes[,1]>0)
result[13,1] <- mean(stop.fut.pos.Bayes[,1]>0 & continue.neg.Bayes[,1]>0)
result[14,1] <- mean(stop.fut.pos.Bayes[,1]>0 & stop.fut.neg.Bayes[,1]>0)

result[6,3]  <- mean(continue.1.Bayes==1 & continue.2.Bayes==0 & apply(stop.eff.pos.Bayes,1,sum)>0 & apply(stop.eff.neg.Bayes,1,sum)>0)
result[7,3]  <- mean(continue.2.Bayes==1 & apply(stop.eff.pos.Bayes[,1:2],1,sum)>0)
result[8,3]  <- mean(continue.1.Bayes==1 & continue.2.Bayes==0 & apply(stop.eff.pos.Bayes,1,sum)>0 & apply(stop.fut.neg.Bayes,1,sum)>0)
result[9,3]  <- mean(continue.2.Bayes==1 & apply(stop.eff.neg.Bayes[,1:2],1,sum)>0)
result[10,3] <- mean(continue.pos.Bayes[,2]>0 & continue.neg.Bayes[,2]>0)
result[11,3] <- mean(continue.2.Bayes==1 & apply(stop.fut.neg.Bayes,1,sum)>0)
result[12,3] <- mean(continue.1.Bayes==1 & continue.2.Bayes==0 & apply(stop.fut.pos.Bayes,1,sum)>0 & apply(stop.eff.neg.Bayes,1,sum)>0)
result[13,3] <- mean(continue.2.Bayes==1 & apply(stop.fut.pos.Bayes,1,sum)>0)
result[14,3] <- mean(continue.1.Bayes==1 & continue.2.Bayes==0 & apply(stop.fut.pos.Bayes,1,sum)>0 & apply(stop.fut.neg.Bayes,1,sum)>0)

result[6,5] <- mean(continue.2.Bayes==1 & (!is.na(prev.Bayes)) & apply(stop.eff.pos.Bayes,1,sum)>0 & apply(stop.eff.neg.Bayes,1,sum)>0)
result[8,5] <- mean(continue.2.Bayes==1 & (!is.na(prev.Bayes)) & apply(stop.eff.pos.Bayes,1,sum)>0 & apply(stop.eff.neg.Bayes,1,sum)==0)
result[12,5] <- mean(continue.2.Bayes==1 & (!is.na(prev.Bayes)) & apply(stop.eff.pos.Bayes,1,sum)==0 & apply(stop.eff.neg.Bayes,1,sum)>0)
result[14,5] <- mean(continue.2.Bayes==1 & (!is.na(prev.Bayes)) & apply(stop.eff.pos.Bayes,1,sum)==0 & apply(stop.eff.neg.Bayes,1,sum)==0)


result[6,2]  <- mean(stop.eff.pos.freq[,1]>0 & stop.eff.neg.freq[,1]>0)
result[7,2]  <- mean(stop.eff.pos.freq[,1]>0 & continue.neg.freq[,1]>0)
result[8,2]  <- mean(stop.eff.pos.freq[,1]>0 & stop.fut.neg.freq[,1]>0)
result[9,2]  <- mean(continue.pos.freq[,1]>0 & stop.eff.neg.freq[,1]>0)
result[10,2] <- mean(continue.pos.freq[,1]>0 & continue.neg.freq[,1]>0)
result[11,2] <- mean(continue.pos.freq[,1]>0 & stop.fut.neg.freq[,1]>0)
result[12,2] <- mean(stop.fut.pos.freq[,1]>0 & stop.eff.neg.freq[,1]>0)
result[13,2] <- mean(stop.fut.pos.freq[,1]>0 & continue.neg.freq[,1]>0)
result[14,2] <- mean(stop.fut.pos.freq[,1]>0 & stop.fut.neg.freq[,1]>0)

result[6,4]  <- mean(continue.1.freq==1 & continue.2.freq==0 & apply(stop.eff.pos.freq,1,sum)>0 & apply(stop.eff.neg.freq,1,sum)>0)
result[7,4]  <- mean(continue.2.freq==1 & apply(stop.eff.pos.freq[,1:2],1,sum)>0)
result[8,4]  <- mean(continue.1.freq==1 & continue.2.freq==0 & apply(stop.eff.pos.freq,1,sum)>0 & apply(stop.fut.neg.freq,1,sum)>0)
result[9,4]  <- mean(continue.2.freq==1 & apply(stop.eff.neg.freq[,1:2],1,sum)>0)
result[10,4] <- mean(continue.pos.freq[,2]>0 & continue.neg.freq[,2]>0)
result[11,4] <- mean(continue.2.freq==1 & apply(stop.fut.neg.freq,1,sum)>0)
result[12,4] <- mean(continue.1.freq==1 & continue.2.freq==0 & apply(stop.fut.pos.freq,1,sum)>0 & apply(stop.eff.neg.freq,1,sum)>0)
result[13,4] <- mean(continue.2.freq==1 & apply(stop.fut.pos.freq,1,sum)>0)
result[14,4] <- mean(continue.1.freq==1 & continue.2.freq==0 & apply(stop.fut.pos.freq,1,sum)>0 & apply(stop.fut.neg.freq,1,sum)>0)

result[6,6] <- mean(continue.2.freq==1 & (!is.na(prev.freq)) & apply(stop.eff.pos.freq,1,sum)>0 & apply(stop.eff.neg.freq,1,sum)>0)
result[8,6] <- mean(continue.2.freq==1 & (!is.na(prev.freq)) & apply(stop.eff.pos.freq,1,sum)>0 & apply(stop.eff.neg.freq,1,sum)==0)
result[12,6] <- mean(continue.2.freq==1 & (!is.na(prev.freq)) & apply(stop.eff.pos.freq,1,sum)==0 & apply(stop.eff.neg.freq,1,sum)>0)
result[14,6] <- mean(continue.2.freq==1 & (!is.na(prev.freq)) & apply(stop.eff.pos.freq,1,sum)==0 & apply(stop.eff.neg.freq,1,sum)==0)


# overall trial decision rate
result[2,7] <- mean(final.Bayes=="eff_all")
result[3,7] <- mean(final.Bayes=="eff_no")
result[8,7] <- mean(final.Bayes=="eff_pos")
result[12,7] <- mean(final.Bayes=="eff_neg")

result[2,8] <- mean(final.freq=="eff_all")
result[3,8] <- mean(final.freq=="eff_no")
result[8,8] <- mean(final.freq=="eff_pos")
result[12,8] <- mean(final.freq=="eff_neg")


# dichotomization rate
result[1,1] <- result[2,1]+result[3,1]+result[4,1]
result[1,3] <- result[2,3]+result[3,3]+result[4,3]
result[1,5] <- result[2,5]+result[3,5]
result[5,1] <- sum(result[6:14,1])
result[5,3] <- sum(result[6:14,3])
result[5,5] <- sum(result[c(6,8,12,14),5])
result[1,7] <- mean(final.Bayes=="eff_all" | final.Bayes=="eff_no")
result[5,7] <- mean(final.Bayes=="eff_pos" | final.Bayes=="eff_neg")

result[1,2] <- result[2,2]+result[3,2]+result[4,2]
result[1,4] <- result[2,4]+result[3,4]+result[4,4]
result[1,6] <- result[2,6]+result[3,6]
result[5,2] <- sum(result[6:14,2])
result[5,4] <- sum(result[6:14,4])
result[5,6] <- sum(result[c(6,8,12,14),6])
result[1,8] <- mean(final.freq=="eff_all" | final.freq=="eff_no")
result[5,8] <- mean(final.freq=="eff_pos" | final.freq=="eff_neg")



# sample size
result[1,9] <- round(mean(sampsize.Bayes[final.Bayes=="eff_all" | final.Bayes=="eff_no"]))
result[5,9] <- round(mean(sampsize.pos.Bayes[final.Bayes=="eff_pos" | final.Bayes=="eff_neg"]))
result[6,9] <- round(mean(sampsize.neg.Bayes[final.Bayes=="eff_pos" | final.Bayes=="eff_neg"]))

result[1,10] <- round(mean(sampsize.freq[final.freq=="eff_all" | final.freq=="eff_no"]))
result[5,10] <- round(mean(sampsize.pos.freq[final.freq=="eff_pos" | final.freq=="eff_neg"]))
result[6,10] <- round(mean(sampsize.neg.freq[final.freq=="eff_pos" | final.freq=="eff_neg"]))


# TPR and TNR
result[1,11] <- round(mean(TPR.all.Bayes[!is.nan(TPR.all.Bayes) & !is.na(TPR.all.Bayes)]),3)
result[1,12] <- round(mean(TPR.all.freq[!is.nan(TPR.all.freq) & !is.na(TPR.all.freq)]),3)

result[5,11] <- round(mean(TPR.sub.Bayes[!is.nan(TPR.sub.Bayes) & !is.na(TPR.sub.Bayes)]),3)
result[5,12] <- round(mean(TPR.sub.freq[!is.nan(TPR.sub.freq) & !is.na(TPR.sub.freq)]),3)

result[1,13] <- round(mean(TNR.all.Bayes[!is.nan(TNR.all.Bayes) & !is.na(TNR.all.Bayes)]),3)
result[1,14] <- round(mean(TNR.all.freq[!is.nan(TNR.all.freq) & !is.na(TNR.all.freq)]),3)

result[5,13] <- round(mean(TNR.sub.Bayes[!is.nan(TNR.sub.Bayes) & !is.na(TNR.sub.Bayes)]),3)
result[5,14] <- round(mean(TNR.sub.freq[!is.nan(TNR.sub.freq) & !is.na(TNR.sub.freq)]),3)




# output results
result2 <- xtable(result, digits=c(4,3,3,3,3,3,3,3,3,0,0,3,3,3,3))
print(result2, include.rownames = TRUE, booktabs = TRUE)