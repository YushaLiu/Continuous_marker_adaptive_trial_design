# Direct to the folder where the output is saved
setwd("./output/setting4")

library(xtable)

# n = 500, threshold for marker dichotomization = 0.1, one interim check
load("trial_results.RData")

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


sum(fail.MCMC)
[1] 0



# calculate true positive rate and true negative rate
fnc.trt <- function(x) {return(log(0.4)*as.numeric(x>=0.5))}

TPR1.all.Bayes <- rep(NA, length(final.Bayes))
TPR2.all.Bayes <- rep(NA, length(final.Bayes))

TPR1.sub.Bayes <- rep(NA, length(final.Bayes))
TPR2.sub.Bayes <- rep(NA, length(final.Bayes))

TPR1.all.freq <- rep(NA, 1000)
TPR2.all.freq <- rep(NA, 1000)

TPR1.sub.freq <- rep(NA, 1000)
TPR2.sub.freq <- rep(NA, 1000)

TNR1.all.Bayes <- rep(NA, length(final.Bayes))
TNR2.all.Bayes <- rep(NA, length(final.Bayes))

TNR1.sub.Bayes <- rep(NA, length(final.Bayes))
TNR2.sub.Bayes <- rep(NA, length(final.Bayes))

TNR1.all.freq <- rep(NA, 1000)
TNR2.all.freq <- rep(NA, 1000)

TNR1.sub.freq <- rep(NA, 1000)
TNR2.sub.freq <- rep(NA, 1000)


for(i in 1:length(final.Bayes)) {
	if(final.Bayes[i]=="eff_all") {
		TPR1.all.Bayes[i] <- 1
		TPR2.all.Bayes[i] <- 1
		TNR1.all.Bayes[i] <- 0
		TNR2.all.Bayes[i] <- 0
	}
	else if(final.Bayes[i]=="eff_no") {
		TPR1.all.Bayes[i] <- 0
		TPR2.all.Bayes[i] <- 0
		TNR1.all.Bayes[i] <- 1
		TNR2.all.Bayes[i] <- 1
	}
	else {
		x.Bayes <- marker.value.Bayes[i,]
		idx1.Bayes <- which( unlist(lapply(x.Bayes,fnc.trt)) <= log(0.8) )
		idx2.Bayes <- which( unlist(lapply(x.Bayes,fnc.trt)) <= log(0.95) )
		idx3.Bayes <- which( unlist(lapply(x.Bayes,fnc.trt)) > log(0.95) )
		idx4.Bayes <- which( unlist(lapply(x.Bayes,fnc.trt)) > log(0.99) )		

			if(final.Bayes[i]=="eff_pos") {
				TPR1.sub.Bayes[i] <- mean(marker.status.Bayes[i,idx1.Bayes]==1)
				TPR2.sub.Bayes[i] <- mean(marker.status.Bayes[i,idx2.Bayes]==1)
				TNR1.sub.Bayes[i] <- mean(marker.status.Bayes[i,idx3.Bayes]==0)
				TNR2.sub.Bayes[i] <- mean(marker.status.Bayes[i,idx4.Bayes]==0)
			}
			else if (final.Bayes[i]=="eff_neg") {
				TPR1.sub.Bayes[i] <- mean(marker.status.Bayes[i,idx1.Bayes]==0)
				TPR2.sub.Bayes[i] <- mean(marker.status.Bayes[i,idx2.Bayes]==0)
				TNR1.sub.Bayes[i] <- mean(marker.status.Bayes[i,idx3.Bayes]==1)
				TNR2.sub.Bayes[i] <- mean(marker.status.Bayes[i,idx4.Bayes]==1)
			}				
	}
}


for(i in 1:length(final.freq)) {
	if(final.freq[i]=="eff_all") {
		TPR1.all.freq[i] <- 1
		TPR2.all.freq[i] <- 1
		TNR1.all.freq[i] <- 0
		TNR2.all.freq[i] <- 0
	}
	else if(final.freq[i]=="eff_no") {
		TPR1.all.freq[i] <- 0
		TPR2.all.freq[i] <- 0
		TNR1.all.freq[i] <- 1
		TNR2.all.freq[i] <- 1
	}
	else {
		x.freq <- marker.value.freq[i,]
		idx1.freq <- which( unlist(lapply(x.freq,fnc.trt)) <= log(0.8) )
		idx2.freq <- which( unlist(lapply(x.freq,fnc.trt)) <= log(0.95) )
		idx3.freq <- which( unlist(lapply(x.freq,fnc.trt)) > log(0.95) )
		idx4.freq <- which( unlist(lapply(x.freq,fnc.trt)) > log(0.99) )		

			if(final.freq[i]=="eff_pos") {
				TPR1.sub.freq[i] <- mean(marker.status.freq[i,idx1.freq]==1)
				TPR2.sub.freq[i] <- mean(marker.status.freq[i,idx2.freq]==1)
				TNR1.sub.freq[i] <- mean(marker.status.freq[i,idx3.freq]==0)
				TNR2.sub.freq[i] <- mean(marker.status.freq[i,idx4.freq]==0)
			}
			else if (final.freq[i]=="eff_neg") {
				TPR1.sub.freq[i] <- mean(marker.status.freq[i,idx1.freq]==0)
				TPR2.sub.freq[i] <- mean(marker.status.freq[i,idx2.freq]==0)
				TNR1.sub.freq[i] <- mean(marker.status.freq[i,idx3.freq]==1)
				TNR2.sub.freq[i] <- mean(marker.status.freq[i,idx4.freq]==1)
			}				
	}
}



# create the result table
result <- matrix(NA, nrow=14, ncol=16)
rownames(result) <- c("No Marker Detected", "Overall Efficacy", "Overall Futility", "Overall Continue",
                      "Marker Detected", "M+ Eff/M- Eff", "M+ Eff/M- Cont", "M+ Eff/M- Fut", "M+ Cont/M- Eff",
                      "M+ Cont/M- Cont", "M+ Cont/M- Fut", "M+ Fut/M- Eff", "M+ Fut/M- Cont", "M+ Fut/M- Fut")
colnames(result) <- c("Interim Bayes", "Interim Freq", "Final Bayes", "Final Freq", "Overall Bayes", "Overall Freq", 
			    "Sample Size Bayes", "Sample Size Freq", 
			    "TPR1 Bayes", "TPR1 Freq", "TPR2 Bayes", "TPR2 Freq", "TNR1 Bayes", "TNR1 Freq", "TNR2 Bayes", "TNR2 Freq")



# Interim and final analysis decision rate

# no marker detected: efficacy, futility and continue
result[2,c(1,3)] <- apply(stop.eff.all.Bayes,2,mean)
result[3,1] <- apply(stop.fut.all.Bayes,2,mean)
result[3,3] <- mean(is.na(prev.Bayes) & apply(stop.eff.all.Bayes,1,sum)==0 & apply(stop.fut.all.Bayes,1,sum)==0)
result[4,1] <- apply(continue.all.Bayes,2,mean)

result[2,c(2,4)] <- apply(stop.eff.all.freq,2,mean)
result[3,2] <- apply(stop.fut.all.freq,2,mean)
result[3,4] <- mean(is.na(prev.freq) & apply(stop.eff.all.freq,1,sum)==0 & apply(stop.fut.all.freq,1,sum)==0)
result[4,2] <- apply(continue.all.freq,2,mean)


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

result[6,3] <- mean((continue.all.Bayes[,1]>0 | continue.pos.Bayes[,1]>0 | continue.neg.Bayes[,1]>0) & (!is.na(prev.Bayes))
			  & apply(stop.eff.pos.Bayes,1,sum)>0 & apply(stop.eff.neg.Bayes,1,sum)>0)
result[8,3] <- mean((continue.all.Bayes[,1]>0 | continue.pos.Bayes[,1]>0 | continue.neg.Bayes[,1]>0) & (!is.na(prev.Bayes))
			  & apply(stop.eff.pos.Bayes,1,sum)>0 & apply(stop.eff.neg.Bayes,1,sum)==0)
result[12,3] <- mean((continue.all.Bayes[,1]>0 | continue.pos.Bayes[,1]>0 | continue.neg.Bayes[,1]>0) & (!is.na(prev.Bayes))
			  & apply(stop.eff.pos.Bayes,1,sum)==0 & apply(stop.eff.neg.Bayes,1,sum)>0)
result[14,3] <- mean((continue.all.Bayes[,1]>0 | continue.pos.Bayes[,1]>0 | continue.neg.Bayes[,1]>0) & (!is.na(prev.Bayes))
			  & apply(stop.eff.pos.Bayes,1,sum)==0 & apply(stop.eff.neg.Bayes,1,sum)==0)


result[6,2]  <- mean(stop.eff.pos.freq[,1]>0 & stop.eff.neg.freq[,1]>0)
result[7,2]  <- mean(stop.eff.pos.freq[,1]>0 & continue.neg.freq[,1]>0)
result[8,2]  <- mean(stop.eff.pos.freq[,1]>0 & stop.fut.neg.freq[,1]>0)
result[9,2]  <- mean(continue.pos.freq[,1]>0 & stop.eff.neg.freq[,1]>0)
result[10,2] <- mean(continue.pos.freq[,1]>0 & continue.neg.freq[,1]>0)
result[11,2] <- mean(continue.pos.freq[,1]>0 & stop.fut.neg.freq[,1]>0)
result[12,2] <- mean(stop.fut.pos.freq[,1]>0 & stop.eff.neg.freq[,1]>0)
result[13,2] <- mean(stop.fut.pos.freq[,1]>0 & continue.neg.freq[,1]>0)
result[14,2] <- mean(stop.fut.pos.freq[,1]>0 & stop.fut.neg.freq[,1]>0)

result[6,4] <- mean((continue.all.freq[,1]>0 | continue.pos.freq[,1]>0 | continue.neg.freq[,1]>0) & (!is.na(prev.freq))
			  & apply(stop.eff.pos.freq,1,sum)>0 & apply(stop.eff.neg.freq,1,sum)>0)
result[8,4] <- mean((continue.all.freq[,1]>0 | continue.pos.freq[,1]>0 | continue.neg.freq[,1]>0) & (!is.na(prev.freq))
			  & apply(stop.eff.pos.freq,1,sum)>0 & apply(stop.eff.neg.freq,1,sum)==0)
result[12,4] <- mean((continue.all.freq[,1]>0 | continue.pos.freq[,1]>0 | continue.neg.freq[,1]>0) & (!is.na(prev.freq))
			  & apply(stop.eff.pos.freq,1,sum)==0 & apply(stop.eff.neg.freq,1,sum)>0)
result[14,4] <- mean((continue.all.freq[,1]>0 | continue.pos.freq[,1]>0 | continue.neg.freq[,1]>0) & (!is.na(prev.freq))
			  & apply(stop.eff.pos.freq,1,sum)==0 & apply(stop.eff.neg.freq,1,sum)==0)


# overall trial decision rate
result[2,5] <- mean(final.Bayes=="eff_all")
result[3,5] <- mean(final.Bayes=="eff_no")
result[8,5] <- mean(final.Bayes=="eff_pos")
result[12,5] <- mean(final.Bayes=="eff_neg")

result[2,6] <- mean(final.freq=="eff_all")
result[3,6] <- mean(final.freq=="eff_no")
result[8,6] <- mean(final.freq=="eff_pos")
result[12,6] <- mean(final.freq=="eff_neg")


# dichotomization rate
result[1,1] <- result[2,1]+result[3,1]+result[4,1]
result[1,3] <- result[2,3]+result[3,3]
result[5,1] <- sum(result[6:14,1])
result[5,3] <- sum(result[c(6,8,12,14),3])
result[1,5] <- mean(final.Bayes=="eff_all" | final.Bayes=="eff_no")
result[5,5] <- mean(final.Bayes=="eff_pos" | final.Bayes=="eff_neg")

result[1,2] <- result[2,2]+result[3,2]+result[4,2]
result[1,4] <- result[2,4]+result[3,4]
result[5,2] <- sum(result[6:14,2])
result[5,4] <- sum(result[c(6,8,12,14),4])
result[1,6] <- mean(final.freq=="eff_all" | final.freq=="eff_no")
result[5,6] <- mean(final.freq=="eff_pos" | final.freq=="eff_neg")



# sample size
result[1,7] <- round(mean(sampsize.Bayes[final.Bayes=="eff_all" | final.Bayes=="eff_no"]))
result[5,7] <- round(mean(sampsize.pos.Bayes[final.Bayes=="eff_pos" | final.Bayes=="eff_neg"]))
result[6,7] <- round(mean(sampsize.neg.Bayes[final.Bayes=="eff_pos" | final.Bayes=="eff_neg"]))

result[1,8] <- round(mean(sampsize.freq[final.freq=="eff_all" | final.freq=="eff_no"]))
result[5,8] <- round(mean(sampsize.pos.freq[final.freq=="eff_pos" | final.freq=="eff_neg"]))
result[6,8] <- round(mean(sampsize.neg.freq[final.freq=="eff_pos" | final.freq=="eff_neg"]))



# TPR and TNR
result[1,9] <- round(mean(TPR1.all.Bayes[!is.nan(TPR1.all.Bayes) & !is.na(TPR1.all.Bayes)]),3)
result[1,10] <- round(mean(TPR1.all.freq[!is.nan(TPR1.all.freq) & !is.na(TPR1.all.freq)]),3)
result[1,11] <- round(mean(TPR2.all.Bayes[!is.nan(TPR2.all.Bayes) & !is.na(TPR2.all.Bayes)]),3)
result[1,12] <- round(mean(TPR2.all.freq[!is.nan(TPR2.all.freq) & !is.na(TPR2.all.freq)]),3)

result[5,9] <- round(mean(TPR1.sub.Bayes[!is.nan(TPR1.sub.Bayes) & !is.na(TPR1.sub.Bayes)]),3)
result[5,10] <- round(mean(TPR1.sub.freq[!is.nan(TPR1.sub.freq) & !is.na(TPR1.sub.freq)]),3)
result[5,11] <- round(mean(TPR2.sub.Bayes[!is.nan(TPR2.sub.Bayes) & !is.na(TPR2.sub.Bayes)]),3)
result[5,12] <- round(mean(TPR2.sub.freq[!is.nan(TPR2.sub.freq) & !is.na(TPR2.sub.freq)]),3)

result[1,13] <- round(mean(TNR1.all.Bayes[!is.nan(TNR1.all.Bayes) & !is.na(TNR1.all.Bayes)]),3)
result[1,14] <- round(mean(TNR1.all.freq[!is.nan(TNR1.all.freq) & !is.na(TNR1.all.freq)]),3)
result[1,15] <- round(mean(TNR2.all.Bayes[!is.nan(TNR2.all.Bayes) & !is.na(TNR2.all.Bayes)]),3)
result[1,16] <- round(mean(TNR2.all.freq[!is.nan(TNR2.all.freq) & !is.na(TNR2.all.freq)]),3)

result[5,13] <- round(mean(TNR1.sub.Bayes[!is.nan(TNR1.sub.Bayes) & !is.na(TNR1.sub.Bayes)]),3)
result[5,14] <- round(mean(TNR1.sub.freq[!is.nan(TNR1.sub.freq) & !is.na(TNR1.sub.freq)]),3)
result[5,15] <- round(mean(TNR2.sub.Bayes[!is.nan(TNR2.sub.Bayes) & !is.na(TNR2.sub.Bayes)]),3)
result[5,16] <- round(mean(TNR2.sub.freq[!is.nan(TNR2.sub.freq) & !is.na(TNR2.sub.freq)]),3)




# output results
result2 <- xtable(result, digits=c(4,3,3,3,3,3,3,0,0,3,3,3,3,3,3,3,3))
print(result2, include.rownames = TRUE, booktabs = TRUE)

\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrrrrrrrrrrrr}
  \toprule
 & Interim Bayes & Interim Freq & Final Bayes & Final Freq & Overall Bayes & Overall Freq & Sample Size Bayes & Sample Size Freq & TPR1 Bayes & TPR1 Freq & TPR2 Bayes & TPR2 Freq & TNR1 Bayes & TNR1 Freq & TNR2 Bayes & TNR2 Freq \\ 
  \midrule
No Marker Detected & 0.019 & 0.037 & 0.000 & 0.000 & 0.141 & 0.156 & 399 & 390 & 0.993 & 1.000 & 0.993 & 1.000 & 0.007 & 0.000 & 0.007 & 0.000 \\ 
  Overall Efficacy & 0.014 & 0.034 & 0.000 & 0.000 & 0.140 & 0.156 &  &  &  &  &  &  &  &  &  &  \\ 
  Overall Futility & 0.001 & 0.000 & 0.000 & 0.000 & 0.001 & 0.000 &  &  &  &  &  &  &  &  &  &  \\ 
  Overall Continue & 0.004 & 0.003 &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\ 
  Marker Detected & 0.981 & 0.963 & 0.502 & 0.296 & 0.859 & 0.844 & 196 & 186 & 0.959 & 0.959 & 0.959 & 0.959 & 0.856 & 0.908 & 0.856 & 0.908 \\ 
  M+ Eff/M- Eff & 0.046 & 0.063 & 0.080 & 0.059 &  &  & 189 & 189 &  &  &  &  &  &  &  &  \\ 
  M+ Eff/M- Cont & 0.497 & 0.292 &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\ 
  M+ Eff/M- Fut & 0.437 & 0.607 & 0.422 & 0.237 & 0.859 & 0.844 &  &  &  &  &  &  &  &  &  &  \\ 
  M+ Cont/M- Eff & 0.000 & 0.000 &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\ 
  M+ Cont/M- Cont & 0.001 & 0.000 &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\ 
  M+ Cont/M- Fut & 0.000 & 0.001 &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\ 
  M+ Fut/M- Eff & 0.000 & 0.000 & 0.000 & 0.000 & 0.000 & 0.000 &  &  &  &  &  &  &  &  &  &  \\ 
  M+ Fut/M- Cont & 0.000 & 0.000 &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\ 
  M+ Fut/M- Fut & 0.000 & 0.000 & 0.000 & 0.000 &  &  &  &  &  &  &  &  &  &  &  &  \\ 
   \bottomrule
\end{tabular}
\end{table}


