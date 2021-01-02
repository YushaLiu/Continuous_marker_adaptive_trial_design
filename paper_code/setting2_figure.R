# Direct to the folder where the output is saved
setwd("./output/setting2")


# n = 500, threshold for marker dichotomization = 0.1, one interim check
load("trial_results.RData")

fail.MCMC <- res$fail.MCMC
final.Bayes <- res$final.Bayes
marker.value.Bayes <- res$marker.value.Bayes
marker.status.Bayes <- res$marker.status.Bayes

final.freq <- res$final.freq
marker.value.freq <- res$marker.value.freq
marker.status.freq <- res$marker.status.freq


sum(fail.MCMC)
[1] 1


# delete the iteration of Bayesian design where MCMC fails 
idx <- which(fail.MCMC==0)

final.Bayes <- final.Bayes[idx]
marker.value.Bayes <- marker.value.Bayes[idx,] 
marker.status.Bayes <- marker.status.Bayes[idx,]


# overall trial decision rate
decision.rate <- matrix(NA, nrow=3, ncol=2)
rownames(decision.rate) <- c("fut_all", "eff_sub", "eff_all")
colnames(decision.rate) <- c("Bayes", "freq")

decision.rate[1,1] <- mean(final.Bayes=="eff_no")	# Bayes overall futility
decision.rate[2,1] <- mean(final.Bayes=="eff_pos" | final.Bayes=="eff_neg")	# Bayes subgroup efficacy
decision.rate[3,1] <- mean(final.Bayes=="eff_all")	# Bayes overall efficacy
decision.rate[1,2] <- mean(final.freq=="eff_no")	# Freq overall futility
decision.rate[2,2] <- mean(final.freq=="eff_pos" | final.freq=="eff_neg")	# Freq subgroup efficacy
decision.rate[3,2] <- mean(final.freq=="eff_all")	# Freq overall efficacy



# true benefit and lack of benefit detection rate

# trt effect as a function of continuous marker
fnc.trt <- function(x) {return(log(0.5))}


# true positive rate and true negative rate per trial iteration
TPR.Bayes <- rep(NA, length(final.Bayes))

TPR.freq <- rep(NA, length(final.freq))


for(i in 1:length(final.Bayes)) {

	if(final.Bayes[i]=="eff_all") {
		TPR.Bayes[i] <- 1
	}

	else if(final.Bayes[i]=="eff_no") {
		TPR.Bayes[i] <- 0
	}

	else {
		x.Bayes <- marker.value.Bayes[i,]
		idx.Bayes <- which( unlist(lapply(x.Bayes,fnc.trt)) <= log(0.8) )		

			if(final.Bayes[i]=="eff_pos") {
				TPR.Bayes[i] <- mean(marker.status.Bayes[i,idx.Bayes]==1)
			}

			else if (final.Bayes[i]=="eff_neg") {
				TPR.Bayes[i] <- mean(marker.status.Bayes[i,idx.Bayes]==0)
			}				
	}
}


for(i in 1:length(final.freq)) {

	if(final.freq[i]=="eff_all") {
		TPR.freq[i] <- 1
	}

	else if(final.freq[i]=="eff_no") {
		TPR.freq[i] <- 0
	}

	else {
		x.freq <- marker.value.freq[i,]
		idx.freq <- which( unlist(lapply(x.freq,fnc.trt)) <= log(0.8) )

			if(final.freq[i]=="eff_pos") {
				TPR.freq[i] <- mean(marker.status.freq[i,idx.freq]==1)
			}

			else if (final.freq[i]=="eff_neg") {
				TPR.freq[i] <- mean(marker.status.freq[i,idx.freq]==0)
			}				
	}
}


# average detection rate multiplied by the corresponding marker prevalence
detection.rate <- matrix(NA, nrow=5, ncol=2)
rownames(detection.rate) <- c("TN", "FP", "NA", "FN", "TP")
colnames(detection.rate) <- c("Bayes", "freq")

detection.rate[1,1] <- 0	# Bayes true negatives
detection.rate[2,1] <- 0	# Bayes false positives
detection.rate[3,1] <- 0	# prevalence of marker with trt effect on the boundary
detection.rate[4,1] <- 1-mean(TPR.Bayes[!is.nan(TPR.Bayes) & !is.na(TPR.Bayes)])	# Bayes false negatives
detection.rate[5,1] <- mean(TPR.Bayes[!is.nan(TPR.Bayes) & !is.na(TPR.Bayes)])	# Bayes true positives


detection.rate[1,2] <- 0	# Freq true negatives
detection.rate[2,2] <- 0	# Freq false positives
detection.rate[3,2] <- 0	# prevalence of marker with trt effect on the boundary
detection.rate[4,2] <- 1-mean(TPR.freq[!is.nan(TPR.freq) & !is.na(TPR.freq)])	# Freq false negatives
detection.rate[5,2] <- mean(TPR.freq[!is.nan(TPR.freq) & !is.na(TPR.freq)])	# Freq true positives


designs <- c("BCM", "FDM")
colors1 <- c("yellow", "orange", "turquoise")	# color code depends on closeness to the truth, which might change across scenarios
colors2 <- c("paleturquoise", "lightsalmon", "lightgray", "darkorange2", "darkturquoise")	# color code is consistent across scenarios

rDataPath <- "./figures/"

save(decision.rate, detection.rate, colors1, colors2, file=file.path(rDataPath,"setting2.RData"))












