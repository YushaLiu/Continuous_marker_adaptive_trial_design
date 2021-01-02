# Direct to the folder where the simulation summaries are saved for all scenarios
setwd("./figures")


# Figure 2
designs <- c("BCM", "FDM")


png(filename = "Figure2a.png", units="in", width=9, height=15, res=300)

par(mfrow=c(4,1), mar=c(4.1, 4.1, 3.1, 14.1), xpd=TRUE)


# setting1
load("setting1.RData")

colors1 <- c("turquoise", "yellow", "orange")

barplot(decision.rate, main = expression(paste("1. No Marker Effect, No Treatment Effect")), names.arg = designs,
        xlab = "", ylab = "", col = colors1, border="black")

legend("topright", inset=c(-0.22,0), legend=c("Overall Futility","Subgroup Efficacy","Overall Efficacy"), bty="n", 
       pch=15, col=colors1, cex=1.2)


# setting 2
load("setting2.RData")

decision.rate.v2 <- decision.rate[c(3,2),]

colors1 <- c("turquoise", "yellow")

barplot(decision.rate.v2, main = expression(paste("2. No Marker Effect, Constant Treatment Effect with Effect Size ", Delta, "=0.5")), names.arg = designs,
        xlab = "", ylab = "", col = colors1, border="black")

legend("topright", inset=c(-0.22,0), legend=c("Overall Efficacy","Subgroup Efficacy"), bty="n", 
       pch=15, col=colors1, cex=1.2)


# setting 3
load("setting3.RData")

colors1 <- c("turquoise", "yellow", "orange")

barplot(decision.rate, main = expression(paste("3. Prognostic Marker Effect, No Treatment Effect")), names.arg = designs,
        xlab = "", ylab = "", col = colors1, border="black")

legend("topright", inset=c(-0.22,0), legend=c("Overall Futility","Subgroup Efficacy","Overall Efficacy"), bty="n", 
       pch=15, col=colors1, cex=1.2)


# setting 4
load("setting4.RData")

decision.rate.v2 <- decision.rate[c(2,3),]

colors1 <- c("turquoise", "yellow")

barplot(decision.rate.v2, main = expression(paste("4. Predictive Marker Effect (Perfectly Dichotomous) with Max Effect Size ", Delta, "=0.4")), names.arg = designs,
        xlab = "", ylab = "", col = colors1, border="black")

legend("topright", inset=c(-0.22,0), legend=c("Subgroup Efficacy","Overall Efficacy"), bty="n", 
       pch=15, col=colors1, cex=1.2)

dev.off()



png(filename = "Figure2b.png", units="in", width=9, height=15, res=300)

par(mfrow=c(4,1), mar=c(4.1, 4.1, 3.1, 14.1), xpd=TRUE)


# setting 5
load("setting5.RData")

decision.rate.v2 <- decision.rate[c(2,3),]

colors1 <- c("turquoise", "yellow")

barplot(decision.rate.v2, main = expression(paste("5. Predictive Marker Effect (Nearly Dichotomous) with Max Effect Size ", Delta, "=0.4")), names.arg = designs,
        xlab = "", ylab = "", col = colors1, border="black")

legend("topright", inset=c(-0.22,0), legend=c("Subgroup Efficacy","Overall Efficacy"), bty="n", 
       pch=15, col=colors1, cex=1.2)


# setting 6
load("setting6.RData")

decision.rate.v2 <- decision.rate[c(2,3),]

colors1 <- c("turquoise", "yellow")

barplot(decision.rate.v2, main = expression(paste("6. Predictive Marker Effect (Linear) with Max Effect Size ", Delta, "=0.4")), names.arg = designs,
        xlab = "", ylab = "", col = colors1, border="black")

legend("topright", inset=c(-0.22,0), legend=c("Subgroup Efficacy","Overall Efficacy"), bty="n", 
       pch=15, col=colors1, cex=1.2)


# setting 7
load("setting7.RData")

decision.rate.v2 <- decision.rate[c(3,2),]

colors1 <- c("turquoise", "yellow")

barplot(decision.rate.v2, main = expression(paste("7. Predictive Marker Effect (Non-linear and Monotone) with Max Effect Size ", Delta, "=0.4")), names.arg = designs,
        xlab = "", ylab = "", col = colors1, border="black")

legend("topright", inset=c(-0.22,0), legend=c("Overall Efficacy","Subgroup Efficacy"), bty="n", 
       pch=15, col=colors1, cex=1.2)


# setting 8
load("setting8.RData")

decision.rate.v2 <- decision.rate[c(2,3,1),]

colors1 <- c("turquoise", "yellow", "orange")

barplot(decision.rate.v2, main = expression(paste("8. Predictive Marker Effect (Non-linear and Monotone) with Max Effect Size ", Delta, "=0.4")), names.arg = designs,
        xlab = "", ylab = "", col = colors1, border="black")

legend("topright", inset=c(-0.22,0), legend=c("Subgroup Efficacy","Overall Efficacy","Overall Futility"), bty="n", 
       pch=15, col=colors1, cex=1.2)


dev.off()



png(filename = "Figure2c.png", units="in", width=9, height=15, res=300)

par(mfrow=c(4,1), mar=c(4.1, 4.1, 3.1, 14.1), xpd=TRUE)


# setting 9
load("setting9.RData")

decision.rate.v2 <- decision.rate[c(2,3,1),]

colors1 <- c("turquoise", "yellow", "orange")

barplot(decision.rate.v2, main = expression(paste("9. Predictive Marker Effect (Non-linear and Non-monotone) with Max Effect Size ", Delta, "=0.4")), names.arg = designs,
        xlab = "", ylab = "", col = colors1, border="black")

legend("topright", inset=c(-0.22,0), legend=c("Subgroup Efficacy","Overall Efficacy","Overall Futility"), bty="n", 
       pch=15, col=colors1, cex=1.2)


# setting 10
load("setting10.RData")

decision.rate.v2 <- decision.rate[c(2,3),]

colors1 <- c("turquoise", "yellow")

barplot(decision.rate.v2, main = expression(paste("10. Predictive Marker Effect (Non-linear and Non-monotone) with Max Effect Size ", Delta, "=0.4")), names.arg = designs,
        xlab = "", ylab = "", col = colors1, border="black")

legend("topright", inset=c(-0.22,0), legend=c("Subgroup Efficacy","Overall Efficacy"), bty="n", 
       pch=15, col=colors1, cex=1.2)


dev.off()


