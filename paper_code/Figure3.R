# Direct to the folder where the simulation summaries are saved for all scenarios
setwd("./figures")


# Figure 3
designs <- c("BCM", "FDM")


png(filename = "Figure3a.png", units="in", width=10, height=15, res=300)

par(mfrow=c(4,1), mar=c(4.1, 4.1, 3.1, 14.1), xpd=TRUE)


# setting 1
load("setting1.RData")

detection.rate.v2 <- detection.rate[c(1,2),]

colors2 <- c("paleturquoise", "lightsalmon")

barplot(detection.rate.v2, main = expression(paste("1. No Marker Effect, No Treatment Effect")), names.arg = designs, space=c(0,0.75),
        xlab = "", ylab = "", col = colors2, border="black")

coord <- par("usr")

segments(x0=0,y0=0, x1=coord[2]*0.8, y1=0)

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1])

text(mean(c(coord[1],coord[2])), coord[4]/2, expression(paste("Treatment Effect: ", italic("HR"), "" > "0.95")), cex=1) 

legend("topright", inset=c(-0.18,0), legend=c("True Negatives","False Positives"), bty="n", 
       pch=15, col=colors2, cex=1.2)



# setting 2
load("setting2.RData")

detection.rate.v2 <- detection.rate[c(5,4),]

colors2 <- c("darkturquoise", "darkorange2")

barplot(detection.rate.v2, main = expression(paste("2. No Marker Effect, Constant Treatment Effect with Effect Size ", Delta, "=0.5")), 
        names.arg = designs, space=c(0,0.75), xlab = "", ylab = "", col = colors2, border="black")

coord <- par("usr")

segments(x0=0,y0=0, x1=coord[2]*0.8, y1=0)

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1])

text(mean(c(coord[1],coord[2])), coord[4]/2, expression(paste("Treatment Effect: ", italic("HR"), "" <= "0.8")), cex=1) 

legend("topright", inset=c(-0.18,0), legend=c("True Positives","False Negatives"), bty="n", 
       pch=15, col=colors2, cex=1.2)



# setting 3
load("setting3.RData")

detection.rate.v2 <- detection.rate[c(1,2),]

colors2 <- c("paleturquoise", "lightsalmon")

barplot(detection.rate.v2, main = expression(paste("3. Prognostic Marker Effect, No Treatment Effect")), names.arg = designs, space=c(0,0.75),
        xlab = "", ylab = "", col = colors2, border="black")

coord <- par("usr")

segments(x0=0,y0=0, x1=coord[2]*0.8, y1=0)

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1])

text(mean(c(coord[1],coord[2])), coord[4]/2, expression(paste("Treatment Effect: ", italic("HR"), "" > "0.95")), cex=1) 

legend("topright", inset=c(-0.18,0), legend=c("True Negatives","False Positives"), bty="n", 
       pch=15, col=colors2, cex=1.2)



# setting 4
load("setting4.RData")

detection.rate.v2 <- detection.rate[c(5,4,1,2),]

colors2 <- c("darkturquoise", "darkorange2", "paleturquoise", "lightsalmon")

barplot(detection.rate.v2, main = expression(paste("4. Predictive Marker Effect (Perfectly Dichotomous) with Max Effect Size ", Delta, "=0.4")), 
        names.arg = designs, space=c(0,0.75), xlab = "", ylab = "", col = colors2, border="black")

coord <- par("usr")

segments(x0=0,y0=0, x1=coord[2]*0.8, y1=0)

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1])

segments(x0=0,y0=coord[4], x1=coord[2]*0.8, y1=coord[4])

text(mean(c(coord[1],coord[2])), detection.rate.v2[1,1]/2+detection.rate.v2[2,1]/2, expression(paste("Treatment Effect: ", italic("HR"), "" <= "0.8")), cex=1) 

text(mean(c(coord[1],coord[2])), coord[4]-detection.rate.v2[3,1]/2-detection.rate.v2[4,1]/2, expression(paste("Treatment Effect: ", italic("HR"), "" > "0.95")), cex=1) 

legend("topright", inset=c(-0.18,0), legend=c("True Positives", "False Negatives", "True Negatives","False Positives"), bty="n", 
       pch=15, col=colors2, cex=1.2)

dev.off()



png(filename = "Figure3b.png", units="in", width=10, height=15, res=300)

par(mfrow=c(4,1), mar=c(4.1, 4.1, 3.1, 14.1), xpd=TRUE)


# setting 5
load("setting5.RData")

detection.rate.v2 <- detection.rate[c(5,4,3,1,2),]

colors2 <- c("darkturquoise", "darkorange2", "black", "paleturquoise", "lightsalmon")

barplot(detection.rate.v2, main = expression(paste("5. Predictive Marker Effect (Nearly Dichotomous) with Max Effect Size ", Delta, "=0.4")), 
        names.arg = designs, space=c(0,0.75), xlab = "", ylab = "", col = colors2, border="black")

segments(x0=0,y0=0, x1=coord[2]*0.8, y1=0)

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1])

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1])

segments(x0=0,y0=coord[4], x1=coord[2]*0.8, y1=coord[4])

text(mean(c(coord[1],coord[2])), detection.rate.v2[1,1]/2+detection.rate.v2[2,1]/2, expression(paste("Treatment Effect: ", italic("HR"), "" <= "0.8")), cex=1) 

text(mean(c(coord[1],coord[2])), detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1]/2, expression(paste("Treatment Effect: ", "0.8" < "", italic("HR"), "" <= "0.95")), cex=1) 

text(mean(c(coord[1],coord[2])), coord[4]-detection.rate.v2[4,1]/2-detection.rate.v2[5,1]/2, expression(paste("Treatment Effect: ", italic("HR"), "" > "0.95")), cex=1) 

legend("topright", inset=c(-0.18,0), legend=c("True Positives", "False Negatives", "Boundary Trt Effect", "True Negatives","False Positives"), bty="n", 
       pch=15, col=colors2, cex=1.2)



# setting 6
load("setting6.RData")

detection.rate.v2 <- detection.rate[c(5,4,3,1,2),]

colors2 <- c("darkturquoise", "darkorange2", "black", "paleturquoise", "lightsalmon")

barplot(detection.rate.v2, main = expression(paste("6. Predictive Marker Effect (Linear) with Max Effect Size ", Delta, "=0.4")), 
        names.arg = designs, space=c(0,0.75), xlab = "", ylab = "", col = colors2, border="black")

segments(x0=0,y0=0, x1=coord[2]*0.8, y1=0)

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1])

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1])

segments(x0=0,y0=coord[4], x1=coord[2]*0.8, y1=coord[4])

text(mean(c(coord[1],coord[2])), detection.rate.v2[1,1]/2+detection.rate.v2[2,1]/2, expression(paste("Treatment Effect: ", italic("HR"), "" <= "0.8")), cex=1) 

text(mean(c(coord[1],coord[2])), detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1]/2, expression(paste("Treatment Effect: ", "0.8" < "", italic("HR"), "" <= "0.95")), cex=1) 

text(mean(c(coord[1],coord[2])), coord[4]-detection.rate.v2[4,1]/2-detection.rate.v2[5,1]/2, expression(paste("Treatment Effect: ", italic("HR"), "" > "0.95")), cex=1) 

legend("topright", inset=c(-0.18,0), legend=c("True Positives", "False Negatives", "Boundary Trt Effect", "True Negatives","False Positives"), bty="n", 
       pch=15, col=colors2, cex=1.2)



# setting 7
load("setting7.RData")

detection.rate.v2 <- detection.rate[c(5,4,3),]

colors2 <- c("darkturquoise", "darkorange2", "black")

barplot(detection.rate.v2, main = expression(paste("7. Predictive Marker Effect (Non-linear and Monotone) with Max Effect Size ", Delta, "=0.4")),
        names.arg = designs, space=c(0,0.75), xlab = "", ylab = "", col = colors2, border="black")

segments(x0=0,y0=0, x1=coord[2]*0.8, y1=0)

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1])

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1])

text(mean(c(coord[1],coord[2])), detection.rate.v2[1,1]/2+detection.rate.v2[2,1]/2, expression(paste("Treatment Effect: ", italic("HR"), "" <= "0.8")), cex=1) 

text(mean(c(coord[1],coord[2])), detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1]/2, expression(paste("Treatment Effect: ", "0.8" < "", italic("HR"), "" <= "0.95")), cex=1) 

legend("topright", inset=c(-0.18,0), legend=c("True Positives", "False Negatives", "Boundary Trt Effect"), bty="n", 
       pch=15, col=colors2, cex=1.2)



# setting 8
load("setting8.RData")

detection.rate.v2 <- detection.rate[c(5,4,3,1,2),]

colors2 <- c("darkturquoise", "darkorange2", "black", "paleturquoise", "lightsalmon")

barplot(detection.rate.v2, main = expression(paste("8. Predictive Marker Effect (Non-linear and Monotone) with Max Effect Size ", Delta, "=0.4")),  
        names.arg = designs, space=c(0,0.75), xlab = "", ylab = "", col = colors2, border="black")

coord <- par("usr")

segments(x0=0,y0=0, x1=coord[2]*0.8, y1=0)

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1])

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1])

segments(x0=0,y0=coord[4], x1=coord[2]*0.8, y1=coord[4])

text(mean(c(coord[1],coord[2])), detection.rate.v2[1,1]/2+detection.rate.v2[2,1]/2, expression(paste("Treatment Effect: ", italic("HR"), "" <= "0.8")), cex=1) 

text(mean(c(coord[1],coord[2])), detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1]/2, expression(paste("Treatment Effect: ", "0.8" < "", italic("HR"), "" <= "0.95")), cex=1) 

text(mean(c(coord[1],coord[2])), coord[4]-detection.rate.v2[4,1]/2-detection.rate.v2[5,1]/2, expression(paste("Treatment Effect: ", italic("HR"), "" > "0.95")), cex=1) 

legend("topright", inset=c(-0.18,0), legend=c("True Positives", "False Negatives", "Boundary Trt Effect", "True Negatives","False Positives"), bty="n", 
       pch=15, col=colors2, cex=1.2)


dev.off()



png(filename = "Figure3c.png", units="in", width=10, height=15, res=300)

par(mfrow=c(4,1), mar=c(4.1, 4.1, 3.1, 14.1), xpd=TRUE)


# setting 9
load("setting9.RData")

detection.rate.v2 <- detection.rate[c(5,4,3,1,2),]

colors2 <- c("darkturquoise", "darkorange2", "black", "paleturquoise", "lightsalmon")

barplot(detection.rate.v2,  main = expression(paste("9. Predictive Marker Effect (Non-linear and Non-monotone) with Max Effect Size ", Delta, "=0.4")),  
        names.arg = designs, space=c(0,0.75), xlab = "", ylab = "", col = colors2, border="black")

coord <- par("usr")

segments(x0=0,y0=0, x1=coord[2]*0.8, y1=0)

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1])

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1])

segments(x0=0,y0=coord[4], x1=coord[2]*0.8, y1=coord[4])

text(mean(c(coord[1],coord[2])), detection.rate.v2[1,1]/2+detection.rate.v2[2,1]/2, expression(paste("Treatment Effect: ", italic("HR"), "" <= "0.8")), cex=1) 

text(mean(c(coord[1],coord[2])), detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1]/2, expression(paste("Treatment Effect: ", "0.8" < "", italic("HR"), "" <= "0.95")), cex=1) 

text(mean(c(coord[1],coord[2])), coord[4]-detection.rate.v2[4,1]/2-detection.rate.v2[5,1]/2, expression(paste("Treatment Effect: ", italic("HR"), "" > "0.95")), cex=1) 

legend("topright", inset=c(-0.18,0), legend=c("True Positives", "False Negatives", "Boundary Trt Effect", "True Negatives","False Positives"), bty="n", 
       pch=15, col=colors2, cex=1.2)



# setting 10
load("setting10.RData")

detection.rate.v2 <- detection.rate[c(5,4,3,1,2),]

colors2 <- c("darkturquoise", "darkorange2", "black", "paleturquoise", "lightsalmon")

barplot(detection.rate.v2,  main = expression(paste("10. Predictive Marker Effect (Non-linear and Non-monotone) with Max Effect Size ", Delta, "=0.4")),  
        names.arg = designs, space=c(0,0.75), xlab = "", ylab = "", col = colors2, border="black")

coord <- par("usr")

segments(x0=0,y0=0, x1=coord[2]*0.8, y1=0)

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1])

segments(x0=0,y0=detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1], x1=coord[2]*0.8, y1=detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1])

segments(x0=0,y0=coord[4], x1=coord[2]*0.8, y1=coord[4])

text(mean(c(coord[1],coord[2])), detection.rate.v2[1,1]/2+detection.rate.v2[2,1]/2, expression(paste("Treatment Effect: ", italic("HR"), "" <= "0.8")), cex=1) 

text(mean(c(coord[1],coord[2])), detection.rate.v2[1,1]+detection.rate.v2[2,1]+detection.rate.v2[3,1]/2, expression(paste("Treatment Effect: ", "0.8" < "", italic("HR"), "" <= "0.95")), cex=1) 

text(mean(c(coord[1],coord[2])), coord[4]-detection.rate.v2[4,1]/2-detection.rate.v2[5,1]/2, expression(paste("Treatment Effect: ", italic("HR"), "" > "0.95")), cex=1) 

legend("topright", inset=c(-0.18,0), legend=c("True Positives", "False Negatives", "Boundary Trt Effect", "True Negatives","False Positives"), bty="n", 
       pch=15, col=colors2, cex=1.2)


dev.off()
