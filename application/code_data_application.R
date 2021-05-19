setwd("./application/")


####################### read in the trial data to run BCM design with age and WBC as biomarkers ###########################
# load in trial data
data <- read.csv("data.csv", header=TRUE)

# check for missing values in the biomarkers
na.idx <- is.na(data$age) | is.na(data$wbc)
sum(na.idx)

# remove the observations with missing biomarker values
data <- data[!na.idx,]

# map age values to their percentiles in the dataset
age.ecdf <- ecdf(data$age)
age2 <- rep(NA, nrow(data))

for(i in 1:length(age2)) {
  age2[i] <- age.ecdf(data$age[i])
}

# map WBC values to their percentiles in the dataset
wbc.ecdf <- ecdf(data$wbc)
wbc2 <- rep(NA, nrow(data))

for(i in 1:length(wbc2)) {
  wbc2[i] <- wbc.ecdf(data$wbc[i])
}

# create the matrix with observed patient data
surdata <- data.frame(trt=data$trt, enter=data$Randomization_Day, status=data$ccr_x2, T_date=data$Event_Censor_Day)



####################### implement the Bayesian continuous marker design with single biomarker ###########################
# load in the function to run BCM design
source("./BCM_data_application.R")

# BCM design with age as the only biomarker
res.age <- BCM_data_application(surdata = surdata, X=data.frame(x=age2), n_events=sum(surdata$status), m.prev.lwr = 0.25, check = c(0.75, 1.00), 
                                maxIntKnots = 9, censor = FALSE, MCMC.specs = list(burnin=30000, B=2000, thin=30, multiplier=1/3), 
                                alpha = 0.1, p_eff_all = rep(0.95,2), p_eff_grp = rep(0.95,2), p_fut_all = 0.05, p_fut_grp = 0.05, stop_eff = TRUE)


# BCM design with WBC as the only biomarker
res.wbc <- BCM_data_application(surdata = surdata, X=data.frame(x=wbc2), n_events=sum(surdata$status), m.prev.lwr = 0.25, check = c(0.75, 1.00), 
                                maxIntKnots = 9, censor = FALSE, MCMC.specs = list(burnin=30000, B=2000, thin=30, multiplier=1/3), 
                                alpha = 0.1, p_eff_all = rep(0.95,2), p_eff_grp = rep(0.95,2), p_fut_all = 0.05, p_fut_grp = 0.05, stop_eff = TRUE)



####################### implement the Bayesian continuous marker design with both biomarkers ###########################
res.both <- BCM_data_application(surdata = surdata, X=data.frame(x1=age2, x2=wbc2), n_events=sum(surdata$status), 
                                 m.prev.lwr = 0.25,  check = c(0.75, 1.00), 
                                 maxIntKnots = 5, censor = FALSE, MCMC.specs = list(burnin=60000, B=2000, thin=120, multiplier=1/5), 
                                 alpha = 0.1, p_eff_all = rep(0.95,2), p_eff_grp = rep(0.95,2), p_fut_all = 0.05,	p_fut_grp = 0.05, stop_eff = TRUE)






