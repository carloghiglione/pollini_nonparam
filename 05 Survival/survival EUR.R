
library(survival)
library(survminer)
library(dplyr) 
library(ggplot2)
library(knitr)
library(broom)
setwd("D:/Marta/Politecnico/Non parametric statistics/Pollini project")
# Take just EUR countries
survival_data <- read.csv("D:/Marta/Politecnico/Non parametric statistics/Pollini project/survival_data.csv", sep="")
country_list <- read.csv("D:/Marta/Politecnico/Non parametric statistics/Pollini project/country_list.csv")
Eco_data_29_10_21 <- read.csv("D:/Marta/Politecnico/Non parametric statistics/Pollini project/Eco_data_29_10_21.csv")
colnames(Eco_data_29_10_21)[1]<-"country"
library(readxl)
Trips_labeled <- read_excel("Trips_labeled.xlsx")
Trips_labeled<-data.frame(country=Trips_labeled$Country,mobility=Trips_labeled$Mob_lab, EU=Trips_labeled$EU_lab)
countries<-country_list$Country
survival_data_EUR<-filter(survival_data, country==countries)
#uni_scores<-data.frame(country=Eco_data_29_10_21$Country,uni_score_norm=Eco_data_29_10_21$uni_score_norm,uni_score=Eco_data_29_10_21$uni_score)

# data with country,ID,Time(years),censored(T/F),all covariates 
data<-merge(survival_data_EUR,Eco_data_29_10_21,by="country")
data<-merge(data,Trips_labeled)
head(data)
# levels(factor(data$Cens)) # F T
data$status_fact=factor(data$Cens, labels = (c('Event', 'Censor')))
head(data)

ggplot(data=head(data),aes(x=ID,y=Time)) + 
  geom_bar(stat='identity',width=0.2) +
  geom_point(aes(color=status_fact,shape=status_fact),size=6) +
  coord_flip()

# head(Surv(data$Time, data$Cens==FALSE))
Surv(data$Time, data$Cens==FALSE)

# KM est
fit <- survfit(Surv(Time, Cens==FALSE) ~ 1, data = data)
fit
summary(fit)
surv_median(fit) # median 16: after 16 years prob of not moving is 50%

plot(fit, conf.int = T, xlab='Time [Years]', ylab = 'No going out Probability', col='red',
     main="Kaplan-Meier Curve for Non outgoing people")
ggsurvplot(fit,
           data=data,
           title="Kaplan-Meier Curve for outgoing prople"
           )

# cumulative failure prob
cumulative_incidence <- 1 - fit$surv

ggsurvplot(fit,
           data=data,
           title="Cumulative Incidence Curve for outgoing people",
           fun='event')

# hazard
H <- fit$cumhaz

ggsurvplot(fit,
           data=data,
           fun='cumhaz',
           title="Cumulative Hazard Curve for outgoing people")

#########################################################################
#GROUPING

### high vs low mobility

fit_mob <- survfit(Surv(Time, Cens==FALSE) ~ mobility, data=data)
# levels(factor(data$mobility))
ggsurvplot(fit_mob, conf.int = T,
           data=data,
           risk.table = TRUE, 
           surv.median.line = "hv", # Specify median survival
           
           legend.labs=c("high","low"), legend.title="mobility",  
           palette=c("darkblue","cyan3"), 
           title="Kaplan-Meier Curves by mobility for outgoing people")
log_rank_test <- survdiff(Surv(Time, Cens==FALSE) ~ mobility, data=data)
log_rank_test
# 0.9 no evidence for difference 

###country 

fit_c <- survfit(Surv(Time, Cens==FALSE) ~ country, data=data)
# levels(factor(data$mobility))
ggsurvplot(fit_c, conf.int = T,
           data=data,
           risk.table = TRUE, 
           surv.median.line = "hv", # Specify median survival
           
           
           title="Kaplan-Meier Curves by country for outgoing people")
log_rank_test <- survdiff(Surv(Time, Cens==FALSE) ~ country, data=data)
log_rank_test
# 1e-05 difference at least in two countries 

###country EU vs non-EU 

levels(factor(data$EU))
fit_eu <- survfit(Surv(Time, Cens==FALSE) ~ EU, data=data)
# levels(factor(data$mobility))
ggsurvplot(fit_eu, conf.int = T,
           data=data,
           risk.table = TRUE, 
           surv.median.line = "hv", # Specify median survival
           legend.labs=c("N","Y"), legend.title="EU",  
           palette=c("darkblue","cyan3"), 
           title="Kaplan-Meier Curves by country for outgoing people")

log_rank_test <- survdiff(Surv(Time, Cens==FALSE) ~ EU, data=data)
log_rank_test
# pval 0.4 no diff

###uniscore<20 vs uniscore>=20

range(data$uni_score)
plot(data$uni_score)
data$uni_score_20 <- cut(data$uni_score, breaks=c(0, 20, Inf), labels=c("low", "high"))
levels(factor(data$uni_score_20))

fit_uni_score_20 <- survfit(Surv(Time, Cens==FALSE) ~ uni_score_20, data=data)

ggsurvplot(fit_uni_score_20, conf.int = T,
           data=data,
           risk.table = TRUE, 
           surv.median.line = "hv", # Specify median survival
           legend.labs=c("low","high"), legend.title="Uni_score",  
           palette=c("darkblue","cyan3"), 
           title="Kaplan-Meier Curves by country for outgoing people")

log_rank_test <- survdiff(Surv(Time, Cens==FALSE) ~ uni_score_20, data=data)
log_rank_test
# 0.3 no evidence 

###uniscorenorm

range(data$uni_score_norm)
plot(data$uni_score_norm)
data$uni_score_norm_15 <- cut(data$uni_score, breaks=c(0, 0.15, Inf), labels=c("low", "high"))
levels(factor(data$uni_score_20))

fit_uni_score_norm_15 <- survfit(Surv(Time, Cens==FALSE) ~ uni_score_norm_15, data=data)

ggsurvplot(fit_uni_score_norm_15, conf.int = T,
           data=data,
           risk.table = TRUE, 
           surv.median.line = "hv", # Specify median survival
           legend.labs=c("low","high"), legend.title="Uni_score_norm",  
           palette=c("darkblue","cyan3"), 
           title="Kaplan-Meier Curves by country for outgoing people")

log_rank_test <- survdiff(Surv(Time, Cens==FALSE) ~ uni_score_norm_15, data=data)
log_rank_test
# 1 no evidence 

#############################################################################
#COX MODELS

###country

cox_country<-coxph(Surv(Time, Cens==FALSE) ~ country, data = data)
cox_country
summary(cox_country)
# significant 

plot(predict(cox_country), residuals(cox_country, type='martingale'),
     xlab='Fitted values', ylab='Martingale residuals', main='Residual Plot', las=1)
# Add a line for residual=0
abline(h=0, col='red')
# Fit a smoother for the points
lines(smooth.spline(predict(cox_country), residuals(cox_country, type='martingale')), col='blue')
ggcoxdiagnostics(cox_country, type = "martingale")

ggcoxdiagnostics(cox_country, type = "deviance")

ggcoxdiagnostics(cox_country, type = "schoenfeld")

test.ph <- cox.zph(cox_country) # H0: Hazards are proportional we dont want to reject 
test.ph
# 0.02


### uni_score

cox_uni_score<-coxph(Surv(Time, Cens==FALSE) ~ uni_score, data = data)
cox_uni_score
summary(cox_uni_score)
# not significant

plot(predict(cox_uni_score), residuals(cox_uni_score, type='martingale'),
     xlab='Fitted values', ylab='Martingale residuals', main='Residual Plot', las=1)
# Add a line for residual=0
abline(h=0, col='red')
# Fit a smoother for the points
lines(smooth.spline(predict(cox_uni_score), residuals(cox_uni_score, type='martingale')), col='blue')
ggcoxdiagnostics(cox_uni_score, type = "martingale")

ggcoxdiagnostics(cox_uni_score, type = "deviance")

ggcoxdiagnostics(cox_uni_score, type = "schoenfeld")

test.ph <- cox.zph(cox_uni_score) # H0: Hazards are proportional we dont want to reject 
test.ph



### LPPI

cox_lppi<-coxph(Surv(Time, Cens==FALSE) ~ LPPI, data = data)
cox_lppi
summary(cox_lppi)
# significant -> welfare makes people leave just fuckin kill it

plot(predict(cox_lppi), residuals(cox_lppi, type='martingale'),
     xlab='Fitted values', ylab='Martingale residuals', main='Residual Plot', las=1)
# Add a line for residual=0
abline(h=0, col='red')
# Fit a smoother for the points
lines(smooth.spline(predict(cox_lppi), residuals(cox_lppi, type='martingale')), col='blue')
ggcoxdiagnostics(cox_lppi, type = "martingale")

ggcoxdiagnostics(cox_lppi, type = "deviance")

ggcoxdiagnostics(cox_lppi, type = "schoenfeld")

test.ph <- cox.zph(cox_lppi) # H0: Hazards are proportional we dont want to reject 
test.ph
#proportional


