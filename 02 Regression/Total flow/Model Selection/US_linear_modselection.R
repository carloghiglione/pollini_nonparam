#we want to build a model that gives a better interpretation of the University score in terms of more quantifiable features
#====
library(dplyr)

data <- full_join(x = full_tab_31_10_21, y = Trips_labeled, by = c("Country" = "Country"))
pat <- data$Patent.applications..nonresidents + data$Patent.applications..residents
#====

#We start by looking for a simple linear regression in terms of education related parameters
#it should be noted that during the Applied Statistics project, the best linear model for uni_score was based on patents, GDP and GDP^2
#Now we are also considering now topic-related variables
us <- lm(data$uni_score ~ data$Citations.per.document + data$num_staff + data$stud_per_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
summary(us)
shapiro.test(us$residuals)


#To assess the validity of this model, since residuals are not gaussian, we need a Permutation test

#====
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ 1)
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us)$fstatistic[1]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$stud_per_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$fstatistic[1]
  pb$tick()
}
sum(T_stat>=T0)/B
#====

#p_val per B = 5000 = 0.0374

#It is a valid model, so we try to reduce it by backward selection

#====
p_oneatthetime <- numeric(length(us$coefficients)-1)

#1
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$num_staff + data$stud_per_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us)$coefficients[1+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$stud_per_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[1+1,3]
  pb$tick()
}
p_oneatthetime[1] <- sum(abs(T_stat)>=abs(T0))/B

#2
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$stud_per_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us)$coefficients[2+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$stud_per_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[2+1,3]
  pb$tick()
}
p_oneatthetime[2] <- sum(abs(T_stat)>=abs(T0))/B

#3
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$num_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us)$coefficients[3+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$stud_per_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[3+1,3]
  pb$tick()
}
p_oneatthetime[3] <- sum(abs(T_stat)>=abs(T0))/B

#4
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$num_staff + data$stud_per_staff + pat + data$GDP + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us)$coefficients[4+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$stud_per_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[4+1,3]
  pb$tick()
}
p_oneatthetime[4] <- sum(abs(T_stat)>=abs(T0))/B

#5
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$num_staff + data$stud_per_staff + data$Reasearch + data$GDP + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us)$coefficients[5+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$stud_per_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[5+1,3]
  pb$tick()
}
p_oneatthetime[5] <- sum(abs(T_stat)>=abs(T0))/B

#6
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$num_staff + data$stud_per_staff + data$Reasearch + pat + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us)$coefficients[6+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$stud_per_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[6+1,3]
  pb$tick()
}
p_oneatthetime[6] <- sum(abs(T_stat)>=abs(T0))/B

#7
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$num_staff + data$stud_per_staff + data$Reasearch + pat + data$GDP)
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us)$coefficients[7+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$stud_per_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[7+1,3]
  pb$tick()
}
p_oneatthetime[7] <- sum(abs(T_stat)>=abs(T0))/B

#====
p_oneatthetime
# 0.2858 0.8142 0.9418 0.2344 0.1084 0.5716 0.5366

#we eliminate stud_per_staff
#Reduced model
us2 <- lm(data$uni_score ~ data$Citations.per.document + data$num_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
summary(us2)

#we iterate the procedure
#====
p_oneatthetime <- numeric(length(us2$coefficients)-1)

#1
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$num_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us2)$coefficients[1+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[1+1,3]
  pb$tick()
}
p_oneatthetime[1] <- sum(abs(T_stat)>=abs(T0))/B

#2
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$Reasearch + pat + data$GDP + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us2)$coefficients[2+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[2+1,3]
  pb$tick()
}
p_oneatthetime[2] <- sum(abs(T_stat)>=abs(T0))/B

#3
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$num_staff + pat + data$GDP + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us2)$coefficients[3+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[3+1,3]
  pb$tick()
}
p_oneatthetime[3] <- sum(abs(T_stat)>=abs(T0))/B

#4
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$num_staff + data$Reasearch + data$GDP + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us2)$coefficients[4+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[4+1,3]
  pb$tick()
}
p_oneatthetime[4] <- sum(abs(T_stat)>=abs(T0))/B

#5
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$num_staff + data$Reasearch + pat + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us2)$coefficients[5+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[5+1,3]
  pb$tick()
}
p_oneatthetime[5] <- sum(abs(T_stat)>=abs(T0))/B

#6
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$num_staff + data$Reasearch + pat + data$GDP)
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us2)$coefficients[6+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$num_staff + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[6+1,3]
  pb$tick()
}
p_oneatthetime[6] <- sum(abs(T_stat)>=abs(T0))/B

#====
p_oneatthetime
#0.2628 0.8130 0.1992 0.1054 0.5084 0.4628

#reduced model: eliminate num_staff
us3 <- lm(data$uni_score ~ data$Citations.per.document + data$Reasearch + pat + data$GDP + I(data$GDP^2))
summary(us3)

#====
p_oneatthetime <- numeric(length(us3$coefficients)-1)

#1
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Reasearch + pat + data$GDP + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us3)$coefficients[1+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[1+1,3]
  pb$tick()
}
p_oneatthetime[1] <- sum(abs(T_stat)>=abs(T0))/B

#2
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + pat + data$GDP + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us3)$coefficients[2+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[2+1,3]
  pb$tick()
}
p_oneatthetime[2] <- sum(abs(T_stat)>=abs(T0))/B

#3
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$Reasearch + data$GDP + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us3)$coefficients[3+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[3+1,3]
  pb$tick()
}
p_oneatthetime[3] <- sum(abs(T_stat)>=abs(T0))/B

#4
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$Reasearch + pat + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us3)$coefficients[4+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[4+1,3]
  pb$tick()
}
p_oneatthetime[4] <- sum(abs(T_stat)>=abs(T0))/B

#5
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$Reasearch + pat + data$GDP)
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us3)$coefficients[5+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$Reasearch + pat + data$GDP + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[5+1,3]
  pb$tick()
}
p_oneatthetime[5] <- sum(abs(T_stat)>=abs(T0))/B

#====
p_oneatthetime
# 0.2596 0.1954 0.0108 0.5080 0.4628

#reduced model: remove GDP
us4 <- lm(data$uni_score ~ data$Citations.per.document + data$Reasearch + pat + I(data$GDP^2))
summary(us4)

#====
p_oneatthetime <- numeric(length(us4$coefficients)-1)

#1
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Reasearch + pat + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us4)$coefficients[1+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$Reasearch + pat + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[1+1,3]
  pb$tick()
}
p_oneatthetime[1] <- sum(abs(T_stat)>=abs(T0))/B

#2
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + pat + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us4)$coefficients[2+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$Reasearch + pat + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[2+1,3]
  pb$tick()
}
p_oneatthetime[2] <- sum(abs(T_stat)>=abs(T0))/B

#3
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$Reasearch + I(data$GDP^2))
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us4)$coefficients[3+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$Reasearch + pat + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[3+1,3]
  pb$tick()
}
p_oneatthetime[3] <- sum(abs(T_stat)>=abs(T0))/B

#4
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$Reasearch + pat)
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us4)$coefficients[4+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$Reasearch + pat + I(data$GDP^2))
  T_stat[i] <- summary(mod_poop)$coefficients[4+1,3]
  pb$tick()
}
p_oneatthetime[4] <- sum(abs(T_stat)>=abs(T0))/B


#====

p_oneatthetime
#0.0098 0.2066 0.0082 0.3990

#reduced model: remove gdp^2
us5 <- lm(data$uni_score ~ data$Citations.per.document + data$Reasearch + pat)
summary(us5)

#====
p_oneatthetime <- numeric(length(us5$coefficients)-1)

#1
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Reasearch + pat)
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us5)$coefficients[1+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$Reasearch + pat)
  T_stat[i] <- summary(mod_poop)$coefficients[1+1,3]
  pb$tick()
}
p_oneatthetime[1] <- sum(abs(T_stat)>=abs(T0))/B

#2
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + pat)
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us5)$coefficients[2+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$Reasearch + pat)
  T_stat[i] <- summary(mod_poop)$coefficients[2+1,3]
  pb$tick()
}
p_oneatthetime[2] <- sum(abs(T_stat)>=abs(T0))/B

#3
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document + data$Reasearch)
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us5)$coefficients[3+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + data$Reasearch + pat)
  T_stat[i] <- summary(mod_poop)$coefficients[3+1,3]
  pb$tick()
}
p_oneatthetime[3] <- sum(abs(T_stat)>=abs(T0))/B

#====

p_oneatthetime
#0.0064 0.1894 0.0068

#reduced model: eliminate research
us6 <- lm(data$uni_score ~ data$Citations.per.document + pat)
summary(us6)
#WE FALL BELOW 50% EXPLAINED VARIABILITY

#====
p_oneatthetime <- numeric(length(us6$coefficients)-1)

#1
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ pat)
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us6)$coefficients[1+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + pat)
  T_stat[i] <- summary(mod_poop)$coefficients[1+1,3]
  pb$tick()
}
p_oneatthetime[1] <- sum(abs(T_stat)>=abs(T0))/B

#2
B <- 5000
set.seed(17121998)
mod0 <- lm(data$uni_score ~ data$Citations.per.document)
res <- mod0$residuals
fitted <- mod0$fitted.values
T0 <- summary(us6)$coefficients[2+1,3]
T_stat <- numeric(B)
library(progress)
pb <- progress_bar$new(total = B)
pb$tick(0)
for(i in 1:B){
  perm <- sample(1:length(res))
  y_pool <- fitted + res[perm]
  mod_poop <- lm(y_pool ~  data$Citations.per.document + pat)
  T_stat[i] <- summary(mod_poop)$coefficients[2+1,3]
  pb$tick()
}
p_oneatthetime[2] <- sum(abs(T_stat)>=abs(T0))/B

#====

p_oneatthetime
#0.0048 0.0138

#final model: linear regression based on citations per document and number of patents

summary(us6)
summary(us5)
summary(us)

#we can however try to do better
logus <- log(data$uni_score + 1)

x11()
plot(data$Citations.per.document, data$uni_score)
plot(log(pat), logus)

lus <- lm(logus ~ log(pat))

summary(lus)
shapiro.test(lus$residuals)
qqnorm(lus$residuals)
qqline(lus$residuals)
#Gaussian!

lus_pot <- lm(logus ~ log(pat) + log(data$Citations.per.document))
summary(lus_pot)
shapiro.test(lus_pot$residuals)
qqnorm(lus_pot$residuals)
qqline(lus_pot$residuals)
#Gaussian!

#So we can assess the positive correlation between the quality of universities in a country and the entrepreneurial research activity
#with an extra focus of not only the size of said universitarial world but also its international recognition (interpreted by number of citations per document)