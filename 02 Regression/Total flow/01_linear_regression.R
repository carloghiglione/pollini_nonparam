rm(list=ls())
library(car)

full_tab <- read.csv('full_tab_31_10_21.csv')

##################
# creo tabelle con variabili di interesse

log_flow <- log(full_tab$tot_in + full_tab$tot_out)    # in + out

uni_score <- full_tab$uni_score

countries <- full_tab$Country

##################
# visualize results
x11()
plot(uni_score, log_flow, main='Total flow vs Uni Score')



###################################################################################
# LINEAR REGRESSION ON FLOW
mod <- lm(log_flow ~ uni_score)
summary(mod)

x11()
par(mfrow=c(2,2))
plot(mod)
shapiro.test(mod$residuals)

x11()
b <- mod$coefficients
xx = seq(min(uni_score), max(uni_score), length = 100)
plot(uni_score, log_flow, main = "Total flow vs Uni Score")
lines(xx, b[1] + b[2]*xx, col='red')

RMSE <- sqrt(mean((mod$residuals)^2))
RMSE

extractAIC(mod)


############################################################################################
# K-FOLD CROSS-VALIDATION
# I use it to estimate the Root-Mean-Square-Error

source('quantile_kfold.R')

n_fold <- 5
set.seed(1)
data.fr <- data.frame(log_flow=log_flow, uni_score=uni_score)

kfolds <- quantile_kfold(n_fold, data.fr, 'uni_score')
RMSE <- numeric(n_fold)

for(i in 1:n_fold){
  curr_data <- kfolds$datasets[[i]]
  curr_miss <- kfolds$folds[[i]]
  mod.curr <- lm(log_flow ~ uni_score, data = curr_data)
  
  RMSE[i] <- sqrt(mean((predict(mod.curr, curr_miss) - curr_miss$log_flow)^2))
  
  x11()
  b <- mod.curr$coefficients
  xx = seq(min(uni_score), max(uni_score), length = 100)
  plot(curr_data$uni_score, curr_data$log_flow, main = "Normalized Flow vs pca", 
       xlim = range(uni_score), ylim = range(log_flow))
  lines(xx, b[1] + b[2]*xx, col='red')
  points(curr_miss$uni_score, curr_miss$log_flow, col='blue', pch=19)
  points(curr_miss$uni_score, predict(mod.curr, curr_miss), col='blue', pch='x')
}
RMSE <- mean(RMSE)
RMSE

graphics.off()
