rm(list=ls())
library(car)
library(ISLR2)
library(np)
library(splines)
library(fda)
library(pbapply)
library(parallel)

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



#######################################################################################
# NATURAL CUBIC B-SPLINES
# SET KNOTS ACCORDING TO EVENLY SPACED QUANTILES
# find optimal number of knots with AIC


# set the grid of dof search
knots.grid <- seq(3, 20, by=1)

# function to find AIC for a certain bandwidth
find.AIC <- function(curr.knots){
  int.knots <- quantile(uni_score, seq(0, 1, length.out=curr.knots+2))[3:(curr.knots)]
  ext.knots <- quantile(uni_score, seq(0, 1, length.out=curr.knots+2))[c(2, curr.knots+1)]
  mod.curr <- lm(log_flow ~ ns(uni_score, knots = int.knots, Boundary.knots = ext.knots))
  return(extractAIC(mod.curr)[2]) 
}


# define cores for parallel computation
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = list("log_flow", "uni_score","find.AIC", "ns"))
AIC_wrapper <- function(curr.knots){find.AIC(curr.knots)} 

# find AIC for all grid points with parallel computing
all_AIC <- pbsapply(knots.grid, AIC_wrapper, cl=cl)


# optimal number of neighbors, span and corresponding RMSE
opt.knots.AIC <- knots.grid[which.min(all_AIC)]
opt.knots.AIC
min(all_AIC)

x11()
plot(knots.grid, all_AIC, type ='l', lwd='2', main='AIC vs #knots')
abline(v=opt.knots.AIC, col='red', lwd=2)


int.knots <- quantile(uni_score, seq(0, 1, length.out=opt.knots.AIC+2))[3:(opt.knots.AIC)]
ext.knots <- quantile(uni_score, seq(0, 1, length.out=opt.knots.AIC+2))[c(2, opt.knots.AIC+1)]

mod.nat.3.opt.AIC <- lm(log_flow ~ ns(uni_score, knots = int.knots, Boundary.knots = ext.knots))
summary(mod.nat.3.opt.AIC)

x11()
par(mfrow=c(2,2))
plot(mod.nat.3.opt.AIC)
shapiro.test(mod.nat.3.opt.AIC$residuals)

# I plot fitted values +- 2 * standard errors)
xx <- seq(min(uni_score), max(uni_score), length.out = 1000)
preds <- predict(mod.nat.3.opt.AIC, data.frame(uni_score = xx), se.fit=T)
x11()
plot(uni_score, log_flow, main = 'Natural Cubic splines, optimal AIC')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)
abline(v=c(int.knots, ext.knots), lty=2, col='gray')


RMSE.AIC <- sqrt(mean((mod.nat.3.opt.AIC$residuals)^2))
RMSE.AIC



############################################################################################
# NATURAL CUBIC B-SPLINES
# I find optimal dof with RMSE computed on STRATIFIED K-FOLD CROSS-VALIDATION

source('quantile_kfold.R')

# define the folds
n_fold <- 5
set.seed(1)
data.fr <- data.frame(log_flow=log_flow, uni_score = uni_score)
kfolds <- quantile_kfold(n_fold, data.fr, 'uni_score')

# function to find RMSE for a certain bandwidth
find.RMSE <- function(curr.knots, n_fold, kfolds){
  require(np)
  RMSE <- numeric(n_fold)
  for(i in 1:n_fold){
    curr_data <- kfolds$datasets[[i]]
    curr_miss <- kfolds$folds[[i]]
    int.knots <- quantile(curr_data$uni_score, seq(0, 1, length.out=curr.knots+2))[3:(curr.knots)]
    ext.knots <- quantile(curr_data$uni_score, seq(0, 1, length.out=curr.knots+2))[c(2, curr.knots+1)]
    mod.curr <- lm(log_flow ~ ns(uni_score, knots = int.knots, Boundary.knots = ext.knots), data = curr_data)
    RMSE[i] <- sqrt(mean((predict(mod.curr, newdata=curr_miss) - curr_miss$log_flow)^2))
  }
  return(mean(RMSE)) 
}


# set the grid of dof search
knots.grid <- seq(3, 20, by=1)

# define cores for parallel computation
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = list("n_fold", "find.RMSE", "kfolds", "ns"))
RMSE_wrapper <- function(curr.knots){find.RMSE(curr.knots, n_fold, kfolds)} 

# find RMSE for all grid points with parallel computing
all_RMSE <- pbsapply(knots.grid, RMSE_wrapper, cl=cl)

# optimal number of neighbors, span and corresponding RMSE
opt.knots.RMSE <- knots.grid[which.min(all_RMSE)]
opt.knots.RMSE
min(all_RMSE)

x11()
plot(knots.grid, all_RMSE, type ='l', lwd='2', main='Root MSE vs #knots')
abline(v=opt.knots.RMSE, col='red', lwd=2)



#######################################################################################
# NATURAL CUBIC B-SPLINES WITH OPTIMAL DOF


int.knots <- quantile(uni_score, seq(0, 1, length.out=opt.knots.RMSE+2))[3:(opt.knots)]
ext.knots <- quantile(uni_score, seq(0, 1, length.out=opt.knots.RMSE+2))[c(2, opt.knots+1)]

mod.nat.3.opt.RMSE <- lm(log_flow ~ ns(uni_score, knots = int.knots, Boundary.knots = ext.knots))
summary(mod.nat.3.opt.RMSE)

x11()
par(mfrow=c(2,2))
plot(mod.nat.3.opt.RMSE)
shapiro.test(mod.nat.3.opt.RMSE$residuals)

# I plot fitted values +- 2 * standard errors)
xx <- seq(min(uni_score), max(uni_score), length.out = 1000)
preds <- predict(mod.nat.3.opt.RMSE, data.frame(uni_score = xx), se.fit=T)
x11()
plot(uni_score, log_flow, main = 'Natural Cubic splines, optimal RMSE')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)
abline(v=c(int.knots, ext.knots), lty=2, col='gray')


RMSE <- sqrt(mean((mod.nat.3.opt.RMSE$residuals)^2))
RMSE

extractAIC(mod.nat.3.opt.RMSE)