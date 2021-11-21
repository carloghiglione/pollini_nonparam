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



###############################################################
# SMOOTHING CUBIC B-SPLINES 
# SELECT OPTIMAL LAMBDA WITH CROSS-VALIDATION

fit.s.spline.loocv <- smooth.spline(uni_score, log_flow, cv=T)
fit.s.spline.gcv <- smooth.spline(uni_score, log_flow, cv=F)

x11()
plot(uni_score, log_flow, main = 'Smoothing cubic splines')
lines(fit.s.spline.loocv, col='red', lwd=2)
#lines(fit.s.spline.gcv, col='blue', lwd=2)
legend('bottomright', legend = c('LOOCV', 'GCV'), fill = c('red', 'blue'))


############################################################################################
# I find RMSE and optimal span with CROSS-VALIDATION

source('quantile_kfold.R')

# define the folds
n_fold <- 5
set.seed(1)
data.fr <- data.frame(log_flow=log_flow, uni_score = uni_score)
kfolds <- quantile_kfold(n_fold, data.fr, 'uni_score')

# function to find RMSE for a certain bandwidth
find.RMSE <- function(curr.lam, n_fold, kfolds){
  require(np)
  RMSE <- numeric(n_fold)
  for(i in 1:n_fold){
    curr_data <- kfolds$datasets[[i]]
    curr_miss <- kfolds$folds[[i]]
    mod.curr <- smooth.spline(curr_data$uni_score, curr_data$log_flow, lambda = curr.lam)
    RMSE[i] <- sqrt(mean((predict(mod.curr, curr_miss$uni_score)$y - curr_miss$log_flow)^2))
  }
  return(mean(RMSE)) 
}


# set the grid of dof search
lam.grid <- seq(1e-7, 1e-3, by=1e-7)

# define cores for parallel computation
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = list("n_fold", "find.RMSE", "kfolds", "smooth.spline"))
RMSE_wrapper <- function(curr.lam){find.RMSE(curr.lam, n_fold, kfolds)} 

# find RMSE for all grid points with parallel computing
all_RMSE <- pbsapply(lam.grid, RMSE_wrapper, cl=cl)

# optimal number of neighbors, span and corresponding RMSE
opt.lam <- lam.grid[which.min(all_RMSE)]
opt.lam
min(all_RMSE)

x11()
plot(lam.grid, all_RMSE, type ='l', lwd='2', main='Root MSE vs span')
abline(v=opt.lam, col='red', lwd=2)



#######################################################################################
# SMOOTHING CUBIC B-SPLINES WITH OPTIMAAL LAMBDA

fit.smooth.opt <- smooth.spline(uni_score, log_flow, lambda = opt.lam)

x11()
plot(uni_score, log_flow, main = 'Smoothing cubic splines')
lines(fit.smooth.opt, col='red', lwd=2)
points(fit.smooth.opt$x, fit.smooth.opt$y, pch='x', col='blue')

RMSE <- sqrt(mean((predict(fit.smooth.opt, uni_score)$y- log_flow)^2))
RMSE
