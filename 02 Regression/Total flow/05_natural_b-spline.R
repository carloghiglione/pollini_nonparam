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

# number of knots
n_knots <- 5

int.knots <- quantile(uni_score, seq(0, 1, length.out=n_knots+2))[3:(n_knots)]
ext.knots <- quantile(uni_score, seq(0, 1, length.out=n_knots+2))[c(2, n_knots+1)]

mod.nat.spline.3 <- lm(log_flow ~ ns(uni_score, knots = int.knots, Boundary.knots = ext.knots))
summary(mod.nat.spline.3)

x11()
par(mfrow=c(2,2))
plot(mod.nat.spline.3)
shapiro.test(mod.nat.spline.3$residuals)

# I plot fitted values +- 2 * standard errors)
xx <- seq(min(uni_score), max(uni_score), length.out = 1000)
preds <- predict(mod.nat.spline.3, data.frame(uni_score = xx), se.fit=T)
x11()
plot(uni_score, log_flow, main = 'Natural Cubic splines')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)
abline(v=c(int.knots, ext.knots), lty=2, col='gray')



############################################################################################
# I find RMSE and optimal span with CROSS-VALIDATION

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
opt.knots <- knots.grid[which.min(all_RMSE)]
opt.knots
min(all_RMSE)

x11()
plot(knots.grid, all_RMSE, type ='l', lwd='2', main='Root MSE vs span')
abline(v=opt.knots, col='red', lwd=2)



#######################################################################################
# NATURAL CUBIC B-SPLINES WITH OPTIMAL DOF


int.knots <- quantile(uni_score, seq(0, 1, length.out=opt.knots+2))[3:(opt.knots)]
ext.knots <- quantile(uni_score, seq(0, 1, length.out=opt.knots+2))[c(2, opt.knots+1)]

mod.nat.3.opt <- lm(log_flow ~ ns(uni_score, knots = int.knots, Boundary.knots = ext.knots))
summary(mod.nat.3.opt)

x11()
par(mfrow=c(2,2))
plot(mod.nat.3.opt)
shapiro.test(mod.nat.3.opt$residuals)

# I plot fitted values +- 2 * standard errors)
xx <- seq(min(uni_score), max(uni_score), length.out = 1000)
preds <- predict(mod.nat.3.opt, data.frame(uni_score = xx), se.fit=T)
x11()
plot(uni_score, log_flow, main = 'Natural Cubic splines')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)
abline(v=c(int.knots, ext.knots), lty=2, col='gray')


RMSE <- sqrt(mean((mod.nat.3.opt$residuals)^2))
RMSE