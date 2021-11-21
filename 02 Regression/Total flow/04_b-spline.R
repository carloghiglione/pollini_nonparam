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
# CUBIC B-SPLINES (ORDER = 3), SET THE DEGREES OF FREEDOM
# KNOTS ARE PLACED AT THE CENTRAL (DOF - ORDER) EVEN QUANTILES (EXCLUDING 0% AND 100%)

# order of splines
order.s <- 3

# degrees of freedom (knots = dof - orders)
dof <- 7

mod.spline.3 <- lm(log_flow ~ bs(uni_score, df=dof, degree=order.s))
summary(mod.spline.3)

# extract the knots
used.knots <- attributes(bs(uni_score, df=dof, degree=order.s))$knots

x11()
par(mfrow=c(2,2))
plot(mod.spline.3)
shapiro.test(mod.spline.3$residuals)

# I plot fitted values +- 2 * standard errors)
xx <- seq(min(uni_score), max(uni_score), length.out = 1000)
preds <- predict(mod.spline.3, data.frame(uni_score = xx), se.fit=T)
x11()
plot(uni_score, log_flow, main = 'Cubic splines')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)
abline(v=used.knots, lty=2, col='gray')



############################################################################################
# I find RMSE and optimal span with CROSS-VALIDATION

source('quantile_kfold.R')

# define the folds
n_fold <- 5
set.seed(1)
data.fr <- data.frame(log_flow=log_flow, uni_score = uni_score)
kfolds <- quantile_kfold(n_fold, data.fr, 'uni_score')

# function to find RMSE for a certain bandwidth
find.RMSE <- function(curr.dof, n_fold, kfolds){
  require(np)
  RMSE <- numeric(n_fold)
  for(i in 1:n_fold){
    curr_data <- kfolds$datasets[[i]]
    curr_miss <- kfolds$folds[[i]]
    mod.curr <- lm(log_flow ~ bs(uni_score, df=curr.dof, degree=3), data = curr_data)
    RMSE[i] <- sqrt(mean((predict(mod.curr, newdata=curr_miss) - curr_miss$log_flow)^2))
  }
  return(mean(RMSE)) 
}


# set the grid of dof search
dof.grid <- seq(4, 20, by=1)

# define cores for parallel computation
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = list("n_fold", "find.RMSE", "kfolds", "bs"))
RMSE_wrapper <- function(curr.dof){find.RMSE(curr.dof, n_fold, kfolds)} 

# find RMSE for all grid points with parallel computing
all_RMSE <- pbsapply(dof.grid, RMSE_wrapper, cl=cl)

# optimal number of neighbors, span and corresponding RMSE
opt.dof <- dof.grid[which.min(all_RMSE)]
opt.dof
min(all_RMSE)

x11()
plot(dof.grid, all_RMSE, type ='l', lwd='2', main='Root MSE vs span')
abline(v=opt.dof, col='red', lwd=2)



#######################################################################################
# CUBIC B-SPLINES (ORDER = 3) WITH OPTIMAL DEGREE OF FREEDOM

mod.spline.3.opt <- lm(log_flow ~ bs(uni_score, df=opt.dof, degree=order.s))
summary(mod.spline.3.opt)

# extract the knots
used.knots <- attributes(bs(uni_score, df=opt.dof, degree=order.s))$knots

x11()
par(mfrow=c(2,2))
plot(mod.spline.3.opt)
shapiro.test(mod.spline.3.opt$residuals)

# I plot fitted values +- 2 * standard errors)
xx <- seq(min(uni_score), max(uni_score), length.out = 1000)
preds <- predict(mod.spline.3.opt, data.frame(uni_score = xx), se.fit=T)
x11()
plot(uni_score, log_flow, main = 'Cubic splines')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)
abline(v=used.knots, lty=2, col='gray')

RMSE <- sqrt(mean((mod.spline.3.opt$residuals)^2))
RMSE
