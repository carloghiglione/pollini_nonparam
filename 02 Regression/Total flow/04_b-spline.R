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
# find optimal dof with AIC

# order of splines
order.s <- 3


# set the grid of dof search (knots = dof - orders)
dof.grid <- seq(4, 20, by=1)

# function to find AIC for a certain bandwidth
find.AIC <- function(curr.dof){
  mod.curr <- lm(log_flow ~ bs(uni_score, df=curr.dof, degree=order.s))
  return(extractAIC(mod.curr)[2]) 
}


# define cores for parallel computation
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = list("log_flow", "uni_score","find.AIC", "bs", "order.s"))
AIC_wrapper <- function(curr.dof){find.AIC(curr.dof)} 

# find AIC for all grid points with parallel computing
all_AIC <- pbsapply(dof.grid, AIC_wrapper, cl=cl)


# optimal number of neighbors, span and corresponding RMSE
opt.dof.AIC <- dof.grid[which.min(all_AIC)]
opt.dof.AIC
min(all_AIC)

x11()
plot(dof.grid, all_AIC, type ='l', lwd='2', main='AIC vs dof')
abline(v=opt.dof.AIC, col='red', lwd=2)


mod.spline.3.opt.AIC <- lm(log_flow ~ bs(uni_score, df=opt.dof.AIC, degree=order.s))
summary(mod.spline.3.opt.AIC)

# extract the knots
used.knots <- attributes(bs(uni_score, df=opt.dof.AIC, degree=order.s))$knots

x11()
par(mfrow=c(2,2))
plot(mod.spline.3.opt.AIC)
shapiro.test(mod.spline.3.opt.AIC$residuals)

# I plot fitted values +- 2 * standard errors)
xx <- seq(min(uni_score), max(uni_score), length.out = 1000)
preds <- predict(mod.spline.3.opt.AIC, data.frame(uni_score = xx), se.fit=T)
x11()
plot(uni_score, log_flow, main = 'Cubic splines, optimal AIC')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)
abline(v=used.knots, lty=2, col='gray')


RMSE.AIC <- sqrt(mean((mod.spline.3.opt.AIC$residuals)^2))
RMSE.AIC



############################################################################################
# CUBIC B-SPLINES (ORDER = 3), SET THE DEGREES OF FREEDOM
# I find optimal dof with RMSE computed on STRATIFIED K-FOLD CROSS-VALIDATION

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
plot(dof.grid, all_RMSE, type ='l', lwd='2', main='Root MSE vs dof')
abline(v=opt.dof, col='red', lwd=2)



#######################################################################################
# CUBIC B-SPLINES (ORDER = 3) WITH OPTIMAL DEGREE OF FREEDOM wrt RMSE

mod.spline.3.opt.RMSE <- lm(log_flow ~ bs(uni_score, df=opt.dof, degree=order.s))
summary(mod.spline.3.opt.RMSE)

# extract the knots
used.knots <- attributes(bs(uni_score, df=opt.dof, degree=order.s))$knots

x11()
par(mfrow=c(2,2))
plot(mod.spline.3.opt.RMSE)
shapiro.test(mod.spline.3.opt.RMSE$residuals)

# I plot fitted values +- 2 * standard errors)
xx <- seq(min(uni_score), max(uni_score), length.out = 1000)
preds <- predict(mod.spline.3.opt.RMSE, data.frame(uni_score = xx), se.fit=T)
x11()
plot(uni_score, log_flow, main = 'Cubic splines, optimal RMSE')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)
abline(v=used.knots, lty=2, col='gray')

RMSE <- sqrt(mean((mod.spline.3.opt.RMSE$residuals)^2))
RMSE

extractAIC(mod.spline.3.opt.RMSE)