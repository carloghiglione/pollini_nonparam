rm(list=ls())
library(car)
library(ISLR2)
library(np)
library(splines)
library(fda)
library(pbapply)
library(parallel)
library(npreg)

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
# SMOOTHING CUBIC B-SPLINES (native R function)
# SELECT OPTIMAL LAMBDA WITH CROSS-VALIDATION

fit.s.spline.loocv <- smooth.spline(uni_score, log_flow, cv=T)
fit.s.spline.gcv <- smooth.spline(uni_score, log_flow, cv=F)

x11()
plot(uni_score, log_flow, main = 'Smoothing cubic splines')
lines(fit.s.spline.loocv, col='red', lwd=2)
#lines(fit.s.spline.gcv, col='blue', lwd=2)
legend('bottomright', legend = c('LOOCV', 'GCV'), fill = c('red', 'blue'))



###############################################################
# SMOOTHING CUBIC B-SPLINES
# find optimal number of knots with AIC

# set the grid of dof search
knots.grid <- seq(3, 15, by=1)

# function to find AIC for a certain bandwidth
find.AIC <- function(curr.knots){
  my_knots <- quantile(uni_score, seq(0, 1, length.out=curr.knots+2))[2:curr.knots+1]
  mod.curr <- ss(uni_score, log_flow, method = 'AIC', knots = my_knots)
  return(mod.curr$aic) 
}


# define cores for parallel computation
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = list("log_flow", "uni_score","find.AIC", "ss"))
AIC_wrapper <- function(curr.knots){find.AIC(curr.knots)} 

# find AIC for all grid points with parallel computing
all_AIC <- pbsapply(knots.grid, AIC_wrapper, cl=cl)


# optimal number of neighbors, span and corresponding RMSE
opt.knots.AIC <- knots.grid[which.min(all_AIC)]
opt.knots.AIC
min(all_AIC)

x11()
plot(knots.grid, all_AIC, type ='l', lwd='2', main='AIC vs n° knots', ylab = 'AIC', xlab='n° knots')
abline(v=opt.knots.AIC, col='red', lwd=2)


my_knots <- quantile(uni_score, seq(0, 1, length.out=opt.knots.AIC+2))[2:opt.knots.AIC+1]
s.spline.opt.aic <- ss(uni_score, log_flow, method = 'AIC', knots = my_knots)
summary(s.spline.opt.aic)
s.spline.opt.aic$aic



x11()
plot(s.spline.opt.aic, xlab='Uni Score', ylab='Log(Flow)', ylim=c(4,12.5))              # plots 95% confidence interval
points(uni_score, log_flow)
#abline(v=s.spline.opt.aic$fit$knot, lty=2, col='gray')

RMSE.AIC <- sqrt(mean((predict(s.spline.opt.aic, uni_score)$y - log_flow)^2))
RMSE.AIC



############################################################################################
# I find optimal dof with RMSE computed on STRATIFIED K-FOLD CROSS-VALIDATION

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
    #mod.curr <- smooth.spline(curr_data$uni_score, curr_data$log_flow, lambda = curr.lam)            # vecchio
    #mod.curr <- ss(curr_data$uni_score, curr_data$log_flow, lambda = curr.lam, knots = used_knots)   # opzione meno smoothing, migliore RMSE
    mod.curr <- ss(curr_data$uni_score, curr_data$log_flow, lambda = curr.lam)                        # opzione più smoothing, peggiore RMSE
    RMSE[i] <- sqrt(mean((predict(mod.curr, curr_miss$uni_score)$y - curr_miss$log_flow)^2))
  }
  return(mean(RMSE)) 
}


# set the grid of dof search
#lam.grid <- seq(1e-6, 1e-3, by=1e-6)      # opzione meno smoothing, migliore RMSE
lam.grid <- 10^(-seq(2, 5, length=1000))   # opzione più smoothing, peggiore RMSE

# set the used knots
used_knots <- s.spline.opt.aic$fit$knot



# define cores for parallel computation
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = list("n_fold", "find.RMSE", "kfolds", "smooth.spline", 'ss', 'used_knots'))
RMSE_wrapper <- function(curr.lam){find.RMSE(curr.lam, n_fold, kfolds)} 

# find RMSE for all grid points with parallel computing
all_RMSE <- pbsapply(lam.grid, RMSE_wrapper, cl=cl)

# optimal number of neighbors, span and corresponding RMSE
opt.lam <- lam.grid[which.min(all_RMSE)]
opt.lam
min(all_RMSE)

x11()
plot(lam.grid, all_RMSE, type ='l', lwd='2', main='Root MSE vs Smoothing Param', log='x', ylab='RMSE', xlab='Lambda')
abline(v=opt.lam, col='red', lwd=2)



#######################################################################################
# SMOOTHING CUBIC B-SPLINES WITH OPTIMAL LAMBDA

#fit.smooth.opt <- ss(uni_score, log_flow, lambda = opt.lam, knots = used_knots)   # opzione meno smoothing, meglio RMSE
fit.smooth.opt <- ss(uni_score, log_flow, lambda = opt.lam)                        # opzione più smoothing, peggiore RMSE


x11()
plot(fit.smooth.opt, xlab='Uni Score', ylab='Log(Flow)', ylim=c(4,12.5))              # plots 95% confidence interval
points(uni_score, log_flow)
#abline(v=fit.smooth.opt$fit$knot, lty=2, col='gray')


RMSE <- sqrt(mean((predict(fit.smooth.opt, uni_score)$y- log_flow)^2))
RMSE
