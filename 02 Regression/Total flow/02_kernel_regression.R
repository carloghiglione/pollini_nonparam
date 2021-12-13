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



##################################################################################
# KERNEL LOCAL AVERAGING REGRESSION
# OPTIMAL BANDWIDTH DEFINED WITH NATIVE CV.AIC PARAMETER

mod.ker <- npreg(log_flow ~ uni_score, ckertype = 'gaussian', bwmethod="cv.aic")
summary(mod.ker)

xx <- seq(min(uni_score), max(uni_score), length.out = 1000)
preds <- predict(mod.ker, newdata = data.frame(uni_score = xx), se.fit=T)
x11()
plot(uni_score, log_flow, main = 'Gaussian kernel regression, AIC model')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)

# bandwidth
mod.ker$bw

# RMSE
sqrt(mod.ker$MSE)



####################################################################################
# I find RMSE with CROSS-VALIDATION

source('quantile_kfold.R')

# define the folds
n_fold <- 5
set.seed(1)
data.fr <- data.frame(log_flow=log_flow, uni_score = uni_score)
kfolds <- quantile_kfold(n_fold, data.fr, 'uni_score')

# function to find RMSE for a certain bandwidth
find.RMSE.ker <- function(curr.bws, n_fold, kfolds){
  require(np)
  RMSE <- numeric(n_fold)
  for(i in 1:n_fold){
    curr_data <- kfolds$datasets[[i]]
    curr_miss <- kfolds$folds[[i]]
    mod.curr <- npreg(log_flow ~ uni_score, ckertype = 'gaussian', bws=curr.bws, data = curr_data)
    RMSE[i] <- sqrt(mean((predict(mod.curr, newdata=curr_miss) - curr_miss$log_flow)^2))
  }
  return(mean(RMSE)) 
}


# set the grid of bandwidth for search
bws.grid <- seq(0.1, 2.5, by=0.01)

# define cores for parallel computation
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = list("n_fold", "find.RMSE.ker", "kfolds"))
RMSE_wrapper <- function(curr.bws){find.RMSE.ker(curr.bws, n_fold, kfolds)}  

# find RMSE for all grid points with parallel computing
all_RMSE <- pbsapply(bws.grid, RMSE_wrapper, cl=cl)

# optimal bandwidth and corresponding RMSE
opt.bws <- bws.grid[which.min(all_RMSE)]
opt.bws
min(all_RMSE)

x11()
plot(bws.grid, all_RMSE, type ='l', lwd='2', main='Root MSE vs bandwidth')
abline(v=opt.bws, col='red', lwd=2)



#########################################################################################
# KERNEL LOCAL AVERAGING REGRESSION WITH OPTIMAL BANDWIDTH

mod.ker.opt <- npreg(log_flow ~ uni_score, ckertype = 'gaussian', bws=opt.bws)
summary(mod.ker.opt)

xx <- seq(min(uni_score), max(uni_score), length.out = 1000)
preds <- predict(mod.ker.opt, newdata = data.frame(uni_score = xx), se.fit=T)
x11()
plot(uni_score, log_flow, main = 'Gaussian kernel regression')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)

# RMSE
sqrt(mod.ker.opt$MSE)

graphics.off()
