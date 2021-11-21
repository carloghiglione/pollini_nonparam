rm(list=ls())
library(car)
library(ISLR2)
library(np)
library(splines)
library(fda)
library(pbapply)
library(parallel)

#################
# load data

full_tab <- read.csv('full_tab_31_10_21.csv')

flow <- full_tab$tot_in - full_tab$tot_out    # in - out
flow.norm <- flow / full_tab$num_ric

countries <- full_tab$Country

Z <- data.frame(flow.norm, full_tab[, c(4,6,8,17)])
x11()
pairs(Z)



##################################################################
# do PCA on the regressors
tab <- Z[,-1]

# standardized data
tab.std <- data.frame(scale(tab))

# PCA
pca.std <- princomp(tab.std, scores = T)
print(summary(pca.std))

# collect scores of first PCA
scores.pc1  <- pca.std$scores[,1]

# remove NA (identified previous script)
flow.norm[4] <- NA    # tolgo Cipro
flow.norm[19] <- NA    # tolgo Lussemburgo

flow.norm <- flow.norm[-c(4,19)]
scores.pc1 <- scores.pc1[-c(4,19)]



###############################################################
# SMOOTHING CUBIC B-SPLINES 
# SELECT OPTIMAL LAMBDA WITH CROSS-VALIDATION

fit.s.spline.loocv <- smooth.spline(scores.pc1, flow.norm, cv=T)
fit.s.spline.gcv <- smooth.spline(scores.pc1, flow.norm, cv=F)

x11()
plot(scores.pc1, flow.norm, main = 'Smoothing cubic splines')
lines(fit.s.spline.loocv, col='red', lwd=2)
#lines(fit.s.spline.gcv, col='blue', lwd=2)
legend('bottomright', legend = c('LOOCV', 'GCV'), fill = c('red', 'blue'))



############################################################################################
# I find RMSE and optimal span with CROSS-VALIDATION

source('quantile_kfold.R')

# define the folds
n_fold <- 5
set.seed(1)
data.fr <- data.frame(flow.norm=flow.norm, scores.pc1 = scores.pc1)
kfolds <- quantile_kfold(n_fold, data.fr, 'scores.pc1')

# function to find RMSE for a certain bandwidth
find.RMSE <- function(curr.lam, n_fold, kfolds){
  require(np)
  RMSE <- numeric(n_fold)
  for(i in 1:n_fold){
    curr_data <- kfolds$datasets[[i]]
    curr_miss <- kfolds$folds[[i]]
    mod.curr <- smooth.spline(curr_data$scores.pc1, curr_data$flow.norm, lambda = curr.lam)
    RMSE[i] <- sqrt(mean((predict(mod.curr, curr_miss$scores.pc1)$y - curr_miss$flow.norm)^2))
  }
  return(mean(RMSE)) 
}


# set the grid of dof search
lam.grid <- seq(0.01, 1, by=0.001)

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

fit.smooth.opt <- smooth.spline(scores.pc1, flow.norm, lambda = opt.lam)

x11()
plot(scores.pc1, flow.norm, main = 'Smoothing cubic splines')
lines(fit.smooth.opt, col='red', lwd=2)
points(fit.smooth.opt$x, fit.smooth.opt$y, pch='x')

RMSE <- sqrt(mean((predict(fit.smooth.opt, scores.pc1)$y- flow.norm)^2))
RMSE
