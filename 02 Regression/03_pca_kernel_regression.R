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



##################################################################################
# KERNEL LOCAL AVERAGING REGRESSION
# OPTIMAL BANDWIDTH DEFINED WITH NATIVE CV.LS PARAMETER

mod.ker <- npreg(flow.norm ~ scores.pc1, ckertype = 'gaussian', bwmethod="cv.ls")
summary(mod.ker)

xx <- seq(min(scores.pc1), max(scores.pc1), length.out = 1000)
preds <- predict(mod.ker, newdata = data.frame(scores.pc1 = xx), se.fit=T)
x11()
plot(scores.pc1, flow.norm, main = 'Gaussian kernel regression')
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
data.fr <- data.frame(flow.norm=flow.norm, scores.pc1 = scores.pc1)
kfolds <- quantile_kfold(n_fold, data.fr, 'scores.pc1')

# function to find RMSE for a certain bandwidth
find.RMSE.ker <- function(curr.bws, n_fold, kfolds){
  require(np)
  RMSE <- numeric(n_fold)
  for(i in 1:n_fold){
    curr_data <- kfolds$datasets[[i]]
    curr_miss <- kfolds$folds[[i]]
    mod.curr <- npreg(flow.norm ~ scores.pc1, ckertype = 'gaussian', bws=curr.bws, data = curr_data)
    RMSE[i] <- sqrt(mean((predict(mod.curr, newdata=curr_miss) - curr_miss$flow.norm)^2))
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

mod.ker.opt <- npreg(flow.norm ~ scores.pc1, ckertype = 'gaussian', bws=opt.bws)
summary(mod.ker.opt)

xx <- seq(min(scores.pc1), max(scores.pc1), length.out = 1000)
preds <- predict(mod.ker.opt, newdata = data.frame(scores.pc1 = xx), se.fit=T)
x11()
plot(scores.pc1, flow.norm, main = 'Gaussian kernel regression')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)

# RMSE
sqrt(mod.ker.opt$MSE)

graphics.off()