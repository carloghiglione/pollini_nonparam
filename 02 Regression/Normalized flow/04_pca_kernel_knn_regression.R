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
# MIXED GAUSSIAN KERNEL & KNN LOCAL LINEAR REGRESSION
# NUMBER OF NEIGHBORS SELECTION WITH CROSS-VALIDATION ON AIC

# direcly select number of neighbors
k <- 5
# select span
n <- length(flow.norm)
span <- 0.05
k <- as.integer(span*n)

mod.knn <- npreg(flow.norm ~ scores.pc1, bwtype='adaptive_nn', bwmethod="cv.aic", ckertype = 'gaussian')
summary(mod.knn)

xx <- seq(min(scores.pc1), max(scores.pc1), length.out = 1000)
preds <- predict(mod.knn, newdata = data.frame(scores.pc1 = xx), se.fit=T)
x11()
plot(scores.pc1, flow.norm, main = 'Gaussian Kernel - KNN regression, AIC')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)

# number of neighbors
mod.knn$bw

# RMSE
sqrt(mod.knn$MSE)



############################################################################################
# I find RMSE and optimal span with CROSS-VALIDATION

source('quantile_kfold.R')

# define the folds
n_fold <- 5
set.seed(1)
data.fr <- data.frame(flow.norm=flow.norm, scores.pc1 = scores.pc1)
kfolds <- quantile_kfold(n_fold, data.fr, 'scores.pc1')

# function to find RMSE for a certain bandwidth
find.RMSE.knn <- function(curr.bws, n_fold, kfolds){
  require(np)
  RMSE <- numeric(n_fold)
  for(i in 1:n_fold){
    curr_data <- kfolds$datasets[[i]]
    curr_miss <- kfolds$folds[[i]]
    mod.curr <- npreg(flow.norm ~ scores.pc1, bwtype='adaptive_nn', ckertype = 'gaussian', bws=curr.bws, data = curr_data)
    RMSE[i] <- sqrt(mean((predict(mod.curr, newdata=curr_miss) - curr_miss$flow.norm)^2))
  }
  return(mean(RMSE)) 
}


# set the grid of nearest neighbors search
neigh.grid <- seq(1, ceiling(n/2), by=1)
span.grid <- neigh.grid/n

# define cores for parallel computation
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = list("n_fold", "find.RMSE.knn", "kfolds"))
RMSE_wrapper <- function(curr.bws){find.RMSE.knn(curr.bws, n_fold, kfolds)} 

# find RMSE for all grid points with parallel computing
all_RMSE <- pbsapply(neigh.grid, RMSE_wrapper, cl=cl)

# optimal number of neighbors, span and corresponding RMSE
opt.neigh <- neigh.grid[which.min(all_RMSE)]
opt.neigh
opt.span <- opt.neigh/n
opt.span
min(all_RMSE)

x11()
plot(span.grid, all_RMSE, type ='l', lwd='2', main='Root MSE vs span')
abline(v=opt.span, col='red', lwd=2)



#########################################################################################
# MIXED GAUSSIAN KERNEL & KNN LOCAL LINEAR REGRESSION WITH OPTIMAL NUMBER OF NEIGHBORS

mod.knn.opt <- npreg(flow.norm ~ scores.pc1, ckertype = 'gaussian', bws=opt.neigh, bwtype='adaptive_nn')
summary(mod.knn)

xx <- seq(min(scores.pc1), max(scores.pc1), length.out = 1000)
preds <- predict(mod.knn.opt, newdata = data.frame(scores.pc1 = xx), se.fit=T)
x11()
plot(scores.pc1, flow.norm, main = 'Gaussian kernel - KNN regression')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)

# RMSE
sqrt(mod.knn.opt$MSE)

graphics.off()
