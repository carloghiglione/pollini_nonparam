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



#######################################################################################
# CUBIC B-SPLINES (ORDER = 3), SET THE DEGREES OF FREEDOM
# KNOTS ARE PLACED AT THE CENTRAL (DOF - ORDER) EVEN QUANTILES (EXCLUDING 0% AND 100%)

# order of splines
order.s <- 3

# degrees of freedom (knots = dof - orders)
dof <- 7

mod.spline.3 <- lm(flow.norm ~ bs(scores.pc1, df=dof, degree=order.s))
summary(mod.spline.3)

# extract the knots
used.knots <- attributes(bs(scores.pc1, df=dof, degree=order.s))$knots

x11()
par(mfrow=c(2,2))
plot(mod.spline.3)
shapiro.test(mod.spline.3$residuals)

# I plot fitted values +- 2 * standard errors)
xx <- seq(min(scores.pc1), max(scores.pc1), length.out = 1000)
preds <- predict(mod.spline.3, data.frame(scores.pc1 = xx), se.fit=T)
x11()
plot(scores.pc1, flow.norm, main = 'Cubic splines')
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
data.fr <- data.frame(flow.norm=flow.norm, scores.pc1 = scores.pc1)
kfolds <- quantile_kfold(n_fold, data.fr, 'scores.pc1')

# function to find RMSE for a certain bandwidth
find.RMSE <- function(curr.dof, n_fold, kfolds){
  require(np)
  RMSE <- numeric(n_fold)
  for(i in 1:n_fold){
    curr_data <- kfolds$datasets[[i]]
    curr_miss <- kfolds$folds[[i]]
    mod.curr <- lm(flow.norm ~ bs(scores.pc1, df=curr.dof, degree=3), data = curr_data)
    RMSE[i] <- sqrt(mean((predict(mod.curr, newdata=curr_miss) - curr_miss$flow.norm)^2))
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

mod.spline.3.opt <- lm(flow.norm ~ bs(scores.pc1, df=opt.dof, degree=order.s))
summary(mod.spline.3.opt)

# extract the knots
used.knots <- attributes(bs(scores.pc1, df=opt.dof, degree=order.s))$knots

x11()
par(mfrow=c(2,2))
plot(mod.spline.3.opt)
shapiro.test(mod.spline.3.opt$residuals)

# I plot fitted values +- 2 * standard errors)
xx <- seq(min(scores.pc1), max(scores.pc1), length.out = 1000)
preds <- predict(mod.spline.3.opt, data.frame(scores.pc1 = xx), se.fit=T)
x11()
plot(scores.pc1, flow.norm, main = 'Cubic splines')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)
abline(v=used.knots, lty=2, col='gray')

RMSE <- sqrt(mean((mod.spline.3.opt$residuals)^2))
RMSE
