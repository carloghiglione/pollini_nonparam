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
# NATURAL CUBIC B-SPLINES
# SET KNOTS ACCORDING TO EVENLY SPACED QUANTILES
# find optimal number of knots with AIC


# set the grid of dof search
knots.grid <- seq(3, 10, by=1)

# function to find AIC for a certain bandwidth
find.AIC <- function(curr.knots){
  int.knots <- quantile(scores.pc1, seq(0, 1, length.out=curr.knots+2))[3:(curr.knots)]
  ext.knots <- quantile(scores.pc1, seq(0, 1, length.out=curr.knots+2))[c(2, curr.knots+1)]
  mod.curr <- lm(flow.norm ~ ns(scores.pc1, knots = int.knots, Boundary.knots = ext.knots))
  return(extractAIC(mod.curr)[2]) 
}


# define cores for parallel computation
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = list("flow.norm", "scores.pc1","find.AIC", "ns"))
AIC_wrapper <- function(curr.knots){find.AIC(curr.knots)} 

# find AIC for all grid points with parallel computing
all_AIC <- pbsapply(knots.grid, AIC_wrapper, cl=cl)


# optimal number of neighbors, span and corresponding RMSE
opt.knots.AIC <- knots.grid[which.min(all_AIC)]
opt.knots.AIC
min(all_AIC)

x11()
plot(knots.grid, all_AIC, type ='l', lwd='2', main='AIC vs #knots')
abline(v=opt.knots.AIC, col='red', lwd=2)


int.knots <- quantile(scores.pc1, seq(0, 1, length.out=opt.knots.AIC+2))[3:(opt.knots.AIC)]
ext.knots <- quantile(scores.pc1, seq(0, 1, length.out=opt.knots.AIC+2))[c(2, opt.knots.AIC+1)]

mod.nat.3.opt.AIC <- lm(flow.norm ~ ns(scores.pc1, knots = int.knots, Boundary.knots = ext.knots))
summary(mod.nat.3.opt.AIC)

x11()
par(mfrow=c(2,2))
plot(mod.nat.3.opt.AIC)
shapiro.test(mod.nat.3.opt.AIC$residuals)

# I plot fitted values +- 2 * standard errors)
xx <- seq(min(scores.pc1), max(scores.pc1), length.out = 1000)
preds <- predict(mod.nat.3.opt.AIC, data.frame(scores.pc1 = xx), se.fit=T)
x11()
plot(scores.pc1, flow.norm, main = 'Natural Cubic splines, optimal AIC')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)
abline(v=c(int.knots, ext.knots), lty=2, col='gray')


RMSE.AIC <- sqrt(mean((mod.nat.3.opt.AIC$residuals)^2))
RMSE.AIC



############################################################################################
# NATURAL CUBIC B-SPLINES
# I find optimal dof with RMSE computed on STRATIFIED K-FOLD CROSS-VALIDATION

source('quantile_kfold.R')

# define the folds
n_fold <- 5
set.seed(1)
data.fr <- data.frame(flow.norm=flow.norm, scores.pc1 = scores.pc1)
kfolds <- quantile_kfold(n_fold, data.fr, 'scores.pc1')

# function to find RMSE for a certain bandwidth
find.RMSE <- function(curr.knots, n_fold, kfolds){
  require(np)
  RMSE <- numeric(n_fold)
  for(i in 1:n_fold){
    curr_data <- kfolds$datasets[[i]]
    curr_miss <- kfolds$folds[[i]]
    int.knots <- quantile(curr_data$scores.pc1, seq(0, 1, length.out=curr.knots+2))[3:(curr.knots)]
    ext.knots <- quantile(curr_data$scores.pc1, seq(0, 1, length.out=curr.knots+2))[c(2, curr.knots+1)]
    mod.curr <- lm(flow.norm ~ ns(scores.pc1, knots = int.knots, Boundary.knots = ext.knots), data = curr_data)
    RMSE[i] <- sqrt(mean((predict(mod.curr, newdata=curr_miss) - curr_miss$flow.norm)^2))
  }
  return(mean(RMSE)) 
}


# set the grid of dof search
knots.grid <- seq(3, 10, by=1)

# define cores for parallel computation
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = list("n_fold", "find.RMSE", "kfolds", "ns"))
RMSE_wrapper <- function(curr.knots){find.RMSE(curr.knots, n_fold, kfolds)} 

# find RMSE for all grid points with parallel computing
all_RMSE <- pbsapply(knots.grid, RMSE_wrapper, cl=cl)

# optimal number of neighbors, span and corresponding RMSE
opt.knots.RMSE <- knots.grid[which.min(all_RMSE)]
opt.knots.RMSE
min(all_RMSE)

x11()
plot(knots.grid, all_RMSE, type ='l', lwd='2', main='Root MSE vs n° knots', xlab = 'n° knots', ylab = 'RMSE')
abline(v=opt.knots.RMSE, col='red', lwd=2)



#######################################################################################
# NATURAL CUBIC B-SPLINES WITH OPTIMAL DOF


int.knots <- quantile(scores.pc1, seq(0, 1, length.out=opt.knots.RMSE+2))[3:(opt.knots)]
ext.knots <- quantile(scores.pc1, seq(0, 1, length.out=opt.knots.RMSE+2))[c(2, opt.knots+1)]

mod.nat.3.opt.RMSE <- lm(flow.norm ~ ns(scores.pc1, knots = int.knots, Boundary.knots = ext.knots))
summary(mod.nat.3.opt.RMSE)

x11()
par(mfrow=c(2,2))
plot(mod.nat.3.opt.RMSE)
shapiro.test(mod.nat.3.opt.RMSE$residuals)

# I plot fitted values +- 2 * standard errors)
xx <- seq(min(scores.pc1), max(scores.pc1), length.out = 1000)
preds <- predict(mod.nat.3.opt.RMSE, data.frame(scores.pc1 = xx), se.fit=T)
x11()
plot(scores.pc1, flow.norm, main = 'Natural Cubic splines, optimal RMSE')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)
abline(v=c(int.knots, ext.knots), lty=2, col='gray')


RMSE <- sqrt(mean((mod.nat.3.opt.RMSE$residuals)^2))
RMSE

extractAIC(mod.nat.3.opt.RMSE)


xx <- seq(min(scores.pc1), max(scores.pc1), length.out = 1000)
preds <- predict(mod.spline.3.opt.RMSE, data.frame(scores.pc1 = xx), se.fit=T)
y1 <- rev(preds$fit + 2*preds$se.fit)
y2 <- preds$fit - 2*preds$se.fit

x11()
plot(scores.pc1, flow.norm, main = 'Natural Cubic B-Spline', ylab='FlowNorm', xlab='CIndex', ylim=c(-0.008,0.014))
polygon(c(xx, rev(xx)), c(y2, y1),
        col = "gray", lty = 0)
points(scores.pc1, flow.norm)
lines(xx, preds$fit, col='black', lwd=2)
