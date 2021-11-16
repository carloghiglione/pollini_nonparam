rm(list=ls())
library(car)
library(ISLR2)
library(np)
library(splines)
library(fda)

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



#######################################################################################
# NATURAL CUBIC B-SPLINES
# SET KNOTS ACCORDING TO EVENLY SPACED QUANTILES

# number of knots
n_knots <- 3

int.knots <- quantile(scores.pc1, seq(0, 1, length.out=n_knots+2))[3:(n_knots)]
ext.knots <- quantile(scores.pc1, seq(0, 1, length.out=n_knots+2))[c(2, n_knots+1)]

mod.nat.spline.3 <- lm(flow.norm ~ ns(scores.pc1, knots = int.knots, Boundary.knots = ext.knots))
summary(mod.nat.spline.3)

x11()
par(mfrow=c(2,2))
plot(mod.nat.spline.3)
shapiro.test(mod.nat.spline.3$residuals)

# I plot fitted values +- 2 * standard errors)
xx <- seq(min(scores.pc1), max(scores.pc1), length.out = 1000)
preds <- predict(mod.nat.spline.3, data.frame(scores.pc1 = xx), se.fit=T)
x11()
plot(scores.pc1, flow.norm, main = 'Natural Cubic splines')
lines(xx, preds$fit, col='red', lwd=2)
matlines(xx, cbind(preds$fit - 2*preds$se.fit , preds$fit + 2*preds$se.fit ), 
         lty = 2, col = 'red', lwd=2)
abline(v=c(int.knots, ext.knots), lty=2, col='gray')



#######################################################################################
# SMOOTHING CUBIC B-SPLINES 
# SELECT OPTIMAL LAMBDA WITH CROSS-VALIDATION

fit.s.spline.loocv <- smooth.spline(scores.pc1, flow.norm, cv=T)
fit.s.spline.gcv <- smooth.spline(scores.pc1, flow.norm, cv=F)

x11()
plot(scores.pc1, flow.norm, main = 'Smoothing cubic splines')
lines(fit.s.spline.loocv, col='red', lwd=2)
#lines(fit.s.spline.gcv, col='blue', lwd=2)
legend('bottomright', legend = c('LOOCV', 'GCV'), fill = c('red', 'blue'))
