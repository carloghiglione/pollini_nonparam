rm(list=ls())
library(ISLR2)
library(car)
library(mgcv)
library(rgl)
library(splines)
library(pbapply)

#################
# load data

full_tab <- read.csv('full_tab_31_10_21.csv')

flow <- full_tab$tot_in - full_tab$tot_out    # in - out
flow.norm <- flow / full_tab$num_ric

countries <- full_tab$Country

Z <- data.frame(flow.norm, full_tab[, c(4,6,8,17)])

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
scores.pc2  <- pca.std$scores[,2]
scores.pc3  <- pca.std$scores[,3]

tab <- data.frame(flow.norm=flow.norm,
                  scores.pc1=scores.pc1,
                  scores.pc2=scores.pc2,
                  scores.pc3=scores.pc3)

# remove NA (Cipro e Lussemburgo)
tab <- tab[-c(4,19),]



################################################################################
# MULTIVARIATE LINEAR MODEL
mod.lin <- lm(flow.norm ~ ., data = tab)             # without interaction
summary(mod.lin)

x11()
par(mfrow=c(2,2))
plot(mod.lin)
shapiro.test(mod.lin$residuals)



################################################################################
# GAM WITH SMOOTHING CUBIC SPLINES
mod.gam <- gam(flow.norm ~ s(scores.pc1, bs='cr') + s(scores.pc2, bs='cr')+ s(scores.pc3, bs='cr'), data = tab)
summary(mod.gam)



# diagnostic
x11()
par(mfrow=c(2,2))
boxplot(mod.gam$residuals, main='Boxplot of residuals')
plot(mod.gam$fitted.values, mod.gam$residuals, main='Residuals vs fitted values')
abline(h=0, col='red')
hist(mod.gam$residuals, breaks = 10, col = 'gray', main = 'Histogram of residuals')
qqnorm(mod.gam$residuals)
qqline(mod.gam$residuals, col='red')
shapiro.test(mod.gam$residuals)


# check the contribution of each regressor to the response
x11()
par(mfrow=c(1,3))
plot(mod.gam)