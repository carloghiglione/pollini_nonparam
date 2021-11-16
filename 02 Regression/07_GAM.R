rm(list=ls())
library(ISLR2)
library(car)
library(mgcv)
library(rgl)
library(splines)
library(pbapply)

#################
# load data (I normalize the covariates)

full_tab <- read.csv('full_tab_31_10_21.csv')

flow <- full_tab$tot_in - full_tab$tot_out    # in - out
flow.norm <- flow / full_tab$num_ric

countries <- full_tab$Country

tab <- data.frame(flow.norm, scale(full_tab[, c(4,6,8,17)]))

# remove NA (Cipro e Lussemburgo)
tab <- tab[-c(4,19),]

x11()
scatterplotMatrix(tab)



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
mod.gam <- gam(flow.norm ~ stud_per_staff +
                           s(Citations.per.document, bs='cr') + 
                           s(Reasearch, bs='cr') +
                           s(uni_score_norm, bs='cr'), data = tab)
summary(mod.gam)


# diagnostic
x11()
par(mfrow=c(2,2))
boxplot(mod.gam$residuals, main='Boxplot of residuals')
plot(mod.gam$fitted.values, mod.gam$residuals, main='Residuals vs fitted values')
abline(h=0, col='red')
hist(mod.gam$residuals, breaks = 10, col = 'gray', main = 'Histogramof residuals')
qqnorm(mod.gam$residuals)
qqline(mod.gam$residuals, col='red')
shapiro.test(mod.gam$residuals)


# check the contribution of each regressor to the response
x11()
par(mfrow=c(1,3))
plot(mod.gam)
