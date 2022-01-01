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
flow.norm <- (flow.norm - min(flow.norm) )/(max(flow.norm) - min(flow.norm))

countries <- full_tab$Country

#tab <- data.frame(flow.norm, scale(full_tab[, c(4,6,8,17)]))

tab <- data.frame(flow.norm, full_tab[, c(4,6,8,17)])
colnames(tab) <- c('FlowNorm', "CitDoc", "StudStaff", "ResGDP", "UniScore")

# remove NA (Cipro e Lussemburgo)
tab <- tab[-c(4,19),]

x11()
scatterplotMatrix(tab)



################################################################################
# MULTIVARIATE LINEAR MODEL
mod.lin <- lm(FlowNorm ~ ., data = tab)             # without interaction
summary(mod.lin)

x11()
par(mfrow=c(2,2))
plot(mod.lin)
shapiro.test(mod.lin$residuals)



################################################################################
# GAM WITH SMOOTHING CUBIC SPLINES
mod.gam <- gam(FlowNorm ~ StudStaff +
                          s(CitDoc, bs='cr') + 
                          s(ResGDP, bs='cr') +
                          s(UniScore, bs='cr'), data = tab)
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

x11()
par(mfrow=c(1,3))
plot(mod.gam, select = 1, ylab = 'NormFlow', main='Component 1 GAM', shade = T, se=2)
plot(mod.gam, select = 2, ylab = 'NormFlow', main='Component 2 GAM', shade = T, se=2)
plot(mod.gam, select = 3, ylab = 'NormFlow', main='Component 3 GAM', shade = T, se=2)
