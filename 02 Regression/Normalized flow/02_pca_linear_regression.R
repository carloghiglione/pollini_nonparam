rm(list=ls())
library(car)

full_tab <- read.csv('full_tab_31_10_21.csv')

##################
# creo tabelle con variabili di interesse

flow <- full_tab$tot_in - full_tab$tot_out    # in - out
flow.norm <- flow / full_tab$num_ric

pat <- full_tab$Patent.applications..nonresidents + full_tab$Patent.applications..residents
#Z <- data.frame(flow, full_tab[, c(4,6,8,17)])
Z <- data.frame(flow.norm, full_tab[, c(4,6,8,17)])
#Z <- data.frame(flow.norm, full_tab[, c(4,6,8,17)], pat=pat)

countries <- full_tab$Country

##################
# visualize results
x11()
pairs(Z)


##########################################################################################
# perform PCA on the regressors
tab <- Z[,-1]

x11()
par(mfrow=c(1,2))
boxplot(scale(tab, center = T, scale = F), main='Boxplot of original data')
barplot(sapply(tab, sd)^2, las=2, main='Original Data Variables')

# very significant difference among covariates, better to normalize them
tab.std <- data.frame(scale(tab))

x11()
boxplot(tab.std, main = 'Boxplot of standardized data')

# I do PCA
pca.std<- princomp(tab.std, scores = T)
print(summary(pca.std))

# visualize PCA performances
x11()
par(mfrow=c(1,3))
boxplot(pca.std$scores, main = "Boxplot of PCA data")
barplot(pca.std$sdev^2/sum(pca.std$sdev^2), ylim = c(0,1), main = "Proportional variances of PCA data")
plot(cumsum(pca.std$sdev^2)/sum(pca.std$sdev^2), ylim = c(0,1), type = "b", main = "Cumulate variance of PCA data")
abline(h=0.8, col = 'red')

# study of loadings
load.pca.std <- pca.std$loadings[,1:4]
x11()
par(mfrow=c(1,2))
barplot(load.pca.std [,1], ylim=c(-1, 1), main = 'PC1 loadings', las=2)   #evidence on variables 1,3,4 , contrast with 2
barplot(load.pca.std [,2], ylim=c(-1, 1), main = 'PC2 loadings', las=2)   #evidence on variable 2

# collect scores of first PCA
scores.pc1  <- pca.std$scores[,1]

x11()
par(mfrow=c(1,2))
plot(scores.pc1, flow, main = "Flow vs pca")
plot(scores.pc1, flow.norm, main = "Flow norm vs pca")



###################################################################################
# LINEAR REGRESSION ON FLOW
mod <- lm(flow ~ scores.pc1)
summary(mod)

x11()
par(mfrow=c(2,2))
plot(mod)
shapiro.test(mod$residuals)

x11()
b <- mod$coefficients
xx = seq(min(scores.pc1), max(scores.pc1), length = 100)
plot(scores.pc1, flow, main = "Flow vs pca")
lines(xx, b[1] + b[2]*xx, col='red')



#################################################################################
# LINEAR REGRESSION ON NORMALIZED FLOW

# remove outliers (4 Cipro e 19 Lussemburgo)
flow.norm <- flow.norm[-c(4,19)]
scores.pc1 <- scores.pc1[-c(4,19)]

mod.norm <- lm(flow.norm ~ scores.pc1)   # è il migliore
summary(mod.norm)

x11()
par(mfrow=c(2,2))
plot(mod.norm)
shapiro.test(mod.norm$residuals)

x11()
b <- mod.norm$coefficients
xx = seq(min(scores.pc1), max(scores.pc1), length = 100)
plot(scores.pc1, flow.norm, main = "Normalized Flow vs pca")
lines(xx, b[1] + b[2]*xx, col='red')

RMSE <- sqrt(mean(mod.norm$fitted.values^2))
RMSE



############################################################################################
# POLYNOMIAL REGRESSION
max_ord <- 8
mod.list <- lapply(1:max_ord, function(deg){lm(flow.norm ~ poly(scores.pc1, degree=deg))})

# I compare consecutive models (deg=i with deg=i+1) to discover the maximal order st it is significant
do.call(anova, mod.list)

# it resuts that linear model is sufficient



############################################################################################
# K-FOLD CROSS-VALIDATION
# I use it to estimate the Root-Mean-Square-Error

source('quantile_kfold.R')

n_fold <- 5
set.seed(1)
data.fr <- data.frame(flow.norm=flow.norm, scores.pc1 = scores.pc1)
kfolds <- quantile_kfold(n_fold, data.fr, 'scores.pc1')
RMSE <- numeric(n_fold)

for(i in 1:n_fold){
  curr_data <- kfolds$datasets[[i]]
  curr_miss <- kfolds$folds[[i]]
  mod.curr <- lm(flow.norm ~ scores.pc1, data = curr_data)
  
  RMSE[i] <- sqrt(mean((predict(mod.curr, curr_miss) - curr_miss$flow.norm)^2))
  
  x11()
  b <- mod.curr$coefficients
  xx = seq(min(scores.pc1), max(scores.pc1), length = 100)
  plot(curr_data$scores.pc1, curr_data$flow.norm, main = "Normalized Flow vs pca", 
       xlim = range(scores.pc1), ylim = range(flow.norm))
  lines(xx, b[1] + b[2]*xx, col='red')
  points(curr_miss$scores.pc1, curr_miss$flow.norm, col='blue', pch=19)
  points(curr_miss$scores.pc1, predict(mod.curr, curr_miss), col='blue', pch='x')
}
RMSE <- mean(RMSE)
RMSE

graphics.off()
