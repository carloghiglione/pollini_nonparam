rm(list=ls())
library(car)

full_tab <- read.csv('full_tab_31_10_21.csv')

##################
# creo tabelle con variabili di interesse
#full_tab <- full_tab[-12,]    #tolgo dato outlier (Gran Bretagna)
#full_tab <- full_tab[-4,]    #tolgo dato outlier Cypro

flow <- full_tab$tot_in - full_tab$tot_out    # in - out
flow.norm <- flow / full_tab$num_ric

#Z <- data.frame(flow, full_tab[, c(4,6,8,17)])
Z <- data.frame(flow, full_tab[, c(4,6,8,17)])

countries <- full_tab$Country


##################
# visualize results
x11()
pairs(Z)


#################
# do PCA on the regressors
tab <- Z[,-1]
x11()
par(mfrow=c(1,2))
boxplot(scale(tab, center = T, scale = F), main='Boxplot of original data')
barplot(sapply(tab, sd)^2, las=2, main='Original Data Variables')
# better to normalize them

tab.std <- data.frame(scale(tab))
x11()
boxplot(tab.std, main = 'Boxplot of standardized data')

pca.std<- princomp(tab.std, scores = T)
print(summary(pca.std))


#visualize total var
x11()
par(mfrow=c(1,3))
boxplot(pca.std$scores, main = "Boxplot of PCA data")
barplot(pca.std$sdev^2/sum(pca.std$sdev^2), ylim = c(0,1), main = "Proportional variances of PCA data")
plot(cumsum(pca.std$sdev^2)/sum(pca.std$sdev^2), ylim = c(0,1), type = "b", main = "Cumulate variance of PCA data")
abline(h=0.8, col = 'red')

#study of loadings
load.pca.std <- pca.std$loadings[,1:4]
x11()
par(mfrow=c(1,2))
barplot(load.pca.std [,1], ylim=c(-1, 1), main = 'PC1 loadings', las=2)   #evidence on variables 1,3,4 , contrast with 2
barplot(load.pca.std [,2], ylim=c(-1, 1), main = 'PC2 loadings', las=2)   #evidence on variable 2

# collect scores of first PCA
scores.pca.std  <- pca.std$scores[,1]

x11()
par(mfrow=c(1,2))
plot(scores.pca.std, flow, main = "Flow vs pca")
plot(scores.pca.std, flow.norm, main = "Flow norm vs pca")


#########################
# linear model on flow
mod <- lm(flow ~ scores.pca.std + I(scores.pca.std^2) + I(scores.pca.std^3))
summary(mod)

x11()
par(mfrow=c(2,2))
plot(mod)
shapiro.test(mod$residuals)

x11()
b <- mod$coefficients
xx = seq(min(scores.pca.std), max(scores.pca.std), length = 100)
plot(scores.pca.std, flow, main = "Flow vs pca")
lines(xx, b[1] + b[2]*xx + b[3]*xx^2 + b[4]*xx^3, col='red')


############################
# linear model on normalized flow
flow.norm[4] <- NA    # tolgo Cipro
flow.norm[19] <- NA    # tolgo Lussemburgo

# mod.norm <- lm(flow.norm ~ scores.pca.std + I(scores.pca.std^2) + I(scores.pca.std^3))
mod.norm <- lm(flow.norm ~ scores.pca.std)   # è il migliore
summary(mod.norm)

x11()
par(mfrow=c(2,2))
plot(mod.norm)
shapiro.test(mod.norm$residuals)

x11()
b <- mod.norm$coefficients
xx = seq(min(scores.pca.std), max(scores.pca.std), length = 100)
plot(scores.pca.std, flow.norm, main = "Normalized Flow vs pca")
#lines(xx, b[1] + b[2]*xx + b[3]*xx^2 + b[4]*xx^3, col='red')
lines(xx, b[1] + b[2]*xx, col='red')

