rm(list=ls())
library(fda)
library(roahd)


set.seed(1234)

inflow <- read.csv('norm_year_inflow.csv')
outflow <- read.csv('norm_year_outflow.csv')

full_tab <- read.csv('full_tab_31_10_21.csv')

n <- dim(inflow)[1]

# (eventually remove last year)
#inflow <- inflow[,1:10]
#outflow <- outflow[,1:10]


# fill NA with the average across the years
k <- which(is.na(inflow), arr.ind=TRUE)
inflow[k] <- rowMeans(inflow[,-1], na.rm=TRUE)[k[,1]]

k <- which(is.na(outflow), arr.ind=TRUE)
outflow[k] <- rowMeans(outflow[,-1], na.rm=TRUE)[k[,1]]



remove_outlier = T
if(remove_outlier){
  inflow <- inflow[-4,]     # tolgo biellorussia che ha valori troppo alti (anche se tenerlo non cambia)
  outflow <- outflow[-4,]
  n <- n-1
}



# time of analysis
tt <- 2006:2015


# plot data by country
x11()
par(mfrow=c(1,2))
matplot(tt, t(inflow[,-1]), type = 'l', main='Inflow', xlab='Year', ylab='Norm. number of people')
matplot(tt, t(outflow[,-1]), type = 'l', main='Outflow', xlab='Year', ylab='Norm. number of people')



# smooth the curves
basis <- create.bspline.basis(rangeval = range(tt), nbasis = 20)     #default order is 4 (degree = 3), cubic spline
inflow.fd <- Data2fd(y=as.matrix(t(inflow[,-1])), argvals = tt, basisobj = basis)
outflow.fd <- Data2fd(y=as.matrix(t(outflow[,-1])), argvals = tt, basisobj = basis)

x11()
par(mfrow=c(1,2))
plot.fd(inflow.fd, xlab='Year', ylab='Norm. number of people', main='ciao')
plot.fd(outflow.fd, xlab='Year', ylab='Norm. number of people')


# create functional data objects
tt.grid <- seq(2006, 2015, length.out = 1000)
inflow.eval <- t(eval.fd(tt.grid, inflow.fd))     # I transpose because in fData I need data by row
outflow.eval <- t(eval.fd(tt.grid, outflow.fd))   # I transpose because in fData I need data by row


fData.in <- fData(tt.grid, inflow.eval)
fData.out <- fData(tt.grid, outflow.eval)


# build bivariate inflow-outflow dataset
bivariate_data <- as.mfData(list(fData.in, fData.out))

x11()
plot(bivariate_data)


# compute Spearman Correlation index (absolute value)
SPC0 <- abs(cor_spearman(bivariate_data, ordering='MEI'))
SPC0


# perform the permutational test
B <- 1000
SPC <- numeric(B)
for(i in 1:B){
  outflow.eval.shuffled <- outflow.eval[sample(n),]                   # shuffle outflow functions
  fData.in.curr <- fData(tt.grid, inflow.eval)
  fData.out.curr <- fData(tt.grid, outflow.eval.shuffled)
  bivariate_data <- as.mfData(list(fData.in.curr, fData.out.curr))
  SPC[i] <- abs(cor_spearman(bivariate_data, ordering='MEI'))
}


x11()
hist(SPC, breaks = 50, xlim = c(0,1), main='Test on Spearman Correlation Index', col = 'gray', xlab='Spearman Correlation Index')
abline(v=SPC0, col='red', lwd=3)

pvalue <- sum(SPC>=SPC0)/B
pvalue



x11()
par(mfrow=c(1,2))
matplot(tt.grid, eval.fd(tt.grid, inflow.fd), xlab='Year', ylab='Normalized Flow', main='InFlow(t)', type = 'l', ylim=c(0, 0.015))
matplot(tt.grid, eval.fd(tt.grid, outflow.fd), xlab='Year', ylab='Normalized Flow', main='OutFlow(t)', type = 'l', ylim=c(0, 0.015))


x11()
hist(SPC, breaks = 50, xlim = c(0,1), main='Independence Test', col = 'gray', xlab='T')
abline(v=SPC0, col='red', lwd=3)