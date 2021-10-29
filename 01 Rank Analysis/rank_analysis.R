rm(list=ls())

tab_flows <- read.table('trips_ext.csv', header=T, sep=';')

out_flow <- tab_flows$tot_out
in_flow <- tab_flows$tot_in


x11()
plot(in_flow, out_flow, main = 'Inflow vs Outflow', pch=19)

library(DepthProc)
tab <- data.frame(log_inflow=log(in_flow), log_outflow=log(out_flow))
tab <- data.frame(inflow=in_flow, outflow=out_flow)

load('mcshapiro.test.RData')
mcshapiro.test(tab)$pvalue
# data look very gaussian if I log trasform them, maha measure looks more appropriate
# data are not gaussian at all if I don't log trasform them, Tukey looks better


#######################################################################
# Tuckey Depth Measure
tukey_depth=depth(u=tab, method='Tukey')    #find depth measure according to Tukey method

#find median
median.tuckey <- depthMedian(tab, depth_params = list(method='Tukey')) 
median.country.tuckey <- tab_flows$Country[which.max(tukey_depth)] 
# without log trasf: CZ
# with log trasf: IE

x11()
plot(tab, pch=16, main='log(Inflow) vs log(Outflow)')
points(tab.median[1], tab.median[2], pch='*', col='red', cex=3)

x11()
depthContour(tab, depth_params = list(method='Tukey'), points=TRUE) 


#######################################################################
# Mahalanobis Depth measure
maha_depth=depth(u=tab, method='Mahalanobis')    #find depth measure according to Tukey method

#find median
median.maha <- depthMedian(tab, depth_params = list(method='Mahalanobis')) 
median.country.maha <- tab_flows$Country[which.max(tukey_depth)] 
# without log trasf: CZ
# with log trasf: IE

x11()
plot(tab, pch=16, main='log(Inflow) vs log(Outflow)')
points(tab.median[1], tab.median[2], pch='*', col='red', cex=3)

x11()
depthContour(tab, depth_params = list(method='Mahalanobis'), points=TRUE) 


########################################################################
#BAGPLOT AND OUTLIER IDENTIFICATION
library(aplpack)

x11()
bagplot(tab, show.whiskers = F, main="Bagplot", pch = '*', cex = 2)  

#outliers
tab.out <- bagplot(tab, create.plot = F)$pxy.outlier
colnames(tab.out) <- colnames(tab)
# with log trasf: there are no outliers
# without log trasf: GB


##################################################################################
#DDPLOT
# nel caso log trasf non ha senso, no outliers
# nel caso no log trasf non ha senso, solo un outlier

outliers.idx <- which(apply(tab, 1, function(x) all(x %in% tab.out)))
tab.clean <- tab[-outliers.idx,]


x11()
ddPlot(x = tab.clean, y = tab.out, depth_params = list(method='Tukey'))
x11()
ddPlot(x = tab, y = tab, depth_params = list(method='Mahalanobis'))

