rm(list=ls())

library(rgl)
library(MASS)
library(rgl)
library(DepthProc)
library(hexbin)
library(packagefinder)
library(aplpack)
library(robustbase)

set.seed(1234)

###################################################
# LOAD DATA
full_tab <- read.csv('full_tab_31_10_21.csv')

pat <- full_tab$Patent.applications..nonresidents + full_tab$Patent.applications..residents
tab <- data.frame(full_tab[, c(4,6,8,17)], pat=log(pat))

countries <- full_tab$Country
######################################################################################
#DEPTH MEASURES

tukey_depth <- depth(u=tab, method='Tukey')    #find depth measure according to Tukey method
hist(tukey_depth, breaks = 20)


#find the median (point with highest depth)
tab.median <- depthMedian(tab, depth_params = list(method='Tukey')) 

x11()
pairs(rbind(tab, tab.median), 
      pch = c(rep(19, dim(tab)[1]),8), 
      col = c(rep('black', dim(tab)[1]), 'red'),
      cex = c(rep(1, dim(tab)[1]),2))



#########################################################################################
#BAGPLOT AND OUTLIER IDENTIFICATION

x11()
bagplot.pairs(tab)


#find the outliers according to bagplot doing pairwise comparisons of covariates
col_names <- colnames(tab)
outliers <- NULL

for(i in 1:4){
  for(j in (i+1):5){
    pxy_out <- bagplot(tab[,c(i,j)], create.plot = F)$pxy.outlier
    if(!is.null(pxy_out)){
      res <- apply(tab[,c(i,j)], MARGIN = 1, FUN = function(x){sum((x-pxy_out)^2)})
      outl <- which(res<1e-4)
      out_count <- as.character(countries[outl])
      outliers <- rbind(outliers, cbind(col_names[i], col_names[j],out_count))
    }else{
      outliers <- rbind(outliers, cbind(col_names[i], col_names[j], 'NaN'))
    }
  }
}
outliers





