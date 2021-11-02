#we want to have a family wise error rate control over testing the difference in center and scaling of our two distributions
#Here the first cluster is the one suspected to have higher center AND higher scale

anova_combination_unilat <- function(data, clust, B, seed = 17121998, prog = T){
  # H0_1: mediana_clust1 == mediana_clust2 vs H11
  # H0_2: MAD_clust1 == MAD_clust2 vs H12 (MAD = median absolute deviation, robust estimator for scaling)
  # H0_12: mediana_clust1 == mediana_clust2 ^ MAD_clust1 == MAD_clust2 vs H0_12
  clus1 <- which(clust == levels(as.factor(clust))[1])
  clus2 <- which(clust == levels(as.factor(clust))[2])
  
  data1 <- data[clus1] #CLUST 1 is the one with HIGHER scale and center
  data2 <- data[clus2]
  n1 <- length(data1)
  n2 <- length(data2)
  n <- n1 + n2
  set.seed(seed)
  #difference in median
  T1_0 <- median(data1) - median(data2)
  #ratio of median absolute deviation
  T2_0 <- median((data1 - median(data1))^2)/median((data2 - median(data2))^2) - 1
  
  med1 <- median(data1)
  med2 <- median(data2)
  medvec <- ifelse(clust == levels(as.factor(clust))[1], med1,med2)
  #maximum between sum of ranks and sum of ranks of absolute deviations in the first cluster
  T12_0 <- max(sum(rank(data)[clus1]), sum(rank((data-medvec)^2)[clus1]))
  
  if(prog){
    library(progress)
    pb = progress_bar$new(total = B)
    pb$tick(0)
    }
  T1 <- numeric(B)
  T2 <- numeric(B)
  T12 <- numeric(B)
  #under the respective H0 the test statistics are invariant wrt the permutation
  #of our observations between the two groups -> if H0 includes equal medians
  #of the "residuals" ->if H0 is equal scaling (regardless of the medians!)
  for(b in  1:B){
    perm <- sample(1:n)
    #base permutation scheme: permute observations (valid under H0_1 and H0_12)
    x_pool_med <- data[perm]  
    
    x1_p_m <- x_pool_med[1:n1]
    x2_p_m <- x_pool_med[(n1+1):n]
    #to test different scaling withouth the need for equal centers we permute the DEVIATIONS from the medians themselves (valid under H0_2)
    x_pool_scale <- medvec + (data - medvec)[perm]
    
    x1_p_s <- x_pool_scale[1:n1]
    x2_p_s <- x_pool_scale[(n1+1):n]
    
    
    T1[b] <- median(x1_p_m) - median(x2_p_m)
    
    T2[b] <- median((x1_p_s - median(x1_p_s))^2)/median((x2_p_s - median(x2_p_s))^2) - 1
    
    med1p <- median(x1_p_m)
    med2p <- median(x2_p_m)
    medvecp <- c(replicate(n1,med1p),replicate(n2,med2p))
    T12[b] <- max(sum(rank(x_pool_med)[1:n1]), sum(rank((x_pool_med-medvecp)^2)[1:n1]))
    if(prog)
      pb$tick()
  }
  
  p1 <- sum(T1>=T1_0)/B
  p2 <- sum(T2>=T2_0)/B
  p12 <- sum(T12>=T12_0)/B
  
  #and now we can adjust the p_values for H0_1 and H0_2 controlling the family wise error rate
  p1tilda <- max(p1,p12)
  p2tilda <- max(p2,p12)
  return(cbind(p1 = p1tilda, p2 = p2tilda))
}

#empirical checks:
B = 100
dat1 <- rnorm(mean = 0, sd = 1, n = 15)
dat2 <- rnorm(mean = 0, sd = 1, n = 25)
dat <- c(dat1,dat2)
clus <- c(replicate(15,1), replicate(25,2))

anova_combination_unilat(dat, clus, B)
anova_combination_bilat(dat, clus, B)
 C <- 500
 p <- cbind(numeric(C), numeric(C), numeric(C))
pp <- progress_bar$new(total = C)
pp$tick(0)
 for(i in 1:C){
   set.seed(i^2)
   dat1 <- rnorm(mean = 0, sd = 1, n = 15)
   dat2 <- rnorm(mean = 0, sd = 1, n = 25)
   pp$tick()
   dat <- c(dat1,dat2)
   p[i,] <- anova_combination_bilat(dat,clus,B)
 }



B = 1000
dat1 <- rnorm(mean = 1, sd = 1, n = 15)
dat2 <- rnorm(mean = 0, sd = 1, n = 25)
dat <- c(dat1,dat2)
clus <- c(replicate(15,1), replicate(25,2))

anova_combination_unilat(dat, clus, B)
anova_combination_bilat(dat, clus, B)


B = 10000
dat1 <- rnorm(mean = 0, sd = 2, n = 15)
dat2 <- rnorm(mean = 0, sd = 1, n = 25)
dat <- c(dat1,dat2)
clus <- c(replicate(15,1), replicate(25,2))

anova_combination_unilat(dat, clus, B)
anova_combination_bilat(dat, clus, B)

B = 1000
dat1 <- rnorm(mean = 1, sd = 2, n = 15)
dat2 <- rnorm(mean = 0, sd = 1, n = 25)
dat <- c(dat1,dat2)
clus <- c(replicate(15,1), replicate(25,2))

anova_combination_unilat(dat, clus, B)



B <- 7500
an <- aov(log(as.numeric(trips_EU$tot_out)) ~ as.factor(trips_EU$High_flow_label))
summary(an)
shapiro.test(an$residuals)
anova_combination_unilat(log(as.numeric(trips_EU$tot_out)), 1 - as.numeric(trips_EU$High_flow_label), B)
boxplot(log(as.numeric(trips_EU$tot_out)) ~ as.factor(trips_EU$High_flow_label))
ansari.test(log(as.numeric(trips_EU$tot_out)) ~ as.factor(trips_EU$High_flow_label))

anova_combination_bilat(log(as.numeric(trips_EU$tot_out)), 1 - as.numeric(trips_EU$High_flow_label), B)



#estimation of power over scaling
s.grid <- seq(from = 0, to = 1, length.out = 100)
cc <- progress_bar$new(total = 100*v)
cc$tick(0)
v <- 100
B <- 200
pp <- numeric(100)
for(j in 1:100){
  p <- numeric(v)
  for(i in 1:v){
    set.seed(i*100 + 15*sqrt(j^3))
    dat1 <- rnorm(mean = 0, sd = 1+s.grid[j], n = 15)
    dat2 <- rnorm(mean = 0, sd = 1, n = 25)
    dat <- c(dat1,dat2)
    clus <- c(replicate(15,1), replicate(25,2))
    
    p[i] <- anova_combination_bilat(dat, clus, B, prog = FALSE)[2]
    cc$tick()
  }
  pp[j] <- sum(p <=0.05)/v
}
