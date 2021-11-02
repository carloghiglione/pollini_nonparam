anova_combination_bilat <- function(data, clust, B, seed = 17121998, prog = TRUE){
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
  T1_0 <- (median(data1) - median(data2))^2
  #ratio of median absolute deviation
  T2_0 <- (median((data1 - median(data1))^2) - median((data2 - median(data2))^2))^2
  
  med1 <- median(data1)
  med2 <- median(data2)
  medvec <- ifelse(clust == levels(as.factor(clust))[1], med1,med2)
  #maximum between sum of ranks and sum of ranks of absolute deviations in the first cluster
  U1  <- sum(rank(data)[clus1]) - n1*(n1+1)/2
  U2  <- sum(rank(data)[clus2]) - n2*(n2+1)/2
  
  V1 <- sum(rank((data-medvec)^2)[clus1]) - n1*(n1+1)/2
  V2 <- sum(rank((data-medvec)^2)[clus2]) - n2*(n2+1)/2
  
  T12_0 <- max(max(U1,U2), max(V1, V2))
  
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
    
    
    T1[b] <- (median(x1_p_m) - median(x2_p_m))^2
    
    T2[b] <- (median((x1_p_s - median(x1_p_s))^2) - median((x2_p_s - median(x2_p_s))^2))^2
    
    U1_p  <- sum(rank(x_pool_med)[1:n1]) - n1*(n1+1)/2
    U2_p  <- sum(rank(x_pool_med)[(n1+1):n]) - n2*(n2+1)/2
    
    med1p <- median(x1_p_m)
    med2p <- median(x2_p_m)
    medvecp <- c(replicate(n1,med1p),replicate(n2,med2p))
    V1_p <- sum(rank((x_pool_med-medvecp)^2)[1:n1]) - n1*(n1+1)/2
    V2_p <- sum(rank((x_pool_med-medvecp)^2)[(n1+1):n]) - n2*(n2+1)/2
    
    T12[b] <- max(max(U1_p,U2_p), max(V1_p, V2_p))
    if(prog)
      pb$tick()
  }
  
  p1 <- sum(T1>=T1_0)/B
  p2 <- sum(T2>=T2_0)/B
  p12 <- sum(T12>=T12_0)/B
  
  #and now we can adjust the p_values for H0_1 and H0_2 controlling the family wise error rate
  p1tilda <- max(p1,p12)
  p2tilda <- max(p2,p12)
  return(cbind(center = p1tilda, scale = p2tilda, schifo = p12))
}

anova_combination_bilat(log(as.numeric(trips_EU$tot_out)), 1 - as.numeric(trips_EU$High_flow_label), B)
