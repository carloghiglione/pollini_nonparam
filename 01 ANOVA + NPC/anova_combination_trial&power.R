# H0_1: location_clust1 == location_clust2 vs H11 (test statistic : square difference of medians)
# H0_2: scale_clust1 == scale_clust2 vs H12 (test statistic: squeare difference of MAD)
# H0_12: location_clust1 == location_clust2 ^ scale_clust1 == scale_clust2 vs H0_12

#We are using Nonparametric Combination to simultaneously test for differences in location and scale between two populations
#the two parallel approaches we are following are:
# - using as test statistic for the joint test a handcrafted one, sensible to both violations
# - using a p-value combining function to automatically have the p-value of said test, without the need to define and compute a new statistic

#_____________________________________________________________________________________
#=====================================================================================
#_____________________________________________________________________________________

anova_combination <- function(data, clust, B, seed = 17121998, prog = TRUE, output = T){

  #first we get all data and quantities needed to define the test structure
  clus1 <- which(clust == levels(as.factor(clust))[1])
  clus2 <- which(clust == levels(as.factor(clust))[2])
  
  data1 <- data[clus1]
  data2 <- data[clus2]
  n1 <- length(data1)
  n2 <- length(data2)
  n <- n1 + n2
  set.seed(seed)
  
  
  #T_0 of the first test (difference in location)
  T1_0 <- (median(data1) - median(data2))^2
  
  #T_0 of the second test (difference in scale)
  T2_0 <- (median((data1 - median(data1))^2) - median((data2 - median(data2))^2))^2
  
  
  #in the first approach we defined ourselves a test statistic sensible to H0_12 under both violation 
  med1 <- median(data1)
  med2 <- median(data2)
  medvec <- ifelse(clust == levels(as.factor(clust))[1], med1,med2)
  #maximum between the test statistics of a rank test between observations and between absolute deviations, in each cluster
  U1  <- sum(rank(data)[clus1]) - n1*(n1+1)/2
  U2  <- sum(rank(data)[clus2]) - n2*(n2+1)/2
  
  V1 <- sum(rank((data-medvec)^2)[clus1]) - n1*(n1+1)/2
  V2 <- sum(rank((data-medvec)^2)[clus2]) - n2*(n2+1)/2

  #This way we can detect both violations by just taking the maximum, without having to worry about different scales of tests statistics based on difference of medians or difference of MADs
  #It's the maximum between two sums of ranks, always comparable
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
  # - of our observations between the two groups -> if H0 includes equal location
  # - of the "residuals" ->if H0 is equal scaling (regardless of the location!)
  for(b in  1:B){
    perm <- sample(1:n)
    #base permutation scheme: permute observations (valid under H0_1 and H0_12)
    x_pool_location <- data[perm]  
    
    x1_p_l <- x_pool_location[1:n1]
    x2_p_l <- x_pool_location[(n1+1):n]
    #to test different scaling withouth the need for equal centers we permute the DEVIATIONS from the medians themselves (valid under H0_2)
    x_pool_scale <- medvec + (data - medvec)[perm]
    
    x1_p_s <- x_pool_scale[clus1]
    x2_p_s <- x_pool_scale[clus2]
    
    #now we have the permuted datasets, so we can compute the new test statistics
    T1[b] <- (median(x1_p_l) - median(x2_p_l))^2
    
    T2[b] <- (median((x1_p_s - median(x1_p_s))^2) - median((x2_p_s - median(x2_p_s))^2))^2
    
    U1_p  <- sum(rank(x_pool_location)[1:n1]) - n1*(n1+1)/2
    U2_p  <- sum(rank(x_pool_location)[(n1+1):n]) - n2*(n2+1)/2
    
    med1p <- median(x1_p_l)
    med2p <- median(x2_p_l)
    medvecp <- c(replicate(n1,med1p),replicate(n2,med2p))
    V1_p <- sum(rank((x_pool_location-medvecp)^2)[1:n1]) - n1*(n1+1)/2
    V2_p <- sum(rank((x_pool_location-medvecp)^2)[(n1+1):n]) - n2*(n2+1)/2
    
    T12[b] <- max(max(U1_p,U2_p), max(V1_p, V2_p))
    if(prog)
      pb$tick()
  }
  
  p1 <- sum(T1>=T1_0)/B
  p2 <- sum(T2>=T2_0)/B
  p12 <- sum(T12>=T12_0)/B
  #and this is it for our first approach
  
  #another test statistic for the joint test is the fisher combination of the p-values of the first two tests
  #so
  
  Tglobal0 <- -2*(log(p1) + log(p2))
  
  #now we have to compute it for all permuted datasets
  p1stars <- numeric(B)
  p2stars <- numeric(B)
  Tglobals <- numeric(B)
  for(b in 1:B){
    p1stars[b] <- sum(T1>=T1[b])/B
    p2stars[b] <- sum(T2>=T2[b])/B
  }
  Tglobals <- -2*(log(p1stars) + log(p2stars))
  p12_fisher <- sum(Tglobals >= Tglobal0)/B
  
  #and now we can adjust the p_values for H0_1 and H0_2 controlling the family wise error rate using both approaches
  p1tilda <- max(p1,p12)
  p2tilda <- max(p2,p12)
  
  p1tildastar <- max(p1,p12_fisher)
  p2tildastar <- max(p2,p12_fisher)
  
  results <- cbind(c(p1tilda,p1tildastar), c(p2tilda, p2tildastar))
  rownames(results) = c("Rank-correction", "Fisher-correction")
  colnames(results) = c("Location", "Scale")
  if(output)
    print(results)
  return(results)
}


#_____________________________________________________________________________________
#=====================================================================================
#_____________________________________________________________________________________
#Trials
seed = 27011999

B = 1000
set.seed(seed)
dat1 <- rnorm(mean = 0, sd = 1.25, n = 50)
dat2 <- rnorm(mean = 0, sd = 1, n = 50)
data <- c(dat1,dat2)
clus <- c(replicate(50,1), replicate(50,2))

x11()
boxplot(data~clus)

tab = anova_combination(data, clus, B)
ansari.test(dat1,dat2)

#We tested multiple combinations at first, here is a more formal power test of the two approaches to detect a difference of scale given equal location
B = 100
delta_grid=seq(.25,2,by=.25)
alpha <- 0.05
clust <- c(replicate(50,1), replicate(50,2))

power.s.r=numeric(length(delta_grid))
power.s.f=numeric(length(delta_grid))
pb=progress_bar$new(total=B*length(delta_grid))
pb$tick(0)
for(ii in 1:length(delta_grid)){
  p.value.s.r <- numeric(B)
  p.value.s.f <- numeric(B)
  delta=delta_grid[ii]
  for(j in 1:B){
    set.seed(j*10 - sqrt(j))
    dat1 <- rnorm(mean = 0, sd = sqrt(1 + delta), n = 50)
    dat2 <- rnorm(mean = 0, sd = 1, n = 50)
    data <- c(dat1,dat2)
    tab <- anova_combination(data= data, clust = clust, B = 1000, prog = F, output = F)
    p.value.s.r[j] <- tab[1,2]
    p.value.s.f[j] <- tab[2,2]
    pb$tick()
  }
  
  estimated.power.r <- sum(p.value.s.r < alpha)/B
  power.s.r[ii]=estimated.power.r
  estimated.power.f <- sum(p.value.s.f < alpha)/B
  power.s.f[ii]=estimated.power.f
  
}      


#and the power to detect different locations given equal scales

power.l.r=numeric(length(delta_grid))
power.l.f=numeric(length(delta_grid))
pb=progress_bar$new(total=B*length(delta_grid))
pb$tick(0)

for(ii in 1:length(delta_grid)){
  p.value.s.r <- numeric(B)
  p.value.s.f <- numeric(B)
  delta=delta_grid[ii]
  for(j in 1:B){
    set.seed(j*10 - sqrt(j))
    dat1 <- rnorm(mean = 0+delta, sd = 1, n = 50)
    dat2 <- rnorm(mean = 0, sd = 1, n = 50)
    data <- c(dat1,dat2)
    tab <- anova_combination(data= data, clust = clust, B = 1000, prog = F, output = F)
    p.value.s.r[j] <- tab[1,1]
    p.value.s.f[j] <- tab[2,1]
    pb$tick()
  }
  
  estimated.power.r <- sum(p.value.s.r < alpha)/B
  power.l.r[ii]=estimated.power.r
  estimated.power.f <- sum(p.value.s.f < alpha)/B
  power.l.f[ii]=estimated.power.f
  
}


x11()
plot(delta_grid, power.l.f,ylim = c(0,1), col = 'red', type = 'o', ylab = 'estimated_power')
points(delta_grid, power.l.r, col = 'blue', type = 'o')
points(delta_grid, power.s.f, col = 'purple', type = 'o')
points(delta_grid, power.s.r, col = 'black', type = 'o')
points(delta_grid, power.ansari, col = 'green', type = 'o')
points(delta_grid, power.ansari.bonf, col = 'yellow', type = 'o')
legend(legend = c("location_fisher","location_rank","scale_fisher","scale_rank"), x = 'topleft', fill = c('red', 'blue', 'purple', 'black'))

power.l.f # 0.09 0.37 0.82 0.96 1.00 1.00 1.00 1.00
power.l.r #0.09 0.43 0.87 1.00 1.00 1.00 1.00 1.00
power.s.f #0.04 0.05 0.08 0.16 0.28 0.38 0.44 0.52
power.s.r #0.04 0.06 0.12 0.18 0.29 0.42 0.49 0.57


#_____________________________________________________________________________________
#=====================================================================================
#_____________________________________________________________________________________


# to get a better glance we can compare the less powerful "branch" with an analogous nonparametric test: ansari test
pb=progress_bar$new(total=B*length(delta_grid))
pb$tick(0)
power.ansari <- numeric(length(delta_grid))
power.ansari.bonf <- numeric(length(delta_grid))
for(ii in 1:length(delta_grid)){
  p.value <- numeric(B)
  delta=delta_grid[ii]
  for(j in 1:B){
    set.seed(j*10 - sqrt(j))
    dat1 <- rnorm(mean = 0, sd = sqrt(1+delta), n = 50)
    dat2 <- rnorm(mean = 0, sd = 1, n = 50)
    p.value[j] <- ansari.test(scale(dat1,scale = F,center = T) ,scale(dat2, scale = F, center = T))$p.value
    pb$tick()
  }
  
  estimated.power <- sum(p.value < alpha)/B
  estimated.power.bonf <- sum(p.value*2 < alpha)/B
  power.ansari[ii]=estimated.power
  power.ansari.bonf[ii]= estimated.power.bonf
}
power.ansari       # 0.08 0.15 0.26 0.39 0.48 0.63 0.77 0.83
power.ansari.bonf  # 0.05 0.11 0.16 0.27 0.39 0.48 0.60 0.72


x11()
plot(delta_grid, power.s.f,ylim = c(0,1), col = 'red', type = 'o', ylab = 'estimated_power')
points(delta_grid, power.s.r, col = 'blue', type = 'o')
points(delta_grid, power.ansari, col = 'green', type = 'o')
points(delta_grid, power.ansari.bonf, col = 'yellow', type = 'o')
legend(legend = c("scale_fisher","scale_rank", "ansari", "ansari with bonferroni correction"), x = 'topleft', fill = c('red', 'blue', 'purple', 'black'))

#so we can see that we do not lose much power, but in this setting we are still more conservative than Bonferroni

#Let's see however what happens if the assumption of the ansari test are not satisfied and we indeed are without the same location

B = 100
delta_grid=seq(.25,2,by=.25)
alpha <- 0.05
clust <- c(replicate(50,1), replicate(50,2))

power.s.r_ld=numeric(length(delta_grid))
power.s.f_ld=numeric(length(delta_grid))
power.ansari_ld <- numeric(length(delta_grid))
power.ansari.bonf_ld <- numeric(length(delta_grid))
pb=progress_bar$new(total=B*length(delta_grid))
pb$tick(0)
for(ii in 1:length(delta_grid)){
  p.value.s.r <- numeric(B)
  p.value.s.f <- numeric(B)
  delta=delta_grid[ii]
  for(j in 1:B){
    set.seed(j*10 - sqrt(j))
    dat1 <- rnorm(mean = 1, sd = sqrt(1 + delta), n = 50)
    dat2 <- rnorm(mean = 0, sd = 1, n = 50)
    data <- c(dat1,dat2)
    tab <- anova_combination(data= data, clust = clust, B = 1000, prog = F, output = F)
    p.value.s.r[j] <- tab[1,2]
    p.value.s.f[j] <- tab[2,2]
    p.value.ans[j] <- ansari.test(scale(dat1, scale = F, center = T),scale(dat2, scale = F, center = T))$p.value
    pb$tick()
  }
  
  estimated.power.r <- sum(p.value.s.r < alpha)/B
  power.s.r_ld[ii]=estimated.power.r
  estimated.power.f <- sum(p.value.s.f < alpha)/B
  power.s.f_ld[ii]=estimated.power.f
  estimated.power.ans <- sum(p.value.ans < alpha)/B
  estimated.power.ans.bonf <- sum(p.value.ans*2 < alpha)/B
  power.ansari_ld[ii]=estimated.power.ans
  power.ansari.bonf_ld[ii]= estimated.power.ans.bonf
}

power.s.f_ld          # 0.04 0.10 0.19 0.30 0.38 0.51 0.56 0.64
power.s.r_ld          # 0.04 0.10 0.19 0.30 0.38 0.49 0.55 0.64
power.ansari_ld       # 0.08 0.15 0.26 0.39 0.48 0.63 0.77 0.83
power.ansari.bonf_ld  # 0.05 0.11 0.16 0.27 0.39 0.48 0.60 0.72

set.seed(seed)
jit <- rnorm(mean = 0, sd = 10^-2, n = length(delta_grid))

x11()
plot(delta_grid, power.s.f_ld+jit,ylim = c(0,1),main = 'N(1, 1+delta) vs N(0,1)', col = 'red', type = 'o', ylab = 'estimated_power')
points(delta_grid, power.s.r_ld, col = 'blue', type = 'o')
points(delta_grid, power.ansari_ld, col = 'purple', type = 'o')
points(delta_grid, power.ansari.bonf_ld, col = 'black', type = 'o')
legend(legend = c("scale_fisher","scale_rank", "ansari", "ansari with bonferroni correction"), x = 'topleft', fill = c('red', 'blue', 'purple', 'black'))
#_____________________________________________________________________________________
#=====================================================================================
#_____________________________________________________________________________________

#So we can see that this procedure is as conservative as applying Bonferroni on a quite powerful test (like the Ansari-Bradley test)
#Let's see if this procedure would be worse than a Bonferroni on the same permutation test

#=====================================================================================
basicperm_scale <- function(data, clust, B, seed = 17121998, prog = TRUE, output = T){
  
  #first we get all data and quantities needed to define the test structure
  clus1 <- which(clust == levels(as.factor(clust))[1])
  clus2 <- which(clust == levels(as.factor(clust))[2])
  
  data1 <- data[clus1]
  data2 <- data[clus2]
  n1 <- length(data1)
  n2 <- length(data2)
  n <- n1 + n2
  set.seed(seed)

  #T_0 of the difference in scale test
  T_0 <- (median((data1 - median(data1))^2) - median((data2 - median(data2))^2))^2
  
  med1 <- median(data1)
  med2 <- median(data2)
  medvec <- ifelse(clust == levels(as.factor(clust))[1], med1,med2)
  
  if(prog){
    library(progress)
    pb = progress_bar$new(total = B)
    pb$tick(0)
  }
  T_perm <- numeric(B)
  for(b in  1:B){
    perm <- sample(1:n)
    
    x_pool_scale <- medvec + (data - medvec)[perm]
    
    x1_p_s <- x_pool_scale[clus1]
    x2_p_s <- x_pool_scale[clus2]
    
    
    T_perm[b] <- (median((x1_p_s - median(x1_p_s))^2) - median((x2_p_s - median(x2_p_s))^2))^2
    
    if(prog)
      pb$tick()
  }
  
  pval <- sum(T_perm>=T_0)/B
  return(pval)
}

#=====================================================================================
power.NPC.r=numeric(length(delta_grid))
power.NPC.f=numeric(length(delta_grid))
power.perm=numeric(length(delta_grid))
power.perm.bonf=numeric(length(delta_grid))
pb=progress_bar$new(total=B*length(delta_grid))
pb$tick(0)

for(ii in 1:length(delta_grid)){
  p.value.NPC.r <- numeric(B)
  p.value.NPC.f <- numeric(B)
  p.value.perm <- numeric(B)
  delta=delta_grid[ii]
  for(j in 1:B){
    set.seed(ii*j*10 - sqrt(j))
    dat1 <- rnorm(mean = 0, sd = 1+delta, n = 50)
    dat2 <- rnorm(mean = 0, sd = 1, n = 50)
    data <- c(dat1,dat2)
    tab <- anova_combination(data= data, clust = clust, B = 1000, prog = F, output = F)
    p.value.NPC.r[j] <- tab[1,2]
    p.value.NPC.f[j] <- tab[2,2]
    p.value.perm[j] <- basicperm_scale(data= data, clust = clust, B = 1000, prog = F, output = F)
    pb$tick()
  }
  estimated.power.r <- sum(p.value.NPC.r < alpha)/B
  power.NPC.r[ii]=estimated.power.r
  estimated.power.f <- sum(p.value.NPC.f < alpha)/B
  power.NPC.f[ii]=estimated.power.f
  estimated.power.perm <- sum(p.value.perm < alpha)/B
  power.perm[ii]=estimated.power.perm
  estimated.power.perm.b <- sum(2*p.value.perm < alpha)/B
  power.perm.bonf[ii]=estimated.power.perm.b
}

x11()
plot(delta_grid, power.NPC.f,ylim = c(0,1), col = 'red', type = 'o', ylab = 'estimated_power', main = 'N(0,(1+delta)^2) vs N(0,1)')
points(delta_grid, power.NPC.r, col = 'blue', type = 'o')
points(delta_grid, power.perm, col = 'purple', type = 'o')
points(delta_grid, power.perm.bonf, col = 'black', type = 'o')
legend(legend = c("scale_fisher","scale_rank","MAD test","MAD test Bonf"), x = 'topleft', fill = c('red', 'blue', 'purple', 'black'))

power.NPC.r      # 0.07 0.44 0.66 0.86 0.95 0.99 0.99 0.99
power.NPC.f      # 0.06 0.37 0.59 0.82 0.90 0.98 0.99 0.99
power.perm       # 0.12 0.48 0.70 0.88 0.95 0.99 0.99 0.99
power.perm.bonf  # 0.05 0.40 0.57 0.83 0.89 0.97 0.99 0.99

#___________________________________________________________________________________
power.NPC.r.d=numeric(length(delta_grid))
power.NPC.f.d=numeric(length(delta_grid))
power.perm.d=numeric(length(delta_grid))
power.perm.bonf.d=numeric(length(delta_grid))
pb=progress_bar$new(total=B*length(delta_grid))
pb$tick(0)

for(ii in 1:length(delta_grid)){
  p.value.NPC.r <- numeric(B)
  p.value.NPC.f <- numeric(B)
  p.value.perm <- numeric(B)
  delta=delta_grid[ii]
  for(j in 1:B){
    set.seed(ii*j*10 - sqrt(j))
    dat1 <- rnorm(mean = 0, sd = 1+delta, n = 50)
    dat2 <- rnorm(mean = 1, sd = 1, n = 50)
    data <- c(dat1,dat2)
    tab <- anova_combination(data= data, clust = clust, B = 1000, prog = F, output = F)
    p.value.NPC.r[j] <- tab[1,2]
    p.value.NPC.f[j] <- tab[2,2]
    p.value.perm[j] <- basicperm_scale(data= data, clust = clust, B = 1000, prog = F, output = F)
    pb$tick()
  }
  estimated.power.r <- sum(p.value.NPC.r < alpha)/B
  power.NPC.r.d[ii]=estimated.power.r
  estimated.power.f <- sum(p.value.NPC.f < alpha)/B
  power.NPC.f.d[ii]=estimated.power.f
  estimated.power.perm <- sum(p.value.perm < alpha)/B
  power.perm.d[ii]=estimated.power.perm
  estimated.power.perm.b <- sum(2*p.value.perm < alpha)/B
  power.perm.bonf.d[ii]=estimated.power.perm.b
}

x11()
plot(delta_grid, power.NPC.f.d,ylim = c(0,1), col = 'red', type = 'o', ylab = 'estimated_power', main = 'N(0,(1+delta)^2) vs N(1,1)')
points(delta_grid, power.NPC.r.d, col = 'blue', type = 'o')
points(delta_grid, power.perm.d, col = 'purple', type = 'o')
points(delta_grid, power.perm.bonf.d, col = 'black', type = 'o')
legend(legend = c("scale_fisher","scale_rank","MAD test","MAD test Bonf"), x = 'topleft', fill = c('red', 'blue', 'purple', 'black'))

power.NPC.r.d      # 0.12 0.48 0.69 0.88 0.95 0.99 0.99 0.99
power.NPC.f.d      # 0.12 0.48 0.68 0.87 0.94 0.98 0.99 0.99
power.perm.d       # 0.12 0.48 0.70 0.88 0.95 0.99 0.99 0.99
power.perm.bonf.d  # 0.05 0.40 0.57 0.83 0.89 0.97 0.99 0.99

#___________________________________________________________________________________
power.NPC.r.exp=numeric(length(delta_grid))
power.NPC.f.exp=numeric(length(delta_grid))
power.perm.exp=numeric(length(delta_grid))
power.perm.bonf.exp=numeric(length(delta_grid))
pb=progress_bar$new(total=B*length(delta_grid))
pb$tick(0)

for(ii in 1:length(delta_grid)){
  p.value.NPC.r <- numeric(B)
  p.value.NPC.f <- numeric(B)
  p.value.perm <- numeric(B)
  delta=delta_grid[ii]
  for(j in 1:B){
    set.seed(ii*j*10 - sqrt(j))
    dat1 <- rexp(1/(1+delta), n = 50)
    dat2 <- rexp(1, n = 50)
    data <- c(dat1,dat2)
    tab <- anova_combination(data= data, clust = clust, B = 1000, prog = F, output = F)
    p.value.NPC.r[j] <- tab[1,2]
    p.value.NPC.f[j] <- tab[2,2]
    p.value.perm[j] <- basicperm_scale(data= data, clust = clust, B = 1000, prog = F, output = F)
    pb$tick()
  }
  estimated.power.r <- sum(p.value.NPC.r < alpha)/B
  power.NPC.r.exp[ii]=estimated.power.r
  estimated.power.f <- sum(p.value.NPC.f < alpha)/B
  power.NPC.f.exp[ii]=estimated.power.f
  estimated.power.perm <- sum(p.value.perm < alpha)/B
  power.perm.exp[ii]=estimated.power.perm
  estimated.power.perm.b <- sum(2*p.value.perm < alpha)/B
  power.perm.bonf.exp[ii]=estimated.power.perm.b
}

x11()
plot(delta_grid, power.NPC.f.exp,ylim = c(0,1), col = 'red', type = 'o', ylab = 'estimated_power', main = 'Exp(1/(1+delta)) vs Exp(1)')
points(delta_grid, power.NPC.r.exp, col = 'blue', type = 'o')
points(delta_grid, power.perm.exp, col = 'purple', type = 'o')
points(delta_grid, power.perm.bonf.exp, col = 'black', type = 'o')
legend(legend = c("scale_fisher","scale_rank","MAD test","MAD test Bonf"), x = 'topleft', fill = c('red', 'blue', 'purple', 'black'))

power.NPC.r.exp      # 0.14 0.27 0.47 0.66 0.80 0.94 0.96 0.95
power.NPC.f.exp      # 0.17 0.38 0.51 0.72 0.80 0.93 0.95 0.96
power.perm.exp       # 0.21 0.39 0.53 0.72 0.82 0.94 0.96 0.96
power.perm.bonf.exp  # 0.18 0.34 0.44 0.61 0.73 0.92 0.91 0.93

#___________________________________________________________________________________
power.NPC.r.csq=numeric(length(delta_grid))
power.NPC.f.csq=numeric(length(delta_grid))
power.perm.csq=numeric(length(delta_grid))
power.perm.bonf.csq=numeric(length(delta_grid))
pb=progress_bar$new(total=B*length(delta_grid))
pb$tick(0)

for(ii in 1:length(delta_grid)){
  p.value.NPC.r <- numeric(B)
  p.value.NPC.f <- numeric(B)
  p.value.perm <- numeric(B)
  delta=delta_grid[ii]
  for(j in 1:B){
    set.seed(ii*j*10 - sqrt(j))
    dat1 <- rchisq(5+4*delta, n = 50)
    dat2 <- rchisq(5, n = 50)
    data <- c(dat1,dat2)
    tab <- anova_combination(data= data, clust = clust, B = 1000, prog = F, output = F)
    p.value.NPC.r[j] <- tab[1,2]
    p.value.NPC.f[j] <- tab[2,2]
    p.value.perm[j] <- basicperm_scale(data= data, clust = clust, B = 1000, prog = F, output = F)
    pb$tick()
  }
  estimated.power.r <- sum(p.value.NPC.r < alpha)/B
  power.NPC.r.csq[ii]=estimated.power.r
  estimated.power.f <- sum(p.value.NPC.f < alpha)/B
  power.NPC.f.csq[ii]=estimated.power.f
  estimated.power.perm <- sum(p.value.perm < alpha)/B
  power.perm.csq[ii]=estimated.power.perm
  estimated.power.perm.b <- sum(2*p.value.perm < alpha)/B
  power.perm.bonf.csq[ii]=estimated.power.perm.b
}

set.seed(seed)
jit1 <- rnorm(mean = 0, sd = 10^-3, n = length(delta_grid))
jit2 <- rnorm(mean = 0, sd = 10^-3, n = length(delta_grid))

x11()
plot(4*delta_grid, power.NPC.f.csq + jit1,ylim = c(0,1), col = 'red', type = 'o', ylab = 'estimated_power', main = 'Chisq(5 + delta) vs Chisq(5')
points(4*delta_grid, power.NPC.r.csq + jit2, col = 'blue', type = 'o')
points(4*delta_grid, power.perm.csq, col = 'purple', type = 'o')
points(4*delta_grid, power.perm.bonf.csq, col = 'black', type = 'o')
legend(legend = c("scale_fisher","scale_rank","MAD test","MAD test Bonf"), x = 'topleft', fill = c('red', 'blue', 'purple', 'black'))

power.NPC.r.csq       # 0.04 0.20 0.28 0.37 0.41 0.45 0.59 0.67
power.NPC.f.csq       # 0.04 0.20 0.28 0.37 0.41 0.45 0.59 0.67
power.perm.csq        # 0.04 0.20 0.28 0.37 0.41 0.45 0.59 0.67
power.perm.bonf.csq   # 0.03 0.15 0.18 0.22 0.31 0.36 0.49 0.58

#=================================================================================
#So we can see that this strategy, altought far from being perfect, allows us to have a control
#over the FWER, while being less conservative than a Bonferroni correction

#Our intent is to use this approach to account for the dependence among the location and scale tests,
#treating them as simultanoeous results; hence being able to identify a difference in the location of two groups
#accounting for the possible eteroschedasticity