#ANOVA
# H0: center(X_1) = center(X_2) vs H1
#we simultaneously conduct the analysis over the permutational f_test and over a robust version wrt outliers, based on the median

anova_nonparam_test <- function(x, lab, B, seed = 27011999){
  clus1 <- which(lab == levels(as.factor(lab))[1])
  clus2 <- which(lab == levels(as.factor(lab))[2])
  x1 <- x[clus1]
  x2 <- x[clus2]
  n1 <- length(x1)
  n2 <- length(x2)
  n <- n1 + n2
  
  #robust wrt outliers
  T0_rob <- (median(x1) - median(x2))^2
  fit_0 <- aov(x ~ lab)
  #permutational f.test, (asymptotically) optimal if assumptions are verified
  T0_f <- summary(fit_0)[[1]][1,4]
  
  T_rob <- numeric(B)
  T_f <- numeric(B)
  
  set.seed(seed)
  library(progress)
  pb <- progress_bar$new(total = B)
  pb$tick(0)
  for(i in 1:B){
    #under H0 the joint distribution of the dataset is invariant wrt permutations of the units
    perm <- sample(1:n)
    x_pool <- x[perm]
    x1_p <- x_pool[clus1]
    x2_p <- x_pool[clus2]
    T_rob[i] <- (median(x1_p) - median(x2_p))^2
    fit_p <- aov(x_pool ~ lab)
    T_f[i] <- summary(fit_p)[[1]][1,4]
    pb$tick()
  }
  p_rob <- sum(T_rob >= T0_rob)/B
  p_f <- sum(T_f >= T0_f)/B
  return(cbind(p_robust = p_rob, p_f.test = p_f))
}
