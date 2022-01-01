#Here we apply the Ad Hoc ANOVA test that we designed in order to study the structure of the groups we identified

#Here is the definition of the test
#====

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


#====

Trips_labeled <- read_excel("Trips_labeled.xlsx")

my_col=c('red', 'blue')
plot(Trips_labeled$tot_in, Trips_labeled$tot_out, pch=19, main='Clustering Structure',
     log='xy', xlab='InFlow', ylab='OutFlow', col=my_col[as.factor(Trips_labeled$Mob_lab)])
legend('topleft', legend = c('High Mobility', 'Low Mobility'), fill = my_col)

B = 1000

anova_combination(data = Trips_labeled$tot_in, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                     Location   Scale                                                                                      
#Rank-correction          0      0.03
#Fisher-correction        0      0.03

anova_combination(data = Trips_labeled$tot_out, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                     Location   Scale                                                                                      
#Rank-correction          0       0
#Fisher-correction        0       0

anova_combination(data = Trips_labeled$tot_in + Trips_labeled$tot_out, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                     Location   Scale                                                                                      
#Rank-correction          0       0
#Fisher-correction        0       0

#indeed
x11()
boxplot(Trips_labeled$tot_in+Trips_labeled$tot_out ~ Trips_labeled$Mob_lab)

#_________________________________________________________________________________
anova_combination(data = Trips_labeled$tot_in, clust = Trips_labeled$EU_lab, B = B, output = F)
#                     Location   Scale                                                                                      
#Rank-correction          1      0.925
#Fisher-correction        1      0.996

anova_combination(data = Trips_labeled$tot_out, clust = Trips_labeled$EU_lab, B = B, output = F)
#                     Location   Scale                                                                                      
#Rank-correction        0.894    0.757
#Fisher-correction      0.939    0.939
anova_combination(data = Trips_labeled$tot_in + Trips_labeled$tot_out, clust = Trips_labeled$EU_lab, B = B, output = F)
#                     Location   Scale                                                                                      
#Rank-correction          1      0.898
#Fisher-correction        1      0.995

#indeed
x11()
boxplot(Trips_labeled$tot_in+Trips_labeled$tot_out ~ Trips_labeled$EU_lab)
summary(as.factor(Trips_labeled$EU_lab))

#====

#So we can see that the grouping structure defined by hierarchical clustering indeed
#separates high mobility countries from low mobility countries


#Let's see if this distinction is reflected in countries features

#====
#full_tab_31_10_21 <- read.csv(".../full_tab_31_10_21.csv")
data <- full_tab_31_10_21

#First with our most used variable:
anova_combination(data = data$uni_score, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction        0     0.02
#Fisher-correction      0     0.02

#Indeed there is a significant difference


#and our other potential regressors
anova_combination(data = data$Citations.per.document, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0.100  0.388
#Fisher-correction    0.163  0.388

#no significant difference

#____
anova_combination(data = data$num_staff, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0.011  0.025
#Fisher-correction    0.011  0.025

#significant difference
x11()
boxplot(data$num_staff ~ Trips_labeled$Mob_lab)

#____
anova_combination(data = data$stud_per_staff, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0.945  0.736
#Fisher-correction    0.945  0.913

#no significant difference

#____
anova_combination(data = data$GDP, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0.097  0.92
#Fisher-correction    0.201  0.92

#no significant difference
x11()
boxplot(data$GDP ~ Trips_labeled$Mob_lab)

#___
anova_combination(data = data$Reasearch, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0.048  0.71
#Fisher-correction    0.142  0.71

#there is a significant difference
x11()
boxplot(data$Reasearch ~ Trips_labeled$Mob_lab)

#____
anova_combination(data = data$Education, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0.759  0.759
#Fisher-correction    0.452  0.452

#No significant difference

#____
anova_combination(data = data$LPPI, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0.045  0.5
#Fisher-correction    0.102  0.5

#significant difference
x11()
boxplot(data$LPPI ~ Trips_labeled$Mob_lab)

#____
anova_combination(data = data$HDI, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0.042  0.262
#Fisher-correction    0.059  0.262

#significant difference
x11()
boxplot(data$HDI ~ Trips_labeled$Mob_lab)

#____
anova_combination(data = data$GGGI, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0.178  0.495
#Fisher-correction    0.127  0.495

#no significant difference

#____
anova_combination(data = data$PS, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0.031  0.242
#Fisher-correction    0.031  0.242

#significant difference
x11()
boxplot(data$PS ~ Trips_labeled$Mob_lab)
#remember that in this indicator, the smaller the better!

#____
anova_combination(data = data$SA, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0.407  0.407
#Fisher-correction    0.201  0.253

#no significant difference

#____
pat <- data$Patent.applications..nonresidents + data$Patent.applications..residents
anova_combination(data = pat, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0.015  0.036
#Fisher-correction    0.015  0.036

#significant difference
x11()
boxplot(pat ~ Trips_labeled$Mob_lab)

#____
anova_combination(data = data$uni_score_norm, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0.005  0.005
#Fisher-correction    0.000  0.002

#significant difference
x11()
boxplot(data$uni_score_norm ~ Trips_labeled$Mob_lab)
#so the significance of the quality of universities is confirmed also on this normalized indicator

#____
anova_combination(data = data$num_ric, clust = Trips_labeled$Mob_lab, B = B, output = F)
#                   Location Scale                                                                                 
#Rank-correction      0      0.006
#Fisher-correction    0      0.006

#significant difference
x11()
boxplot(data$uni_score_norm ~ Trips_labeled$Mob_lab)
#so the significance of the quality of universities is confirmed also on this normalized indicator

#____

#====

#Hence we can say that the prototype of a High Mobility country is one with:
# - A good academic enviroment: high score of its universities, with have a lot of staff
# - Spends more than other in R&D, proportionally to its GDP
# - has a good quality of living: good public services, human developement index and purchasing power
# - has an active entrepreneurial enviroment: a high number of patent applications