#The goal is to determine wether there is an evolution of the impact of Uni_Score wrt time


#IN_FLOW
#data preparation
#====

year_inflow <- read.csv(".../year_inflow.csv", header = T)
full_tab_31_10_21 <- read.csv(".../full_tab_31_10_21.csv")

rnames <- year_inflow$country
cnames <- colnames(year_inflow)[-1]
flow_in <- as.vector(t(year_inflow[,-1]))
n <- length(data)
nc <- length(rnames)
ny <- length(cnames)
names <- c()
time <- as.vector(replicate(nc, cnames))
uni_score <- c()
k = 1
for (i in 1:nc){
  for (j in 1:ny){
    names[k] <- paste(rnames[i],cnames[j],sep = '')
    uni_score[k] <- full_tab_31_10_21$uni_score[which(full_tab_31_10_21$Country==rnames[i])]
    k = k+1
  }
}

time <- as.factor(time)
#====

#We first start by still contemplating a time-dependent intercept
#base model
g_base_in <- lm(flow_in ~ time + uni_score + uni_score:time)
summary(g_base_in)

#Under H0, the phenomenon is independent from time, up to an intercept
#We build a permutation test for this
g_H0_in <- lm(flow_in ~ time + uni_score)
summary(g_H0_in)

library(progress)
B = 10000
pb = progress_bar$new(total = B)
pb$tick(0)
set.seed(27011999)
fit_in <- g_H0_in$fitted.values
res_in <- g_H0_in$residuals
T_0_in <- sum(g_base_in$coefficients[12:20]^2)
T_stat_in <- numeric(B)
for(b in 1:B){
  perm <- numeric(n)
  for (i in 1:nc){
    perm[(1+(i-1)*ny):(ny+(i-1)*ny)] <- sample(ny)+(i-1)*ny
  }
  flow_pool <- fit_in + res_in[perm]
  g_perm <- lm(flow_pool ~ time + uni_score + uni_score:time)
  T_stat_in[b] <- sum(g_perm$coefficients[12:20]^2)
  pb$tick()
}
p_val_in <- sum(T_stat_in>=T_0_in)/B
p_val_in
# 0.8304 -> we cannot reject the null hypothesis

#We can even try to test the complet irrelevance of time factor
summary(g_base_in)

#Under H0, the phenomenon is independent from time
#We build a permutation test for this
g_H0_in_2 <- lm(flow_in ~ uni_score)
summary(g_H0_in_2)

library(progress)
B = 10000
pb = progress_bar$new(total = B)
pb$tick(0)
set.seed(27011999)
fit_in_2 <- g_H0_in_2$fitted.values
res_in_2 <- g_H0_in_2$residuals
T_0_in_2 <- sum(g_base_in$coefficients[-c(1,11)]^2)
T_stat_in_2 <- numeric(B)
for(b in 1:B){
  perm <- numeric(n)
  for (i in 1:nc){
    perm[(1+(i-1)*ny):(ny+(i-1)*ny)] <- sample(ny)+(i-1)*ny
  }
  flow_pool <- fit_in_2 + res_in_2[perm]
  g_perm <- lm(flow_pool ~ time + uni_score + uni_score:time)
  T_stat_in_2[b] <- sum(g_perm$coefficients[-11]^2)
  pb$tick()
}
p_val_in_2 <- sum(T_stat_in_2>=T_0_in_2)/B
p_val_in_2
# 1 -> we cannot reject the null hypothesis
#====

#we now do the same procedure for OUT_FLOW
#data preparation
year_outflow <- read.csv(".../year_outflow.csv", header=T)
flow_out <- as.vector(t(year_outflow[,-1]))


#We first start by still contemplating a time-dependent intercept
#base model
g_base_out <- lm(flow_out ~ time + uni_score + uni_score:time)
summary(g_base_out)

#Under H0, the phenomenon is independent from time, up to an intercept
#We build a permutation test for this
g_H0_out <- lm(flow_out ~ time + uni_score)
summary(g_H0_out)

library(progress)
B = 10000
pb = progress_bar$new(total = B)
pb$tick(0)
set.seed(27011999)
fit_out <- g_H0_out$fitted.values
res_out <- g_H0_out$residuals
T_0_out <- sum(g_base_out$coefficients[12:20]^2)
T_stat_out <- numeric(B)
for(b in 1:B){
  perm <- numeric(n)
  for (i in 1:nc){
    perm[(1+(i-1)*ny):(ny+(i-1)*ny)] <- sample(ny)+(i-1)*ny
  }
  flow_pool <- fit_out + res_out[perm]
  g_perm <- lm(flow_pool ~ time + uni_score + uni_score:time)
  T_stat_out[b] <- sum(g_perm$coefficients[12:20]^2)
  pb$tick()
}
p_val_out <- sum(T_stat_out>=T_0_out)/B
p_val_out
# 0.1776 -> we cannot reject the null hypothesis

#We can even try to test the complet irrelevance of time factor
summary(g_base)

#Under H0, the phenomenon is independent from time
#We build a permutation test for this
g_H0_out_2 <- lm(flow_out ~ uni_score)
summary(g_H0_out_2)

library(progress)
B = 10000
pb = progress_bar$new(total = B)
pb$tick(0)
set.seed(27011999)
fit_out_2 <- g_H0_out_2$fitted.values
res_out_2 <- g_H0_out_2$residuals
T_0_out_2 <- sum(g_base_out$coefficients[-c(1,11)]^2)
T_stat_out_2 <- numeric(B)
for(b in 1:B){
  perm <- numeric(n)
  for (i in 1:nc){
    perm[(1+(i-1)*ny):(ny+(i-1)*ny)] <- sample(ny)+(i-1)*ny
  }
  flow_pool <- fit_out_2 + res_out_2[perm]
  g_perm <- lm(flow_pool ~ time + uni_score + uni_score:time)
  T_stat_out_2[b] <- sum(g_perm$coefficients[-11]^2)
  pb$tick()
}
p_val_out_2 <- sum(T_stat_out_2>=T_0_out_2)/B
p_val_out_2
# 0.5946 -> we cannot reject the null hypothesis
#====
#Now for the total flow:
flow_tot = flow_in+flow_out
#We first start by still contemplating a time-dependent intercept
#base model
g_base_tot <- lm(flow_tot ~ time + uni_score + uni_score:time)
summary(g_base_tot)

#Under H0, the phenomenon is independent from time, up to an intercept
#We build a permutation test for this
g_H0_tot <- lm(flow_tot ~ time + uni_score)
summary(g_H0_tot)

library(progress)
B = 10000
pb = progress_bar$new(total = B)
pb$tick(0)
set.seed(27011999)
fit_tot <- g_H0_tot$fitted.values
res_tot <- g_H0_tot$residuals
T_0_tot <- sum(g_base_tot$coefficients[12:20]^2)
T_stat_tot <- numeric(B)
for(b in 1:B){
  perm <- numeric(n)
  for (i in 1:nc){
    perm[(1+(i-1)*ny):(ny+(i-1)*ny)] <- sample(ny)+(i-1)*ny
  }
  flow_pool <- fit_tot + res_tot[perm]
  g_perm <- lm(flow_pool ~ time + uni_score + uni_score:time)
  T_stat_tot[b] <- sum(g_perm$coefficients[12:20]^2)
  pb$tick()
}
p_val_tot <- sum(T_stat_tot>=T_0_tot)/B
p_val_tot
# 1 -> we cannot reject the null hypothesis

#We can even try to test the complet irrelevance of time factor
summary(g_base_tot)

#Under H0, the phenomenon is independent from time
#We build a permutation test for this
g_H0_tot_2 <- lm(flow_tot ~ uni_score)
summary(g_H0_tot_2)

library(progress)
B = 10000
pb = progress_bar$new(total = B)
pb$tick(0)
set.seed(27011999)
fit_tot_2 <- g_H0_tot_2$fitted.values
res_tot_2 <- g_H0_tot_2$residuals
T_0_tot_2 <- sum(g_base_tot$coefficients[-11]^2)
T_stat_tot_2 <- numeric(B)
for(b in 1:B){
  perm <- numeric(n)
  for (i in 1:nc){
    perm[(1+(i-1)*ny):(ny+(i-1)*ny)] <- sample(ny)+(i-1)*ny
  }
  flow_pool <- fit_tot_2 + res_tot_2[perm]
  g_perm <- lm(flow_pool ~ time + uni_score + uni_score:time)
  T_stat_tot_2[b] <- sum(g_perm$coefficients[-11]^2)
  pb$tick()
}
p_val_tot_2 <- sum(T_stat_tot_2>=T_0_tot_2)/B
p_val_tot_2
# 1 -> we cannot reject the null hypothesis
#====
#lastly for the net flow

flow_net = flow_in-flow_out
#We first start by still contemplating a time-dependent intercept
#base model
g_base_net <- lm(flow_net ~ time + uni_score + uni_score:time)
summary(g_base_net)

#Under H0, the phenomenon is independent from time, up to an intercept
#We build a permutation test for this
g_H0_net <- lm(flow_net ~ time + uni_score)
summary(g_H0_net)

library(progress)
B = 10000
pb = progress_bar$new(total = B)
pb$tick(0)
set.seed(27011999)
fit_net <- g_H0_net$fitted.values
res_net <- g_H0_net$residuals
T_0_net <- sum(g_base_net$coefficients[12:20]^2)
T_stat_net <- numeric(B)
for(b in 1:B){
  perm <- numeric(n)
  for (i in 1:nc){
    perm[(1+(i-1)*ny):(ny+(i-1)*ny)] <- sample(ny)+(i-1)*ny
  }
  flow_pool <- fit_net + res_net[perm]
  g_perm <- lm(flow_pool ~ time + uni_score + uni_score:time)
  T_stat_net[b] <- sum(g_perm$coefficients[12:20]^2)
  pb$tick()
}
p_val_net <- sum(T_stat_net>=T_0_net)/B
p_val_net
# 0.8839 -> we cannot reject the null hypothesis

#We can even try to test the complet irrelevance of time factor
summary(g_base_net)

#Under H0, the phenomenon is independent from time
#We build a permutation test for this
g_H0_net_2 <- lm(flow_net ~ uni_score)
summary(g_H0_net_2)

library(progress)
B = 10000
pb = progress_bar$new(total = B)
pb$tick(0)
set.seed(27011999)
fit_net_2 <- g_H0_net_2$fitted.values
res_net_2 <- g_H0_net_2$residuals
T_0_net_2 <- sum(g_base_net$coefficients[-11]^2)
T_stat_net_2 <- numeric(B)
for(b in 1:B){
  perm <- numeric(n)
  for (i in 1:nc){
    perm[(1+(i-1)*ny):(ny+(i-1)*ny)] <- sample(ny)+(i-1)*ny
  }
  flow_pool <- fit_net_2 + res_net_2[perm]
  g_perm <- lm(flow_pool ~ time + uni_score + uni_score:time)
  T_stat_net_2[b] <- sum(g_perm$coefficients[-11]^2)
  pb$tick()
}
p_val_net_2 <- sum(T_stat_net_2>=T_0_net_2)/B
p_val_net_2
# 0.7045 -> we cannot reject the null hypothesis
#====

#In our test there seems to be no effect at all of time over this phenomenon, it should however be noted that:
#- we have a VERY limited time horizon
#- we are assuming that the researchers signed in ORCID represent an unbiased sample of the total research world, and the unbiasedness wrt time is a very strong assumption
#- this test does not account for autocorrelation of the time series, indeed it considers the random errors to be iid country-wise, this assumption too is very strong and is likely leading to a very weak test
