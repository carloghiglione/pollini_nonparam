# Prepare dataframes
#===== 
setwd("~/Documents/Polimi/NONPARAMETRIC/pollini_nonparam/03 Analisi lungo il tempo/")
df.uni_score <- read.csv("../02 Regression/Total flow/full_tab_31_10_21.csv")
colnames(df)
df.uni_score <- df.uni_score[, c(1,18)]
df.uni_score

df <- read.csv("Flows_per_year.csv")

df.out <- aggregate(count~startC + year, df[,-2], FUN = sum)
df.out
df.in <- aggregate(count~endC + year, df[,-1], FUN = sum)
df.in
country <- unique(df$endC)
df.ok <- data.frame(country)

# assign values
for (y in unique(df.out$year)){
  print(y)
  df.ok <- merge(x=df.ok, y=df.out[df.out$year==y,c(1,3)], 
                 by.x = "country",  by.y = "startC")
  names(df.ok)[ncol(df.ok)] = paste(y, "-out") 
  df.ok <- merge(x=df.ok, y=df.in[df.in$year==y,c(1,3)], 
                 by.x = "country",  by.y = "endC")
  names(df.ok)[ncol(df.ok)] = paste(y, "-in") 
  # in - out ## TODO CHANGE
  df.ok$new =  df.ok[,ncol(df.ok)] - (df.ok[,ncol(df.ok)-1])
  names(df.ok)[ncol(df.ok)] = paste(y, "-total")
}
rownames(df.ok) <- df.ok$country
df.ok
df.ok.tot <- df.ok[df.ok$country %in% unique(df.out$startC),sapply(names(df.ok), grepl, pattern="total")]

df.ok.tot$Country <- rownames(df.ok.tot)

df.aux <- df.ok.tot[,11]

df.aux <- cbind(df.aux, as.numeric(df.ok.tot[,1]))

library(reshape)
df.re <- melt(df.ok.tot, id.vars="Country", var="time")

df.re$time <- as.factor(sapply(df.re$time, substr, start=1, stop=4))
colnames(df.re)[3] = "net_flow"

df.re <- merge(df.re, df.uni_score, by="Country")


###########
# Permutation test of mixed effects linear model
#====
library(mgcv)
library(nlme)
# H0 model is the one already fitted
# H1: interaction of uni_score with time is significant

# Test statistic
lmm.full <-  lme(fixed = net_flow ~ 1 + uni_score + uni_score:time,
                 random = ~ 1 + as.integer(time)|Country, # random effect ~N(time, Sigma)
                 correlation= NULL,  # e_i are iid , ~N(0,sigma^2 I)
                 data=df.re, method = "REML")
summ.full <- summary(lmm.full)
gam.full <- gam(log(net_flow) ~  s(uni_score, bs="cr") + 
                  uni_score:time + s(Country, bs="re"),
                data=df.re)

getVarCov(lmm.full, df.re$Country,  type="random.effects")
getVarCov(lmm.full, df.re$Country, type = "conditional") ## e_i term
getVarCov(lmm.full, "IT", type = "marginal")  # Random term
#getVarCov(lmm.full, df.re$Country, type = "marginal") 

# t value
T0 <- max(abs(summ.full$tTable[-c(1,2),4]))  # extract t value of uni_score:time

# reduced model (under H0)
lmm.red <-  lme(fixed = net_flow ~ 1 + uni_score, 
                random = ~ 1 + as.integer(time)|Country, # random effect ~N(time, Sigma)
                correlation= NULL,  # e_i are iid , ~N(0,sigma^2 I)
                data=df.re, method = "REML" )
getVarCov(lmm.red, df.re$Country,  type="random.effects")
getVarCov(lmm.red, "IT", type = "conditional") ## e_i term
getVarCov(lmm.red, "IT", type = "marginal")  # Random term
y_hat_h0 <- predict(lmm.red, data=df.re[,-3]) #lmm.red$fitted[,2]
# residuals (asymptotic estimation)
e.h0 <- df.re$net_flow -  y_hat_h0  # lmm.red$residuals[,2]
n <- length(e.h0)
B <- 1e4
T.perms <- NULL
set.seed(41703192)
for (i in seq(1,B)){
  # residuals asymptotically exchangeable
  tryCatch( {  # in case it does not converge and throws an error
    e_b <- e.h0[sample(1:n)]
    y_perm <- y_hat_h0 + e_b
    lmm.perm <- lme(fixed = y_perm ~ 1 + uni_score + uni_score:time,
                    random = ~ 1 + as.integer(time)|Country,
                    correlation= NULL, 
                    data=df.re, method = "REML")
    T.perms= c(T.perms, max(abs(summary(lmm.perm)$tTable[-c(1,2),4])))
  }, error=function(e){ print(e)})
  
}

sapply(T.perms, abs)
sum(T.perms >= T0)/length(T.perms)
layout(1)
hist(T.perms,xlim=range(c(T.perms,T0)),breaks=30,
     main="Distribution of max{t_val1, t_val2...} under H0")
abline(v=T0,col=2,lwd=2)
hist(lmm.red$residuals, breaks=30, main="Histogram of residuals")

#=====
### use a robust linear mixed effects model
# install.packages("robustlmm")
#=====
library(robustlmm)
library(lme4)
l.reg0 <- lmer(formula = net_flow ~ 1+  uni_score + 
                 (I(as.numeric(time)) | Country),
               data=df.re, REML=T)

# full model (H1)
r.reg.full <- rlmer(net_flow ~   1 + uni_score + uni_score:time +
                  (I(as.numeric(time)) | Country),
                   data=df.re,
               rho.e = smoothPsi,
               rho.b = cPsi,
               method="DASvar")
vcov(r.reg, full=T)
lme4::VarCorr(r.reg)
summ.robust.full <- summary(r.reg.full)


# t value
T0 <- max(abs(summ.robust.full$coefficients[-c(1,2),3]))

# reduced model (under H0)
r.reg.red <-  rlmer(net_flow ~ 1 + uni_score + (as.integer(time)|Country), # random effect ~N(time, Sigma)
                 rho.e = smoothPsi,
                 rho.b = cPsi,
                data=df.re, method = "DASvar" )
y_hat_h0 <- fitted.values(r.reg.red) 
# residuals (asymptotic estimation)
e.h0 <- residuals(r.reg.red)
n <- length(e.h0)
B <- 1e3
T.perms <- NULL
set.seed(41703192)
#=====
library(parallel)
library(pbapply)
detectCores()
cl=makeCluster(16)

  T.perms <- NULL
  for (i in seq(1,B)){
    # residuals asymptotically exchangeable
    tryCatch( {  # in case it does not converge and throws an error
      e_b <- e.h0[sample(1:n)]
      y_perm <- y_hat_h0 + e_b
      robustlmm.perm <- rlmer(y_perm ~ 1 + uni_score + uni_score:time+ (as.numeric(time)|Country),
                              rho.e = smoothPsi,
                              rho.b = cPsi,  data=df.re, method = "DASvar")
      T.perms= c(T.perms, max(abs(summary(robustlmm.perm)$coefficients[-c(1,2),3]) ) )
    }, error=function(e){ print(e)})
    
  }
  sapply(T.perms, abs)
 # return ( sum(T.perms >= T0)/length(T.perms) )

T.perms <- permutation_it()
clusterExport(cl,varlist=list("df.re", "permutation_it"))
Ts <- NULL
Tps=pbsapply(1, permutation_it,cl=cl)

layout(1)
hist(T.perms,xlim=range(c(T.perms,T0)),breaks=30,
     main="Distribution of max{t_val1, t_val2^...} under H0")
abline(v=T0,col=3,lwd=2)
###########
# extract covariance
#=====
var.d <- crossprod(getME(r.reg,"Lambdat"))
Zt <- getME(r.reg,"Zt")
vr <- sigma(r.reg)^2
var.b <- vr*(t(Zt) %*% var.d %*% Zt)
sI <- vr * Diagonal(nrow(df.re))
var.y <- var.b + sI
cor(ranef(r.reg))

##########
# Addendum: learnigng to use lme
#====

# l <- lm(flow_in ~ 1 + uni_score)
g <- lme(fixed = net_flow ~ 1 + uni_score + uni_score:time, random = ~ time|Country,
         correlation= NULL, #corARMA(form=~ time, p=1, q=1),
         data=df.re, method = "REML")
summary(g)
VarCorr(g)
getVarCov(g, df.re$Country,  type="random.effects")
getVarCov(g, "IT", type = "conditional")
getVarCov(g, "IT", type = "marginal")

plot(ACF(g))
getVarCov(g, individuals=df.re$Country, type = "marginal")
g$residuals
VarCorr(g)
g$contrasts
cm <- VarCorr(g)
cm <- extract.lme.cov(g, data=df.re[,-3])
shapiro.test(g$residuals)  
hist(g$residuals)

