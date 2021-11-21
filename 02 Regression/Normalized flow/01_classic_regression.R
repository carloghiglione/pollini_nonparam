rm(list=ls())

full_tab <- read.csv('full_tab_31_10_21.csv')

##################
# creo tabelle con variabili di interesse
full_tab <- full_tab[-4,]

flow <- full_tab$tot_out - full_tab$tot_in    # out - in
flow.norm <- flow / full_tab$num_ric

countries <- full_tab$Country
Z <- data.frame(flow, full_tab[,-c(1,2,3,5,19)])
Z.norm <- data.frame(flow.norm, full_tab[,-c(1,2,3,5,19)])


##################
# visulize data
x11()
pairs(Z)

x11()
pairs(Z.norm)


##################
# build linear model
mod <- lm(flow ~., data = Z)
mod.norm <- lm(flow.norm ~., data = Z.norm)

summary(mod)
summary(mod.norm)

x11()
par(mfrow=c(2,2))
plot(mod)
shapiro.test(mod$residuals)

x11()
par(mfrow=c(2,2))
plot(mod.norm)
shapiro.test(mod.norm$residuals)


###################
# collinearity analysis
library(car)

vif(mod)

Z <- Z[, -c(7,8,10,13)]
Z.norm <- Z.norm[, -c(7,8,10,13)]

mod <- lm(flow ~., data = Z)
mod.norm <- lm(flow.norm ~., data = Z.norm)

summary(mod)
summary(mod.norm)

x11()
par(mfrow=c(2,2))
plot(mod)
shapiro.test(mod$residuals)

x11()
par(mfrow=c(2,2))
plot(mod.norm)
shapiro.test(mod.norm$residuals)

x11()
pairs(Z)

x11()
pairs(Z.norm)


#############################
# Conclusioni
# la non normalizzata mi sembra meglio