rm(list=ls())

tab_flows <- read.table('trips_EU.csv', header=T, sep=';')
tab_eco <- read.table('Eco_data_UR_11_10_21.csv', header=T, sep=';')

out_flow <- tab_flows$tot_out
in_flow <- tab_flows$tot_in
uni_score <- tab_eco$Uni_rank

cols <- c('red', 'blue')
cols_lab <- cols[tab_flows$High_flow_label+1]

x11()
par(mfrow=c(2,2))
plot(uni_score, in_flow, main = 'Inflow vs Uni Rank', pch=19, col=cols_lab)
legend('topleft', legend=c('low flow', 'high flow'), fill=cols)
plot(uni_score, out_flow, main = 'Outflow vs Uni Rank', pch=19, col=cols_lab)
legend('topleft', legend=c('low flow', 'high flow'), fill=cols)
plot(uni_score, log(in_flow), main = 'log(Inflow) vs Uni Rank', pch=19, col=cols_lab)
legend('bottomright', legend=c('low flow', 'high flow'), fill=cols)
plot(uni_score, log(out_flow), main = 'log(Outflow) vs Uni Rank', pch=19, col=cols_lab)
legend('bottomright', legend=c('low flow', 'high flow'), fill=cols)


#######################################################################
# Parametric model
Y <- in_flow
mod.p <- lm(Y ~ uni_score + I(uni_score^2))
print(summary(mod.p))

x11()
par(mfrow=c(2,2))
plot(mod.p)
print(shapiro.test(mod.p$residuals)) 

x11()
plot(uni_score, Y, main = 'Y vs Uni Rank', pch=19)
lines(uni_score[order(uni_score)], mod.p$fitted.values[order(uni_score)], col="red")


########################################################################