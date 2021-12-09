#We want to find a GAM for expressing the total number of people moving in and out of each country

#====
#importing the needed libraries
library(ISLR2)
library(car)
library(mgcv)
library(rgl)
library(splines)
library(pbapply)

#importing data
data <- full_tab_31_10_21
pat <- data$Patent.applications..nonresidents + data$Patent.applications..residents #we consider the total number of patents deposited in a country, ignoring if they were deposited by residents or nonresidents in said country


#First attempt: a nonparametric model using as regressor the one that our groupwise linear model was based on
g <- gam(data$tot_in + data$tot_out ~ s(data$uni_score, bs='cr'))
summary(g)

x11()
plot(g)

#however we do not have gaussianity
shapiro.test(g$residuals)
x11()
qqnorm(g$residuals)
qqline(g$residuals)

data$Country[which(abs(g$residuals)>sd(g$residuals))]

#let's try to see if it is just an outlier issue
shapiro.test(g$residuals[-which(abs(g$residuals)>sd(g$residuals))])
x11()
qqnorm(g$residuals[-which(abs(g$residuals)>sd(g$residuals))])
qqline(g$residuals[-which(abs(g$residuals)>sd(g$residuals))])
#it is not

#====
#we try to express mobility in a log scale, as was needed in the parametric approach
g1 <- gam(log(data$tot_in + data$tot_out) ~ s(data$uni_score, bs='cr'))
summary(g1)
shapiro.test(g1$residuals) #now we have Gaussianity
x11()
plot(g1)
points(data$uni_score, log(data$tot_in + data$tot_out) - g1$coefficients[1], col = 'red')

#So this is our base model, we want to see if this choice is meaningful with a variable selection approach

#====
#backward selection is not an option, we apply forward

#g_big <- gam(log(data$tot_in + data$tot_out) ~ 
#               s(data$Citations.per.document, bs='cr') +
#               s(data$num_staff, bs='cr')+
#               s(data$stud_per_staff, bs='cr')+
#               s(data$GDP, bs='cr')+
#               s(data$Reasearch, bs='cr')+
#               s(data$Education, bs='cr')+
#               s(data$LPPI, bs='cr')+
#               s(data$HDI, bs='cr')+
#               s(data$GGGI, bs='cr')+
#               s(data$PS, bs='cr')+
#               s(data$SA, bs='cr')+
#               s(I(data$Patent.applications..nonresidents + data$Patent.applications..residents), bs='cr')+
#               s(data$uni_score_norm, bs='cr')+
#               s(data$uni_score, bs='cr'))

#So we start by building a GAM model with only one nonparametric variable, for each regressor
#we keep the best one in term of Adjusted R.sq
#
#====

g_cit <- gam(log(data$tot_in + data$tot_out) ~ s(data$Citations.per.document, bs='cr'))
summary(g_cit) #0.412, significant
shapiro.test(g_cit$residuals)
x11()
plot(g_cit)
points(data$Citations.per.document, log(data$tot_in + data$tot_out) - g_cit$coefficients[1])
data$Country[which(data$Citations.per.document == max(data$Citations.per.document))] #ISLANDA


#______
g_num <- gam(log(data$tot_in + data$tot_out) ~ s(data$num_staff, bs='cr'))
summary(g_num) #0.577, significant
shapiro.test(g_num$residuals)

x11()
plot(g_num)
points(data$num_staff, log(data$tot_in + data$tot_out) - g_num$coefficients[1])


#______
g_sps <- gam(log(data$tot_in + data$tot_out) ~ s(data$stud_per_staff, bs='cr'))
summary(g_sps) #0, non significant
l_sps <- lm(log(data$tot_in + data$tot_out) ~ data$stud_per_staff)
summary(l_sps)


#______
g_GDP <- gam(log(data$tot_in + data$tot_out) ~ s(data$GDP, bs='cr'))
summary(g_GDP) #0.451, significant
shapiro.test(g_GDP$residuals)

x11()
plot(g_GDP)
points(data$GDP, log(data$tot_in + data$tot_out) - g_GDP$coefficients[1])


#______
g_res <- gam(log(data$tot_in + data$tot_out) ~ s(data$Reasearch, bs='cr'))
summary(g_res) #0.372, significant
shapiro.test(g_res$residuals)

x11()
plot(g_res)
points(data$Reasearch, log(data$tot_in + data$tot_out) - g_res$coefficients[1])


#______
g_edu <- gam(log(data$tot_in + data$tot_out) ~ s(data$Education, bs='cr'))
summary(g_edu) #0.04, non significant


#______
g_LPPI <- gam(log(data$tot_in + data$tot_out) ~ s(data$LPPI, bs='cr'))
summary(g_LPPI) #0.199, significant
shapiro.test(g_LPPI$residuals)

x11()
plot(g_LPPI)
points(data$LPPI, log(data$tot_in + data$tot_out) - g_LPPI$coefficients[1])


#______
g_hdi <- gam(log(data$tot_in + data$tot_out) ~ s(data$HDI, bs='cr'))
summary(g_hdi) #0.226, significant
shapiro.test(g_hdi$residuals)

x11()
plot(g_hdi)
points(data$HDI, log(data$tot_in + data$tot_out) - g_hdi$coefficients[1])


#______
g_gg <- gam(log(data$tot_in + data$tot_out) ~ s(data$GGGI, bs='cr'))
summary(g_gg) #0.14, non significant


#______
g_ps <- gam(log(data$tot_in + data$tot_out) ~ s(data$PS, bs='cr'))
summary(g_ps) #0.471, significant
shapiro.test(g_ps$residuals)

x11()
plot(g_ps)
points(data$PS, log(data$tot_in + data$tot_out) - g_ps$coefficients[1])


#______
g_sa <- gam(log(data$tot_in + data$tot_out) ~ s(data$SA, bs='cr'))
summary(g_sa) #0.04, non significant


#______
g_pat <- gam(log(data$tot_in + data$tot_out) ~ s(pat, bs='cr'))
summary(g_pat) #0.545, significant
shapiro.test(g_pat$residuals)

x11()
plot(g_pat)
points(data$Patent.applications..nonresidents+data$Patent.applications..residents, log(data$tot_in + data$tot_out) - g_pat$coefficients[1])


#______
g_usnorm <- gam(log(data$tot_in + data$tot_out) ~ s(data$uni_score_norm, bs='cr'))
summary(g_usnorm) #0.598, significant
shapiro.test(g_usnorm$residuals)

x11()
plot(g_usnorm)
points(data$uni_score_norm, log(data$tot_in + data$tot_out) - g_usnorm$coefficients[1])


#====
g_us <- gam(log(data$tot_in/data$tot_out) ~ s(data$uni_score, bs='cr'))
summary(g_us) #0.869, significant
shapiro.test(g_us$residuals)

x11()
plot(g_us)
points(data$uni_score, log(data$tot_in) - g_us$coefficients[1])

#university score is by far the best
#we try to see if we can get a better model by considering two regressors

#====
gbig_ns <-  gam(log(data$tot_in + data$tot_out) ~ 
                  s(data$num_staff, bs='cr')+
                  s(data$uni_score, bs='cr'))
summary(gbig_ns) #no significant improvement


#______
gbig_ss <-  gam(log(data$tot_in + data$tot_out) ~ 
                  s(data$stud_per_staff, bs='cr')+
                  s(data$uni_score, bs='cr'))
summary(gbig_ss) #no significant improvement


#______
gbig_GDP <-  gam(log(data$tot_in + data$tot_out) ~ 
                   s(data$GDP, bs='cr')+
                   s(data$uni_score, bs='cr'))
summary(gbig_GDP) #no significant improvement


#______
gbig_res <-  gam(log(data$tot_in + data$tot_out) ~ 
                   s(data$Reasearch, bs='cr')+
                   s(data$uni_score, bs='cr'))
summary(gbig_res) #there may be a significant improvement
shapiro.test(gbig_res$residuals)

x11()
plot(gbig_res) #however interpretability of the effect of R&D expense leaves a lot to be desired


#______
gbig_edu <-  gam(log(data$tot_in + data$tot_out) ~ 
                   s(data$Education, bs='cr')+
                   s(data$uni_score, bs='cr'))
summary(gbig_edu) #no significant improvement


#______
gbig_LPPI <-  gam(log(data$tot_in + data$tot_out) ~ 
                    s(data$LPPI, bs='cr')+
                    s(data$uni_score, bs='cr'))
summary(gbig_LPPI) #no significant improvement


#______
gbig_HDI <-  gam(log(data$tot_in + data$tot_out) ~ 
                   s(data$HDI, bs='cr')+
                   s(data$uni_score, bs='cr'))
summary(gbig_HDI) #no significant improvement


#______
gbig_gg <-  gam(log(data$tot_in + data$tot_out) ~ 
                  s(data$GGGI, bs='cr')+
                  s(data$uni_score, bs='cr'))
summary(gbig_gg) #no significant improvement


#______
gbig_ps <-  gam(log(data$tot_in + data$tot_out) ~ 
                  s(data$PS, bs='cr')+
                  s(data$uni_score, bs='cr'))
summary(gbig_ps) #no significant improvement


#______
gbig_sa <-  gam(log(data$tot_in + data$tot_out) ~ 
                  s(data$SA, bs='cr')+
                  s(data$uni_score, bs='cr'))
summary(gbig_sa) #no significant improvement


#______
gbig_pat <- gam(log(data$tot_in + data$tot_out) ~ 
                  s(pat, bs='cr')+
                  s(data$uni_score, bs='cr'))
summary(gbig_pat) #there may be a significant improvement
shapiro.test(gbig_pat$residuals)

x11()
plot(gbig_pat) #indeed we do not have a strong evidence toward the significance of the number of patents

#we then consider an appropriate test
gbase <-  gam(log(data$tot_in + data$tot_out) ~s(data$uni_score, bs='cr'))
#Test to see if Pat is significant
anova(gbase, gbig_pat, test = 'F')
#we cannot accept at level 5%


#====
#______
gbig_usn <-  gam(log(data$tot_in + data$tot_out) ~ 
                   s(data$uni_score_norm, bs='cr')+
                   s(data$uni_score, bs='cr'))
summary(gbig_usn) #no significant improvement


#______
gbig_nr <-  gam(log(data$tot_in + data$tot_out) ~ 
                  s(data$num_ric, bs='cr')+
                  s(data$uni_score, bs='cr'))
summary(gbig_nr) #no significant improvement
#====

#___________________________________________________________________________________
#so we confirm our initial choice for a model
gmod <- gam(log(data$tot_in + data$tot_out) ~ s(data$uni_score, bs='cr'))
summary(gmod)

x11()
plot(gmod)

#___________________________________________________________________________________

#We then try to do the same for having an interpretaion of our regressor

#usmod <- gam(data$uni_score ~ 
#               s(data$Citations.per.document, bs = 'cr')+
#               s(data$num_staff, bs = 'cr')+
#               s(data$stud_per_staff, bs = 'cr')+
#               s(data$Reasearch, bs = 'cr')+
#               s(data$Education, bs = 'cr')+
#               s(pat, bs = 'cr')+
#               s(data$num_ric, bs = 'cr'))

#====
#again forward selection, too many variables and too little data to apply backward selection
usmod_cd <- gam(data$uni_score ~ s(data$Citations.per.document ,bs = 'cr'))
summary(usmod_cd)#0.192, non sign

usmod_ns <- gam(data$uni_score ~ s(data$num_staff ,bs = 'cr'))
summary(usmod_ns)#0.319, sign
shapiro.test(usmod_ns$residuals) #non gaussian
x11()
plot(usmod_ns)
points(data$num_staff, data$uni_score - usmod_ns$coefficients[1])

usmod_ss <- gam(data$uni_score ~ s(data$stud_per_staff ,bs = 'cr'))
summary(usmod_ss)#0
l_ss <- lm(data$uni_score~data$stud_per_staff)
summary(l_ss)

usmod_re <- gam(data$uni_score ~ s(data$Reasearch ,bs = 'cr'))
summary(usmod_re)#0.08 non sig

usmod_ed <- gam(data$uni_score ~ s(data$Education ,bs = 'cr'))
summary(usmod_ed)#0.05 non sig

usmod_pat <- gam(data$uni_score ~ s(pat ,bs = 'cr'))
summary(usmod_pat)#0.563 sig
shapiro.test(usmod_pat$residuals)#nongaussian
x11()
plot(usmod_pat)
points(pat, data$uni_score - usmod_pat$coefficients[1])

usmod_nr <- gam(data$uni_score ~ s(data$num_ric ,bs = 'cr'))
summary(usmod_nr)#0.691 sig
shapiro.test(usmod_nr$residuals)#nongaussian


usmod_big <- gam(data$uni_score ~
                   s(data$num_ric, bs = 'cr')+
                   s(data$Reasearch, bs = 'cr')+
                   s(pat, bs = 'cr'))
summary(usmod_big)
shapiro.test(usmod_big$residuals)
x11()
plot(usmod_big)

usmod_big2 <- gam(data$uni_score ~
                    s(data$num_ric, bs = 'cr')+
                    s(pat, bs = 'cr'))
summary(usmod_big2)
shapiro.test(usmod_big2$residuals)
x11()
plot(usmod_big2)
#these models do not satisfy us

x11()
plot(data$uni_score, log(pat))
#there may be a logarithmic behaviour

uslog_pat <- gam(data$uni_score ~ s(log(pat), bs = 'cr'))
summary(uslog_pat)#0.478
x11()
plot(uslog_pat)

uslog_ns <- gam(data$uni_score ~ s(log(data$num_staff), bs = 'cr'))
summary(uslog_ns)#0.299

uslog_ss <- gam(data$uni_score ~ s(log(data$stud_per_staff), bs = 'cr'))
summary(uslog_ss)#0

uslog_res <- gam(data$uni_score ~ s(log(data$Reasearch), bs = 'cr'))
summary(uslog_res)#0.08

uslog_edu <- gam(data$uni_score ~ s(log(data$Education), bs = 'cr'))
summary(uslog_edu)#0.04 non sig

uslog_nr <- gam(data$uni_score ~ s(log(data$num_ric), bs = 'cr'))
summary(uslog_nr)#0.503

#====

uslog_big <- gam(data$uni_score ~
                   s(data$num_ric, bs = 'cr')+
                   s(data$Reasearch, bs = 'cr'))
summary(uslog_big)
shapiro.test(uslog_big$residuals)
x11()
plot(uslog_big)
#these models are not satisfying, we stick with the linear regression approach

