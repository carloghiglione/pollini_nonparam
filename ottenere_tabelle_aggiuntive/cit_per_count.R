rm(list=ls())

tab <- read.csv('citations_per_country.csv', sep=';')

iso2 <- read.csv('iso2codes.csv', sep = ',')
iso2 <- iso2[,c(1,2)]
iso2$Alpha.2.code <- gsub( " ", "", iso2$Alpha.2.code) 

#attacco i codici iso
library(dplyr)
tab_new <- tab %>%
  inner_join(iso2, by='Country')

tab_new <- tab_new[, c(10, 8)]
colnames(tab_new)[1] <- 'Country'


Eco_data_UR <- read.csv('Eco_data_UR_12_10_21.csv', header = T)

Eco_data_UR_new <- tab_new %>%
  inner_join(Eco_data_UR, by='Country')

write.csv(Eco_data_UR_new, file='Eco_data_UR_12_10_21.csv', row.names = F, quote = F)
