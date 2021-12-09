rm(list=ls())

library(dplyr)
library(tidyr)


# read data
tab <- read.csv('Flows_per_year.csv')


# list of countries
countries <- read.csv('country_list.csv')



# extract inflow by country and by year
inflow <- tab %>%
  group_by(across(all_of(c('endC', 'year')))) %>%
  summarise(inflow=sum(count)) %>%
  filter(endC %in% countries$Country)

# extract outflow by country and by year
outflow <- tab %>%
  group_by(across(all_of(c('startC', 'year')))) %>%
  summarise(outflow=sum(count)) %>%
  filter(startC %in% countries$Country)


# reshape data with the countries as rows
inflow <- inflow %>%
  group_by(endC) %>%
  pivot_wider(names_from = year, values_from = inflow) %>%
  rename(country = endC)


# reshape data with the countries as rows
outflow <- outflow %>%
  group_by(startC) %>%
  pivot_wider(names_from = year, values_from = outflow) %>%
  rename(country = startC)


write.csv(inflow, file = 'year_inflow.csv', row.names = F)
write.csv(outflow, file = 'year_outflow.csv', row.names = F)
