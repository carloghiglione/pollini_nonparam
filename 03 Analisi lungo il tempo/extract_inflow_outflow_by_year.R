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


####################################################################
# NORMALIZED FLOWS

full_tab <- read.csv('full_tab_31_10_21.csv')


ratio_elements <- colnames(inflow[,-1])
num_ric <- data.frame(country=full_tab$Country, num_ric=full_tab$num_ric)

# normalize inflow
norm_inflow <- num_ric %>%
  inner_join(inflow, by='country')

norm_inflow <- norm_inflow[,-1] %>%
  mutate(across(everything()), . / num_ric) %>%          # divide each column by num_ric
  select(ratio_elements) %>%                             # hold columns in ratio_elements
  mutate(country = inflow$country, .before = '2006')     # add column for country


# normalize outflow
norm_outflow <- num_ric %>%
  inner_join(outflow, by='country')

norm_outflow <- norm_outflow[,-1] %>%
  mutate(across(everything()), . / num_ric) %>%          # divide each column by num_ric
  select(ratio_elements) %>%                             # hold columns in ratio_elements
  mutate(country = outflow$country, .before = '2006')    # add column for country

write.csv(norm_inflow, file = 'norm_year_inflow.csv', row.names = F)
write.csv(norm_outflow, file = 'norm_year_outflow.csv', row.names = F)

