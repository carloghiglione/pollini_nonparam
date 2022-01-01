rm(list=ls())

tab <- read.csv('QS_ranking2018.csv', header=T)
rank_tab <- data.frame(rank=tab$X2018, Country=tab$Country, uni=tab$Institution.Name )

iso2 <- read.csv('iso2codes.csv', header=T)
iso2 <- iso2[,c(1,2)]

library(dplyr)
rank_tab <- rank_tab %>%
  inner_join(iso2, by='Country')

rank_tab <- data.frame(rank = rank_tab$rank, Country=rank_tab$Alpha.2.code, uni=rank_tab$uni)
rank_tab$rank <- as.numeric(row.names(rank_tab))
rank_tab$Country <- gsub( " ", "", rank_tab$Country)

country_list <- read.csv('Eco_data_EXT_12_10_21.csv', header=T)
country_list <- data.frame(Country = country_list$Country)

rank_tab <- rank_tab %>%
  semi_join(country_list, by='Country')

by_fact <- 100
separ <- seq(0,1000, by=by_fact)
rank_tab$class <- numeric(379)

for(i in 1:(1000/by_fact)){
  rank_tab$class[which(rank_tab$rank > separ[i] & rank_tab$rank <= separ[i+1])] <- i
}

###############################################
#group the data

group_cols <- c("Country", 'class')

group_rank <- rank_tab[,c(2,4)] %>%
  group_by(across(all_of(group_cols))) %>%
  count()

##############################################
# table of number of universities for each country
#fonte uniRank

num_uni <- data.frame(Country = country_list$Country,
                      num_uni = c(
                        144, 90, 80, 92, 82, 37, 60, 32, 19, 24, 59, 35, 33, 35, 25, 29,
                        70, 31, 19, 30, 20, 41, 36, 9, 17, 10, 1, 5, 20, 6, 20, 190, 410, 32, 161
                      ))


##############################################
#compute the university scores

uni_score <- numeric(35)
uni_score_norm_top <- numeric(35)
uni_score_norm_all <- numeric(35)
uni_inv <- numeric(35)
uni_inv_norm <- numeric(35)
for(i in 1:35){
  curr_country <- country_list$Country[i]
  if(curr_country %in% unique(group_rank$Country)){
    curr_scores <- group_rank[which(group_rank$Country == curr_country),]
    uni_score_norm_top[i] <- sum(((1000/by_fact+1-curr_scores$class))*curr_scores$n)/sum(curr_scores$n)
    uni_score[i] <- sum((1000/by_fact+1-curr_scores$class)*curr_scores$n)
    uni_score_norm_all[i] <- sum((1000/by_fact+1-curr_scores$class)*curr_scores$n)/num_uni$num_uni[i]
    uni_inv_norm[i] <- sum((1/curr_scores$class)*curr_scores$n)/num_uni$num_uni[i]
    uni_inv[i]<- sum((1/curr_scores$class)*curr_scores$n)
  } else{
    uni_score[i] <- 0
  }
}

tab_uni_score <- data.frame(Country = country_list$Country,
                            uni_score_norm = uni_inv_norm,
                            uni_score = uni_inv)
################################################################################
# Aggregate to tabel
old_tab <- read.csv('ECO_data_EXT_12_10_21.csv', header = T)

new_tab <- old_tab %>%
  inner_join(tab_uni_score, by='Country')

new_tab <- new_tab[,-15]

write.csv(new_tab, file='Eco_data_EXT_13_10_21.csv', row.names = F, quote = F)
