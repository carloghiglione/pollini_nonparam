rm(list=ls())

tab1 <- read.table('academic_staff.tsv', sep='\t', header=T)
write.table(tab1[[1]], file='temp.txt', row.names = F, quote = F, col.names = F)

tab2 <- read.table('temp.txt', sep=',')
colnames(tab2) <- c('unit_measure', 'ISCED11', 'age', 'sex', 'country')

tab <- cbind(tab2, tab1[,2:8])
tab <- lapply(tab, gsub, pattern='p', replacement='')
tab <- lapply(tab, gsub, pattern='b', replacement='')
tab <- lapply(tab, gsub, pattern='e', replacement='')
tab <- lapply(tab, gsub, pattern=':', replacement=NA)
tab <- data.frame(tab)

tab <- tab[which(tab$sex=='T' & tab$age=='TOTAL'),]
tab <- cbind(tab[,1:5], sapply(sapply(tab[,6:12], as.character), as.numeric)) 


library(dplyr)
tab_t <- as_tibble(tab)

group_cols <- c("country", 'ISCED11')
num_cols <- colnames(tab)[6:12]

ris <- tab_t %>% 
  group_by(across(all_of(group_cols))) %>%
  select(num_cols) %>%
  summarise_all(sum)
  
#prendo solo gli accademici ED5-ED8
acad_staff <- ris[which(ris$ISCED11=='ED5-8'),]

#sistemo gli evidenti errori
acad_staff[which(acad_staff$country=='CZ'),] <- ris[which(ris$country=='CZ' & ris$ISCED11=='ED5'),]
acad_staff[which(acad_staff$country=='PL'),] <- ris[which(ris$country=='PL' & ris$ISCED11=='ED5'),]
acad_staff <- data.frame(acad_staff)


final_ris <- acad_staff %>% 
  select(num_cols) %>%
  transmute(acad_staff_mean = rowMeans(., na.rm = T))
final_ris <- data.frame(cbind(country = acad_staff$country, final_ris))
final_ris$acad_staff_mean[is.nan(final_ris$acad_staff_mean)]<-NA



# to visualize the groups
ris %>% group_keys()
ris %>% group_indices()
