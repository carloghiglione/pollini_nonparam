rm(list=ls())

tab1 <- read.table('academic_staff.tsv', sep='\t', header=T)
write.table(tab1[[1]], file='temp.txt', row.names = F, quote = F, col.names = F)

tab2 <- read.table('temp.txt', sep=',')
colnames(tab2) <- c('unit_measure', 'ISCED11', 'age', 'sex', 'country')

tab <- cbind(tab2, tab1[,2:8])
tab <- lapply(tab, gsub, pattern='p', replacement='')
tab <- lapply(tab, gsub, pattern='b', replacement='')
tab <- lapply(tab, gsub, pattern='e', replacement='')
tab <- lapply(tab, gsub, pattern='d', replacement='')
tab <- lapply(tab, gsub, pattern=':', replacement=NA)
tab <- data.frame(tab)

tab <- tab[which(tab$sex=='T' & tab$age=='TOTAL' & tab$ISCED11=='ED5-8'),]

num_tab <- tab[,6:12]
num_tab <- sapply(num_tab, as.character)
num_tab <- matrix(data=sapply(num_tab, as.numeric), ncol=7, nrow=37)

acad_staff <- data.frame(cbind(country = as.character(tab$country),
                    num_staff = rowMeans(num_tab, na.rm = T)))

###############################################################################################

tab1 <- read.table('stud_aca_staff_ratio.tsv', sep='\t', header=T)
write.table(tab1[[1]], file='temp.txt', row.names = F, quote = F, col.names = F)

tab2 <- read.table('temp.txt', sep=',')
colnames(tab2) <- c('unit_measure', 'ISCED11', 'country')

tab <- cbind(tab2, tab1[,2:8])
tab <- lapply(tab, gsub, pattern='p', replacement='')
tab <- lapply(tab, gsub, pattern='b', replacement='')
tab <- lapply(tab, gsub, pattern='e', replacement='')
tab <- lapply(tab, gsub, pattern='d', replacement='')
tab <- lapply(tab, gsub, pattern=':', replacement=NA)
tab <- data.frame(tab)


tab <- tab[which(tab$ISCED11=='ED5-8'),]

num_tab <- tab[,4:10]
num_tab <- sapply(num_tab, as.character)
num_tab <- matrix(data=sapply(num_tab, as.numeric), ncol=7, nrow=36)


stud_per_staff <- cbind(country = as.character(tab$country),
                                   stud_per_staff = rowMeans(num_tab, na.rm = T))

#sistema errori
stud_per_staff <- rbind(stud_per_staff, c('CH', 300000/34775))
stud_per_staff[2,2] <- 80000/10122
stud_per_staff[18,2] <- 220000/9345.5
stud_per_staff <- data.frame(stud_per_staff)


###################################################################################
library(dplyr)
ris <- acad_staff%>%
  right_join(stud_per_staff, by='country')

colnames(ris)[1] <- 'Country'
#sistemo i nomi sbagliati
ris$Country <- as.character(ris$Country)
ris$Country[37] <- 'GB'
ris$Country[11] <- 'GR'

###################################################################################
eco_tab <- read.csv('Eco_data_UR.csv', header = T, sep=';')


eco_tab_complete <- ris %>%
  inner_join(eco_tab, by='Country')

#write.csv(eco_tab_complete, file='Eco_data_UR_12_10_21.csv', row.names = F, quote = F)