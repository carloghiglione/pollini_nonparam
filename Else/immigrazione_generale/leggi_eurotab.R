rm(list=ls())

tab1 <- read.table('migr_imm5prv.tsv', header = T, sep = '\t')

write.table(tab1[[1]], file='temp.txt', row.names = F, col.names = F, quote = F)
tab2 <- read.table('temp.txt', sep = ',', header = F)

tab <- cbind(tab2, tab1[,-1])
colnames(tab)[1:6]<- c('from', 'age_def', 'age_class', 'unit_measure', 'sex', 'to')

tab <- tab[which(tab$age_class=='TOTAL' & tab$sex=='T' & tab$age_def=='COMPLET'), ]
rownames(tab) <- 1:dim(tab)[1]
tab <- tab[, -c(2,3,4,5)]
tab <- lapply(tab, gsub, pattern='p', replacement='')
tab <- lapply(tab, gsub, pattern='b', replacement='')
tab <- lapply(tab, gsub, pattern='e', replacement='')
tab <- lapply(tab, gsub, pattern=':', replacement=NA)
tab <- as.data.frame(tab)

to_all <- unique(tab$to)


#Visualizza plot from -> to anno per anno
from_c = 'BE'
to_c = 'FR'
sum(!is.na(tab[which(tab$from==from_c & tab$to==to_c),]))-2

