rm(list=ls())
full_tab <-  read.csv('migration_data.csv', header=T)

#select only data of inflow of foreing population by country

tab <- full_tab[which(full_tab$VAR=='B11'), ]
tab <- tab[, c(1, 2, 7, 8, 10, 11)]
colnames(tab)[c(1,2,3,4)] <- c('from', 'from_long', 'to', 'to_long' )

#Visualizza plot from -> to anno per anno
from_c = 'ITA'
to_c = 'FRA'
plot(tab[which(tab$from==from_c & tab$to==to_c),]$Year,
     tab[which(tab$from==from_c & tab$to==to_c),]$Value, 
     type='b', xlab='year', ylab='number', main=paste('from', from_c, 'to', to_c))


from_c <- unique(tab$from)
to_c <- unique(tab$to)
#avail_mat <- matrix(0, length(from_c), length(to_c))
avail_mat <- matrix(0, length(to_c), length(to_c))
colnames(avail_mat) <- to_c
#rownames(avail_mat) <- from_c
rownames(avail_mat) <- to_c


#trovo tabella di quanti anni di dati di ogni country to ho per ogni country from
#come from considero solo i dti dei paesi to che sono quelli importanti
for(i in to_c){
   for(j in to_c){
     avail_mat[i, j] <- length(which(tab$from==i & tab$to==j))
   }
}

#vedo che ci sono nazioni grosse con molti dati mancanti