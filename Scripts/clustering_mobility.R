#CLUSTERING

#We are looking for a clustering structure amongst our countries based on the difference in mobility
#Each country is represented in R2 by (total number of people going in, total number of people going out) during our time horizon
#we first of all build the distance matrix (euclidian distance)
diss <- dist(Trips_ext$tot_in, Trips_ext$tot_out, method = 'euclidian')

#Ward is the only clusterubg structure that does not just separate UK from the rest of Europe
clus <- hclust(diss, method = 'ward.D2')
x11()
plot(clus)

lab <- cutree(clus, 2)
flowlab <- cbind(Country = Trips_ext$Country, label = lab)
flowlab 
Trips_ext$Country[which(lab == "2")]

x11()
plot(Trips_ext$tot_in, Trips_ext$tot_out, main = "Clustering structure", xlab = "Ingoing flow", ylab = "Outgoing flow")
points(Trips_ext$tot_in[which(lab == "1")], Trips_ext$tot_out[which(lab == "1")], col = 'green')
points(Trips_ext$tot_in[which(lab == "2")], Trips_ext$tot_out[which(lab == "2")], col = 'red')
legend(fill = c('green', 'red'),legend = c("First cluster", "Second cluster"), "topleft")


x11()
plot(log(Trips_ext$tot_in), log(Trips_ext$tot_out), main = "Clustering structure", xlab = "Ingoing flow", ylab = "Outgoing flow")
points(log(Trips_ext$tot_in[which(lab == "1")]), log(Trips_ext$tot_out[which(lab == "1")]), col = 'green')
points(log(Trips_ext$tot_in[which(lab == "2")]), log(Trips_ext$tot_out[which(lab == "2")]), col = 'red')
legend(fill = c('green', 'red'),legend = c("First cluster", "Second cluster"), "topleft")


#So we can interpret the resulting clustering structure as "High mobility countries" vs "Low mobility countries"
