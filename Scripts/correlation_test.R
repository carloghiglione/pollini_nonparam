correlation <- function(x, y, B = 2500, seed = 17121998){
  #H0: X ind Y vs H1
  n <- length(x)
  if(n != length(y)){
    print("Different lengths")
    return(-1)
  }
  T0 = cor(x,y)^2
  set.seed(seed)
  library(progress)
  pb <- progress_bar$new(total = B)
  pb$tick(0)
  T_stat <- numeric(B)
  for(i in 1:B){
    #under H0 we can break couples freely
    perm <- sample(1:n)
    y_p <- y[perm]
    T_stat[i] <- cor(x,y_p)^2
    pb$tick()
  }
  p_val <- sum(T_stat >= T0)/B
  return(p_val)
}
