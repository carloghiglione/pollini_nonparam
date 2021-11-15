# perform K-Fold crossvalidation composing the folds with an equal number of data
# coming from the subsets of a given regressor defined by its evenly spaced quantiles

quantile_kfold <- function(n_fold, df, col_name){
  require(dplyr)
  
  n_data <- dim(df)[1]                       # number of observations in the dataframe
  n_data_per_base <- floor(n_data/n_fold)    # number of base sets to build the folds
  
  delimeters <- quantile(scores.pc1, seq(0,1, length.out = n_data_per_base+1))         # evenly spaced quantiles to delimit the base sets ranges
  df$groups <- cut(sapply(df[col_name], as.numeric), delimeters, include.lowest = T)   # factor assigning the base set to each observation
  
  # build the base sets
  base_sets <- df %>%
    group_split(groups)
  
  folds <- vector('list', n_fold)     # initialize the folds
  for(i in 1:length(base_sets)){
    n_curr <- as.numeric(base_sets[[i]] %>% summarise(n=n()))    # current number of data in the base set
    if(n_curr == n_fold){
      fold_assign <- sample(n_fold, replace = F)                         # assign the fold number for each element of the base set
    }else{
      fold_assign <- c(sample(n_fold, replace = F), sample(n_fold, 1))   # assign the fold number for each element of the base set
    }
    for(j in 1:n_curr){
      folds[[fold_assign[j]]] <- rbind(folds[[fold_assign[j]]], base_for_folds[[i]] %>%                    # add the current element in the base set to the corresponding fold
                                                                                    slice(j) %>% 
                                                                                    select(-groups) %>% 
                                                                                    as.data.frame)
    }
  }
  
  # aggregate the folds in the n_fold datasets
  full_folds <- vector('list', n_fold)
  for(i in 1:n_fold){
    full_folds[[i]] <- bind_rows(folds[seq(1,n_fold)[-i]])
  }
  return(full_folds)
}