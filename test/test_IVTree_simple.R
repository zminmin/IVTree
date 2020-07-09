# This program tests out the features of the IVTree package


library(IVTree)


IVTree_mod <- function(form, data_tr, data_es, data_te, if_prune){
  tree <- honest.IVTree(form, data = data_tr, treatment = data_tr$T1, treatment1 = data_tr$T, IV = data_tr$IV,
                            split.Honest = T, cv.option = T, cv.Honest = T, 
                            minsize = 25, split.alpha = 0.5, cv.alpha = 0.5, split.Bucket = F, 
                            est_data = data_es, est_treatment = data_es$T1, est_treatment1 = data_es$T, est_IV = data_es$IV)

  if (if_prune) {
    opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
    tree_prune <- prune(tree, opcp) #pruning
  }
  else{
    tree_prune <- tree #no pruning
  }

  return(tree_prune)
}


load("data/data_all.RData")
load("data/data_te.RData")

# ============ if the original data is not centered, center ===================
# there are four X, one Y, one instrumental variable IV and one W named as T1
# data_all = data.centering(data_all[,1:4], data_all$Y, data_all$T1, data_all$IV) 
# ============ ******************************************** ===================



n <- nrow(data_all)
train_ind <- sample(1:(n), size = n/2)
data_tr = data_all[train_ind,]
data_es = data_all[-train_ind,]

form <- paste('Y ~', paste(paste0('X', 1:(4)), collapse = ' + '))

for(i in 1:200){
  ivtree = IVTree_mod(form, data_tr, data_es, data_te, if_prune = TRUE)
}


# rpart.plot(ivtree)