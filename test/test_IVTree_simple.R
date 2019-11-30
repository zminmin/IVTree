# This program tests out the features of the IVTree package


library(IVTree)


IVTree_mod <- function(form, data_tr, data_es, data_te){
  tree <- honest.IVTree(form, data = data_tr, treatment = data_tr$T1, treatment1 = data_tr$T, IV = data_tr$IV,
                            split.Honest = T, cv.option = T, cv.Honest = T, 
                            minsize = 25, split.alpha = 0.5, cv.alpha = 0.5, split.Bucket = F, 
                            est_data = data_es, est_treatment = data_es$T1, est_treatment1 = data_es$T, est_IV = data_es$IV)

  opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
#  tree_prune <- prune(tree, opcp)
  tree_prune <- tree #no pruning
  leaf=sum(tree_prune$frame$var=='<leaf>')  
  pd=predict(tree_prune, newdata=data_te, type="vector")   
  
  pd_es=predict(tree_prune, newdata=data_es, type="vector")
  mincov=Inf
  maxcov=-Inf
  for (val in unique(pd_es)){
    temp=data_es[pd_es==val,]
    if(cov(temp$T,temp$IV)<mincov)
      mincov=cov(temp$T,temp$IV)
    if(cov(temp$T,temp$IV)>maxcov)
      maxcov=cov(temp$T,temp$IV)
  }
  split_freq = matrix(0, 30, ncol(data_tr)-6)
  nodes <- as.numeric(rownames(tree_prune$frame))
  for (index in 1:nrow(tree_prune$frame)){
   for (variable in 1:ncol(split_freq)){
     for (depth in 1:min(max(rpart:::tree.depth(nodes)), 30)){
	   if (rpart:::tree.depth(nodes)[index] == depth-1 & tree_prune$frame$var[index] == paste0("X", variable)){
	   split_freq[depth, variable] = split_freq[depth, variable] + 1 
	   }
	 }  
   }
  }
  list(nodes=leaf,prediction=pd,tree=tree_prune,mincov=mincov,maxcov=maxcov, split_freq=split_freq)
}



load("data/data_all.RData")
load("data/data_te.RData")

# ============ if the original data is not centered, center ===================
# there are five X, one Y, one instrumental variable IV and one W named as T1
# data_all = data.centering(data_all[,1:4], data_all$Y, data_all$T1, data_all$IV) 
# ============ ******************************************** ===================


report = data.frame()
report.all=data.frame()
prediction = data.frame()
freq_table = matrix(0, 30, ncol(data_all)-6)

for (i in 1:100){
  n <- nrow(data_all)
	train_ind <- sample(1:(n), size = n/2)
	data_tr = data_all[train_ind,]
	data_es = data_all[-train_ind,]

  form <- paste('Y ~', paste(paste0('X', 1:(4)), collapse = ' + '))

	temp = IVTree_mod(form, data_tr, data_es, data_te)
	mse = mean((temp$prediction-data_te$kappa)^2)
	report = rbind(report, c(temp$nodes, mse, temp$mincov, temp$maxcov))

	prediction = rbind(prediction, temp[[2]])
	freq_table = freq_table + temp$split_freq
}

names(report) = c('no_leaf', 'mse', 'mincov', 'maxcov')
temp2 = apply(prediction, 2, function(x)quantile(x,c(0.025,0.975)))
report$cover_prob = mean(data_te$kappa<=temp2[2,] & data_te$kappa>=temp2[1,])
temp3 = apply(prediction, 2, function(x)mean(x))
report$mse_all = mean((temp3-data_te$kappa)^2)
# report$mse_all
report.all = rbind(report.all, colMeans(report))
names(report.all)=c('no_leaf', 'mse', 'mincov', 'maxcov', 'cover_prob', 'mse')
print(report.all)