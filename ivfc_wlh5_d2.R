##presetting.T
library(causalTree)
library(LaplacesDemon)
library(dplyr)

library(grf)
centering <- function(X, Y, W, Z){

forest.Y <- regression_forest(X, Y)
Y.hat = predict(forest.Y)$predictions
        
forest.W <- regression_forest(X, W)
W.hat = predict(forest.W)$predictions

forest.Z <- regression_forest(X, Z)
Z.hat = predict(forest.Z)$predictions
        
data <- cbind(X, Y = Y - Y.hat, T = W - W.hat, IV = Z - Z.hat)

data 
}

data_build_0 <- function(n){
  dt <- rbinom(n * 3, 1, 0.5)
  dim(dt) <- c(n, 3)
  dt <- data.frame(dt)
  names(dt) <- c(paste0('X', 1:3))
  dt$eta <- dt$X1 + dt$X2
  dt$kappa <- dt$X2 
  dt$phi <- rnorm(n, mean=0, sd=sqrt(1))
  dt$IV <- rnorm(n)
  dt$T <- cp1*dt$IV + cp2*dt$phi + sqrt(2/pi-cp1^2-cp2^2)*rnorm(n)
  dt$T <- (sign(dt$T)+1)/2
  dt$Y <- dt$eta + dt$T * dt$kappa + dt$phi + rnorm(n, sd = 0.1)
  dt1 <- centering(dt[,1:3], dt$Y, dt$T, dt$IV)
  dt1$kappa = dt$kappa
  dt1$T1 = dt$T
  dt1
}
data_build_1 <- function(n){
  dt <- rbinom(n * 5, 1, 0.5)
  dim(dt) <- c(n, 5)
  dt <- data.frame(dt)
  names(dt) <- c(paste0('X', 1:5))
  dt$eta <- dt$X1 + dt$X2 + dt$X3
  dt$kappa <- dt$X3
  dt$phi <- rnorm(n,mean=0,sd=sqrt(1)) #dt$X5
  dt$IV <- rnorm(n)
  dt$T <- cp1*dt$IV + cp2*dt$phi + sqrt(2/pi-cp1^2-cp2^2)*rnorm(n)
  dt$T <- (sign(dt$T)+1)/2
  dt$Y <- dt$eta + dt$T * dt$kappa + dt$phi + rnorm(n, sd = 0.1)
  dt1 <- centering(dt[,1:5], dt$Y, dt$T, dt$IV)
  dt1$kappa = dt$kappa
  dt1$T1 = dt$T
  dt1
}
data_build_2 <- function(n){
  dt <- rbinom(n * 10, 1, 0.5)
  dim(dt) <- c(n, 10)
  dt <- data.frame(dt)
  names(dt) <- c(paste0('X', 1:10))
  dt$eta <-  dt$X1 + dt$X2 + dt$X3 + dt$X4 + dt$X5 + dt$X6
  dt$kappa <- dt$X5 + dt$X6 
  dt$phi <- rnorm(n,mean=0,sd=sqrt(1)) #dt$X10
  dt$IV <- rnorm(n)
  dt$T <- cp1*dt$IV + cp2*dt$phi + sqrt(2/pi-cp1^2-cp2^2)*rnorm(n)
  dt$T <- (sign(dt$T)+1)/2
  dt$Y <- dt$eta + dt$T * dt$kappa + dt$phi + rnorm(n, sd = 0.1) #runif(n, min=0, max=0.1) # #rexp(n, rate = 10)
  dt1 <- centering(dt[,1:10], dt$Y, dt$T, dt$IV)
  dt1$kappa = dt$kappa
  dt1$T1 = dt$T
  dt1
}
data_build_3 <- function(n){
  dt <- rbinom(n * 20, 1, 0.5)
  dim(dt) <- c(n, 20)
  dt <- data.frame(dt)
  names(dt) <- c(paste0('X', 1:20))
  dt$eta <- dt$X1 + dt$X2 + dt$X3 + dt$X4 + dt$X5 + dt$X6 + dt$X7 + dt$X8 + dt$X9 + dt$X10 + dt$X11 + dt$X12
  dt$kappa <- dt$X9 + dt$X10 + dt$X11 + dt$X12
  dt$phi <-rnorm(n, mean=0, sd=sqrt(1)) #dt$X20
  dt$IV <- rnorm(n)
  dt$T <- cp1*dt$IV + cp2*dt$phi + sqrt(2/pi-cp1^2-cp2^2)*rnorm(n)
  dt$T <- (sign(dt$T)+1)/2
  dt$Y <- dt$eta + dt$T * dt$kappa + dt$phi + rnorm(n, sd = 0.1)
  dt1 <- centering(dt[,1:20], dt$Y, dt$T, dt$IV)
  dt1$kappa = dt$kappa
  dt1$T1 = dt$T
  dt1
}
data_build_4 <- function(n){
  dt <- rbinom(n * 30, 1, 0.5)
  dim(dt) <- c(n, 30)
  dt <- data.frame(dt)
  names(dt) <- c(paste0('X', 1:30))
  dt$eta <- dt$X1 + dt$X2 + dt$X3 + dt$X4 + dt$X5 + dt$X6 + dt$X7 + dt$X8 + dt$X9 + dt$X10 + dt$X11 + dt$X12 + dt$X13 + dt$X14 + dt$X15 + dt$X16
          + dt$X17 + dt$X18 
  dt$kappa <- dt$X13 + dt$X14 + dt$X15 + dt$X16 + dt$X17 + dt$X18 
  dt$phi <- rnorm(n, mean=0, sd=sqrt(1)) #dt$X30
  dt$IV <- rnorm(n)
  dt$T <- cp1*dt$IV + cp2*dt$phi + sqrt(2/pi-cp1^2-cp2^2)*rnorm(n)
  dt$T <- (sign(dt$T)+1)/2
  dt$Y <- dt$eta + dt$T * dt$kappa + dt$phi + rnorm(n, sd = 0.1)
  dt1 <- centering(dt[,1:30], dt$Y, dt$T, dt$IV)
  dt1$kappa = dt$kappa
  dt1$T1 = dt$T
  dt1
}

IV_causalTree_mod=function(form,data_tr,data_es,data_te){
  tree <- honest.causalTree(form, data=data_tr, treatment=data_tr$T1, treatment1=data_tr$T, IV=data_tr$IV,
                            split.Rule='CT', split.Honest=T, cv.option='CT', cv.Honest=T, 
                            minsize = 25, split.alpha = 0.5, cv.alpha = 0.5, split.Bucket = F, 
                            est_data=data_es, est_treatment=data_es$T1, est_treatment1=data_es$T, est_IV=data_es$IV)
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
# SA_causalTree_mod=function(form, data_gen, n){
#   tree <- honest.causalTree(form, data=data_tr, treatment=data_tr$T,
#                             split.Rule='CT', split.Honest=T, cv.option='CT', cv.Honest=T, 
#                             minsize = 25, split.alpha = .5, cv.alpha = .5,
#                             est_data=data_es, est_treatment=data_es$T)
#   opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
#   tree_prune <- prune(tree, opcp)
  
#   leaf=sum(tree_prune$frame$var=='<leaf>')
#   pd=predict(tree_prune, newdata=data_te, type="vector")
#   mse=mean((pd - data_te$kappa)^2)
#   tree2 <- honest.causalTree(form, data=data_tr, treatment=data_tr$T,
#                              split.Rule='CT', split.Honest=T,cv.option='CT', cv.Honest=T,
#                              minsize = 25, split.alpha = .5, cv.alpha = .5,
#                              est_data=data_te, est_treatment=data_te$T)
#   opcp <- tree2$cptable[,1][which.min(tree2$cptable[,4])]
#   tree2_prune <- prune(tree2, opcp)

#   mse2=mean((predict(tree2_prune, newdata=data_te, type="vector") - pd)^2)
#   list(leaf,pd,mse,mse2)
# }


##smlt_ivf.R

cp1=0.6
omitted=T

report.all=data.frame()
for (s in 1:10){
for (cp2 in c(0.5)){
for (k in c(5)){#
for (n in c(1000, 2000, 3000, 4000, 5000)){#

  start_time <- proc.time()


  set.seed(s) 
  print(paste("CONFIG", s, cp1, cp2, k, n, sep="-"))
  if(k==5){
  form=paste('Y ~', paste(paste0('X', 1:(5-omitted)), collapse = ' + '))
  data_gen=data_build_1
  }
  if(k==6){
  form=paste('Y ~', paste(paste0('X', 1:(5-omitted)), collapse = ' + '))
  data_gen=data_build_6
  }
  if(k==10){
  form=paste('Y ~', paste(paste0('X', 1:(10-omitted)), collapse = ' + '))
  data_gen=data_build_2
  }
  if(k==20){
  form=paste('Y ~', paste(paste0('X', 1:(20-omitted)), collapse = ' + '))
  data_gen=data_build_3
  }

  data_te=data_gen(5000)
  data_all = data_gen(n)
  report=data.frame()
  prediction=data.frame()
  freq_table=matrix(0, 30, ncol(data_all)-6)

  print('Running time before tree construct:')
  print(proc.time() - start_time)

  for (i in 1:100){
    #print(paste(k,n,i, collapse = ','))
    train_ind <- sample(1:(n), size = n/2)
    data_tr=data_all[train_ind,]
    data_es=data_all[-train_ind,]

    temp=IV_causalTree_mod(form, data_tr, data_es, data_te)
    mse=mean((temp$prediction-data_te$kappa)^2)
    report=rbind(report, c(temp$nodes, mse, temp$mincov, temp$maxcov))

    prediction=rbind(prediction, temp[[2]])
    freq_table=freq_table + temp$split_freq
  }
  names(report)=c('no_leaf', 'mse', 'mincov', 'maxcov')
  temp2=apply(prediction, 2, function(x)quantile(x,c(0.025,0.975)))
  report$cover_prob=mean(data_te$kappa<=temp2[2,] & data_te$kappa>=temp2[1,])
  temp3=apply(prediction, 2, function(x)mean(x))
  report$mse_all = mean((temp3-data_te$kappa)^2)
  report$mse_all
  #print(freq_table)
  #print(colMeans(report))
  #stop()
  report.all = rbind(report.all, colMeans(report))
  names(report.all)=c('no_leaf', 'mse', 'mincov', 'maxcov', 'cover_prob', 'mse')
  print(report.all)

  print('Total running time for one case is:')
  print(proc.time() - start_time)
}
}
}
}

print(report.all)
report.all

# target beat time is about 20 seconds