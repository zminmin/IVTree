## This is to generate simulation data


library(grf)


centering <- function(X, Y, W, Z){

  forest.Y <- regression_forest(X, Y)
  Y.hat = predict(forest.Y)$predictions
          
  forest.W <- regression_forest(X, W)
  W.hat = predict(forest.W)$predictions

  forest.Z <- regression_forest(X, Z)
  Z.hat = predict(forest.Z)$predictions
          
  data <- cbind(X, Y = Y - Y.hat, T = W - W.hat, IV = Z - Z.hat)
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




cp1 = 0.6
omitted = T


for (s in 1:10){
  for (cp2 in c(0.5)){
    for (k in c(5)){
      for (n in c(1000, 2000, 3000, 4000, 5000)){

        set.seed(s)


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

        filename <- paste("CONFIG", s, cp1, cp2, k, n, sep="-"))

        save(data_te, data_all, file = "study1.RData")



        # report=data.frame()
        # prediction=data.frame()
        # freq_table=matrix(0, 30, ncol(data_all)-6)


        # for (i in 1:100){
        #   train_ind <- sample(1:(n), size = n/2)
        #   data_tr=data_all[train_ind,]
        #   data_es=data_all[-train_ind,]

        #   temp=IV_causalTree_mod(form, data_tr, data_es, data_te)
        #   mse=mean((temp$prediction-data_te$kappa)^2)
        #   report=rbind(report, c(temp$nodes, mse, temp$mincov, temp$maxcov))

        #   prediction=rbind(prediction, temp[[2]])
        #   freq_table=freq_table + temp$split_freq
        # }
        # names(report)=c('no_leaf', 'mse', 'mincov', 'maxcov')
        # temp2=apply(prediction, 2, function(x)quantile(x,c(0.025,0.975)))
        # report$cover_prob=mean(data_te$kappa<=temp2[2,] & data_te$kappa>=temp2[1,])
        # temp3=apply(prediction, 2, function(x)mean(x))
        # report$mse_all = mean((temp3-data_te$kappa)^2)
        # report$mse_all


      }
    }
  }
}