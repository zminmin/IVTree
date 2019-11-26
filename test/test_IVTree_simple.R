# This program tests out the features of the IVTree package

library(devtools)
library(rpart)
library(rpart.plot)
library(reshape2)
library(plyr)
library(grf)


library(IVTree)


centering <- function(X, Y, W, Z) {

  forest.Y <- regression_forest(X, Y)
  Y.hat = predict(forest.Y)$predictions
          
  forest.W <- regression_forest(X, W)
  W.hat = predict(forest.W)$predictions

  forest.Z <- regression_forest(X, Z)
  Z.hat = predict(forest.Z)$predictions
          
  data <- cbind(X, Y = Y - Y.hat, T = W - W.hat, IV = Z - Z.hat)
}

