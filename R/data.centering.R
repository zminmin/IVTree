
#' data.centering function
#' 
#' library(grf)
#' @importFrom grf regression_forest
#' @importFrom stats predict

data.centering <- function(X, Y, W, Z){

	forest.Y <- regression_forest(X, Y)
	Y.hat = predict(forest.Y)$predictions
	        
	forest.W <- regression_forest(X, W)
	W.hat = predict(forest.W)$predictions

	forest.Z <- regression_forest(X, Z)
	Z.hat = predict(forest.Z)$predictions
        
	data <- cbind(X, Y = Y - Y.hat, T = W - W.hat, IV = Z - Z.hat, T1 = W)

	return(data)
}