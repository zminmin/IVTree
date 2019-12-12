#'
#' 
#' 
#'
#' @title ...
#' @description ...

#' @param ... ....
#' @param ... ....
#' @return ....
#' @export
#' @useDynLib IVForest, .registration = TRUE, .fixes = "C_"


honest.IVTree <- function(formula, data, weights, treatment, treatment1, IV, subset, data_te, 
							  est_data, est_weights, est_treatment, est_treatment1, est_IV, est_subset,
							  na.action = na.IVTree, split.Honest, forest_size, if_center,
							  HonestSampleSize, split.Bucket, bucketNum = 10,
							  bucketMax = 40, cv.Honest, cv.option, minsize = 2L, model = FALSE,
							  x = FALSE, y = TRUE, propensity, control, split.alpha = 0.5, 
							  cv.alpha = 0.5, cv.gamma = 0.5, split.gamma = 0.5, cost, ...) {




	if (missing(if_center)) {
		if_center <- TRUE
		warning("The default if_center = TRUE for your forest.")
	}

	## check the if_center == T/F
	if_center.num <- pmatch(if_center, c(T, F))
	if(is.na(if_center.num)) 
		stop("Invalid if_center input, if_center can be only TRUE or FALSE.")

	if (if_center) {
		#============ if the original data is not centered, center ===================
		# there are four X, one Y, one instrumental variable IV and one W named as T1
		data_all = data.centering(data_all[,1:4], data_all$Y, data_all$T1, data_all$IV) 
		#============ ******************************************** ===================
	}






}



# ============ if the original data is not centered, center ===================
# there are four X, one Y, one instrumental variable IV and one W named as T1
# data_all = data.centering(data_all[,1:4], data_all$Y, data_all$T1, data_all$IV) 
# ============ ******************************************** ===================




