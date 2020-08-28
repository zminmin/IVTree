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


data_build_01 <- function(n, m, cp1, cp2){
# output form: (X1 - Xm, IV, T, Y)

	dt <- rnorm(n * m, mean=0, sd=sqrt(1))
	dim(dt) <- c(n, m)
	dt <- data.frame(dt)
	names(dt) <- c(paste0('X', 1:m))	
	dt$eta <- rowSums(dt[,1:round(0.8*m)])

# linear
	# dt$kappa <- dt$X2
# symmetric 
	# dt$kappa <- (dt$X2>0)*(-dt$X2+1) + (dt$X2<0)*(dt$X2+1)
# asymmetrical
	# dt$kappa <- (dt$X2>0)*(-dt$X2+1) + (dt$X2<0)*(dt$X2+2)
# sin
	# dt$kappa <- sin(pi*dt$X2)
# cos
	# dt$kappa <- cos(pi*dt$X2)
# cos2
	dt$kappa <- cos(2*pi*dt$X2)

	dt$phi <- rnorm(n, mean=0, sd=sqrt(1))
	dt$IV <- rnorm(n)
	dt$T <- cp1*dt$IV + cp2*dt$phi + sqrt(2/pi-cp1^2-cp2^2)*rnorm(n)
	dt$T <- (sign(dt$T)+1)/2
	dt$Y <- dt$eta + dt$T * dt$kappa + dt$phi + rnorm(n, sd = 0.1)
	dt1 <- centering(dt[,1:m], dt$Y, dt$T, dt$IV)
	dt1$kappa = dt$kappa
	dt1$T1 = dt$T
	dt1$IV = dt$IV
	dt1

	return(dt1)
}


data_build_02 <- function(n, m, cp1, cp2){
# output form: (X1 - Xm, IV, T, Y)

	dt <- rnorm(n * m, mean=0, sd=sqrt(1))
	dim(dt) <- c(n, m)
	dt <- data.frame(dt)
	names(dt) <- c(paste0('X', 1:m))	
	dt$eta <- rowSums(dt[,1:round(0.8*m)])
	
# linear
	# dt$kappa <- dt$X2
# symmetric 
	# dt$kappa <- (dt$X2>0)*(-dt$X2+1) + (dt$X2<0)*(dt$X2+1)
# asymmetrical
	# dt$kappa <- (dt$X2>0)*(-dt$X2+1) + (dt$X2<0)*(dt$X2+2)
# sin
	# dt$kappa <- sin(pi*dt$X2)
# cos
	# dt$kappa <- cos(pi*dt$X2)
# cos2
	dt$kappa <- cos(2*pi*dt$X2)

	dt$phi <- rnorm(n, mean=0, sd=sqrt(1))
	dt$IV <- rnorm(n)
	# only half IV works
	rdom <- rbinom(n, 1, 0.5)
	dt$T <- cp1*dt$IV*rdom + cp2*dt$phi + sqrt(2/pi-cp1^2-cp2^2)*rnorm(n)
	dt$T <- (sign(dt$T)+1)/2
	dt$Y <- dt$eta + dt$T * dt$kappa + dt$phi + rnorm(n, sd = 0.1)
	dt1 <- centering(dt[,1:m], dt$Y, dt$T, dt$IV)
	dt1$kappa = dt$kappa
	dt1$T1 = dt$T
	dt1$IV = dt$IV
	dt1

	return(dt1)
}


data_build_03 <- function(n, m, cp1, cp2){
# output form: (X1 - Xm, IV, T, Y)

	dt <- rnorm(n * m, mean=0, sd=sqrt(1))
	dim(dt) <- c(n, m)
	dt <- data.frame(dt)
	names(dt) <- c(paste0('X', 1:m))	
	dt$eta <- rowSums(dt[,1:round(0.8*m)])

# linear
	# dt$kappa <- dt$X2
# symmetric 
	# dt$kappa <- (dt$X2>0)*(-dt$X2+1) + (dt$X2<0)*(dt$X2+1)
# asymmetrical
	# dt$kappa <- (dt$X2>0)*(-dt$X2+1) + (dt$X2<0)*(dt$X2+2)
# sin
	# dt$kappa <- sin(pi*dt$X2)
# cos
	# dt$kappa <- cos(pi*dt$X2)
# cos2
	dt$kappa <- cos(2*pi*dt$X2)

	dt$phi <- rnorm(n, mean=0, sd=sqrt(1))
	dt$IV <- rnorm(n)

	# IV and T correlated if X2 > 0
	# dt$T <- (dt$X2>0)*cp1*dt$IV + cp2*dt$phi + sqrt(2/pi-cp1^2-cp2^2)*rnorm(n)
	# IV and T correlated if X1 > 0
	dt$T <- (dt$X1>0)*cp1*dt$IV + cp2*dt$phi + sqrt(2/pi-cp1^2-cp2^2)*rnorm(n)
	
	dt$T <- (sign(dt$T)+1)/2
	dt$Y <- dt$eta + dt$T * dt$kappa + dt$phi + rnorm(n, sd = 0.1)
	dt1 <- centering(dt[,1:m], dt$Y, dt$T, dt$IV)
	dt1$kappa = dt$kappa
	dt1$T1 = dt$T
	dt1$IV = dt$IV
	dt1

	return(dt1)
}



grf_relabeling <- function(dt, reduced_form_weight){
# Require: dt = (X1 - Xm, IV, T, Y)
# 			reduced_form_weight \in (0,1)

	average_outcome <- mean(dt$Y)
	average_treatment <- mean(dt$T)
	average_IV <- mean(dt$IV)
	average_regularized_IV <- (1-reduced_form_weight)*average_IV + reduced_form_weight*average_treatment

	regularized_IV <- (1-reduced_form_weight) * dt$IV + reduced_form_weight * dt$T

	numerator <- sum((regularized_IV - average_regularized_IV)*(dt$Y - average_outcome))
	denominator <- sum((regularized_IV - average_regularized_IV)*(dt$T - average_treatment))

	if(abs(denominator)<0.0001){
		print("zero denominator in grf_relabeling")
		return(0)
	}

	local_average_treatment_effect <- numerator / denominator
	residual <- (dt$Y - average_outcome) - local_average_treatment_effect * (dt$T - average_treatment)
	dt$rho <- (regularized_IV - average_regularized_IV) * residual

	return(dt)
}


cutting_point <- function(dt, m){
# Require: dt = (X1 - Xm, IV, T, Y)
	

	dt <- grf_relabeling(dt, 0)

	# get CART parent node MSE
	parent_average_rho <- mean(dt$rho)
	parent_MSE <- sum((dt$rho - parent_average_rho)^2)

	# get IVT parent node treatment effect
	parent_average_outcome <- mean(dt$Y)
	parent_average_treatment <- mean(dt$T)
	parent_average_IV <- mean(dt$IV)
	numerator <- sum((dt$IV - parent_average_IV)*(dt$Y - parent_average_outcome))
	denominator <- sum((dt$IV - parent_average_IV)*(dt$T - parent_average_treatment))
	if(abs(denominator)<0.0001){
		print("zero denominator in parent node")
		return(0)
	}
	parent_treatment_effect <- numerator / denominator


	grf_reduced_value <- 0
	grf_splitting_feature <- -1
	grf_splitting_value <- -1
	ivt_reduced_value <- 0
	ivt_splitting_feature <- -1
	ivt_splitting_value <- -1
	for (i in 1:m){
		pair <- dt[order(dt[,i]),]
		unique_value <- unique(pair[,i])
		break_number <- length(unique_value) - 1
		if(break_number == 0){
			next
		}
		for(j in 1:break_number){
			left <- dt[which(dt[,i] <= unique_value[j]),]
			right <- dt[which(dt[,i] > unique_value[j]),]
			# minimal node is 10
			if(nrow(left) <= 10 || nrow(right) <= 10){
				next
			}

			# get CART children node MSE
			left_average <- mean(left$rho)
			right_average <- mean(right$rho)
			left_MSE <- sum((left$rho - left_average)^2)
			right_MSE <- sum((right$rho - right_average)^2)

			# get CART improvement
			improve <- parent_MSE - left_MSE - right_MSE  
			if(improve > grf_reduced_value){
				grf_reduced_value <- improve
				grf_splitting_feature <- i
				grf_splitting_value <- unique_value[j]
			}

			# get IVT parent node treatment effect
			left_average_outcome <- mean(left$Y)
			left_average_treatment <- mean(left$T)
			left_average_IV <- mean(left$IV)
			numerator <- sum((left$IV - left_average_IV)*(left$Y - left_average_outcome))
			denominator <- sum((left$IV - left_average_IV)*(left$T - left_average_treatment))
			if(abs(denominator)<0.0001){
				print("zero denominator in left node")
				next
			}
			left_treatment_effect <- numerator / denominator

			right_average_outcome <- mean(right$Y)
			right_average_treatment <- mean(right$T)
			right_average_IV <- mean(right$IV)
			numerator <- sum((right$IV - right_average_IV)*(right$Y - right_average_outcome))
			denominator <- sum((right$IV - right_average_IV)*(right$T - right_average_treatment))
			if(abs(denominator)<0.0001){
				print("zero denominator in right node")
				next
			}
			right_treatment_effect <- numerator / denominator			

			# get IVT improvement
			improve <- nrow(left)*(left_treatment_effect)^2 + nrow(right)*(right_treatment_effect)^2 - nrow(dt)*(parent_treatment_effect)^2
			if(improve > ivt_reduced_value){
				ivt_reduced_value <- improve
				ivt_splitting_feature <- i
				ivt_splitting_value <- unique_value[j]
			}
		}
	}


	return(list(grf_splitting_feature, grf_splitting_value, ivt_splitting_feature, ivt_splitting_value))
}



n = 1000
m = 3
cp1 = 0.6
cp2 = 0.5
times = 100
result <- matrix(0, times, 4)
for(i in 1:times){
	print(i)
	te1 = data_build_03(n, m, cp1, cp2)
	res = cutting_point(te1, m)
	for (j in 1:4) {
		result[i,j] = res[[j]]
	}
}
colSums(result[,c(1,3)]==2)



























