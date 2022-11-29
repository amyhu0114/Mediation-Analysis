#install.packages("SGL")
#install.packages("glmnet")
#install.packages("gglasso")
library("SGL")
#library("glmnet")
#library("gglasso")
library("Matrix")
library("MASS")

set.seed(2)

########################################################
# Simulation with Multivariate Normal Distribution - cvSGL
########################################################
# n = Number of Observations; p = Number of Predictors; size.groups = 10 means dividing 100 predictors into 10 groups.
n = 50; p = 100; size.groups = 10  
index <- ceiling(1:p / size.groups)


# Covariance Matrix
#-----Case 1 - No correlation between predictors-------#
#10 by 10 block with 1 on diagonal and 0 off diagonals -> Identity Matrix
blocki_10by10_case1 <- diag(10)
blocki_10by10_case1

#-----Case 2 - Weak correlation between predictors-----#
#10 by 10 block with 1 on diagonal and 0.2 off diagonals -> Identity Matrix
blocki_10by10_case2 <- diag(10)
blocki_10by10_case2[blocki_10by10_case2 == 0] <- 0.2
blocki_10by10_case2

#---Case 3 - Moderate correlation between predictors---#
#10 by 10 block with 1 on diagonal and 0.5 off diagonals -> Identity Matrix
blocki_10by10_case3 <- diag(10)
blocki_10by10_case3[blocki_10by10_case3 == 0] <- 0.5
blocki_10by10_case3

#---Case 4 - Strong correlation between predictors---#
#10 by 10 block with 1 on diagonal and 0.8 off diagonals -> Identity Matrix
blocki_10by10_case4 <- diag(10)
blocki_10by10_case4[blocki_10by10_case4 == 0] <- 0.8
blocki_10by10_case4


simulation <- function(blocki_10by10, model_type, a){
	
	#model_type='Lasso'
	# cov is the block digonl matrix with 10 10x10 blocks on the main diagonal
	#blocki_10by10 <- blocki_10by10_case1
	cov <- bdiag(replicate(10, blocki_10by10, simplify=FALSE))
	cov <- as.matrix(cov)
	cov
	cov[cov==NA] <- 0
	mode(cov)

	# X is the nxp design matrix that is generated from normal distribution -> rnorm(x*p, mean=0, sd=#)
	set.seed(2)
	X = matrix(mvrnorm(n, mu=rep(0,100), cov), ncol = p, nrow = n, byrow=TRUE) 

	# Assume only the first 5 predictors actually have non-zero Beta coefficient
	beta = c(1:10, 1:5)

	matEstimates_case <- matrix(NA, nrow=1000, ncol=100)


	optimalLambdas <- rep(NA, 100)
	for (i in 1:1000){
		# y = Response vector with the first 5 non-zero beta coefficients and error
		set.seed(i)
		y = X[,1:15] %*% beta + 0.1*rnorm(n)   

		# data with X matrix and Y vector to enter into the SGL function
		data = list(x = X, y = y)
		
		if(model_type=='cvSGL'){
			cvFit = cvSGL(data, index, alpha=a)
			#names(cvFit)
			# lldiff - An nlam vector of cross validated negative log likelihoods (squared error loss in the linear case, along the regularization path)
			# This gives the average error for each lambda choice -> average cross validation score based on the Beta choices
			# Lower cross-validation score = Less Error <- Desirable
			cvs <- cvFit$lldiff
			cvs
			id.optimal <- which.min(cvs)
			optimalLambdas[i] <- id.optimal
			#id.optimal
			#cvFit$lambdas
			cvFit$lambdas[id.optimal]
			cvFit$fit$beta[,id.optimal]
			#length(cvFit$fit$beta[,id.optimal])
			matEstimates_case[i,] <- cvFit$fit$beta[,id.optimal]
		}
		if(model_type=='GL'){
			#??cvSGL
			cvFit = cvSGL(data, index, type="linear", alpha=0)
			cvs <- cvFit$lldiff
			id.optimal <- which.min(cvs)
			optimalLambdas[i] <- id.optimal
			matEstimates_case[i,] <- cvFit$fit$beta[,id.optimal]
		}
		if(model_type=='Lasso'){
			cvFit = cvSGL(data, index, type="linear", alpha=1)
			cvs <- cvFit$lldiff

			id.optimal <- which.min(cvs)
			optimalLambdas[i] <- id.optimal
			cvFit$lambdas[id.optimal]
			cvFit$fit$beta[,id.optimal]
			matEstimates_case[i,] <- cvFit$fit$beta[,id.optimal]
		}
		

	}
	print(head(matEstimates_case))
	return(matEstimates_case)
}

# cvSGL
# alpha=0.95
matEstimates_case1 <- simulation(blocki_10by10_case1, 'cvSGL', 0.95)
matEstimates_case2 <- simulation(blocki_10by10_case2, 'cvSGL', 0.95)
matEstimates_case3 <- simulation(blocki_10by10_case3, 'cvSGL', 0.95)
matEstimates_case4 <- simulation(blocki_10by10_case4, 'cvSGL', 0.95)

#alpha=0.5
matEstimates_a05_case1 <- simulation(blocki_10by10_case1, 'cvSGL', 0.5)
matEstimates_a05_case2 <- simulation(blocki_10by10_case2, 'cvSGL', 0.5)
matEstimates_a05_case3 <- simulation(blocki_10by10_case3, 'cvSGL', 0.5)
matEstimates_a05_case4 <- simulation(blocki_10by10_case4, 'cvSGL', 0.5)

#alpha=0.25
matEstimates_a025_case1 <- simulation(blocki_10by10_case1, 'cvSGL', 0.25)
matEstimates_a025_case2 <- simulation(blocki_10by10_case2, 'cvSGL', 0.25)
matEstimates_a025_case3 <- simulation(blocki_10by10_case3, 'cvSGL', 0.25)
matEstimates_a025_case4 <- simulation(blocki_10by10_case4, 'cvSGL', 0.25)


# GL
matEstimates_gl_case1 <- simulation(blocki_10by10_case1, 'GL', 0)
matEstimates_gl_case2 <- simulation(blocki_10by10_case2, 'GL', 0)
matEstimates_gl_case3 <- simulation(blocki_10by10_case3, 'GL', 0)
matEstimates_gl_case4 <- simulation(blocki_10by10_case4, 'GL', 0)

## Lasso
matEstimates_lasso_case1 <- simulation(blocki_10by10_case1, 'Lasso', 1)
matEstimates_lasso_case2 <- simulation(blocki_10by10_case2, 'Lasso', 1)
matEstimates_lasso_case3 <- simulation(blocki_10by10_case3, 'Lasso', 1)
matEstimates_lasso_case4 <- simulation(blocki_10by10_case4, 'Lasso', 1)



##########Performance Metrics Evaluation#########
# Evaluates the performance of the prediction by counting the number of correctly identified zero and non-zero predictors and calculating the corresponding percentages.
# metrics is a 1000x6 matrix with...
#		column 1 = Number of Correct Non-Zero Predictors Identified
#		column 2 = Percent of Correct Non-Zero Predictors Identified
#		column 3 = Number of Correct Zero Predictors Identified
#		column 4 = Percent of Correct Zero Predictors Identified
#		column 5 = Number of predictors that are correctly identified
#		column 6 = Percent of predictors that are correctly identified
# Each row of metric represent one simulation (out of 1000)
eval_performance <- function(matEstimates_case){
	metrics <- matrix(NA, nrow=1000, ncol=6)
	numCorrectZero <- 0
	percentCorrectZero <- 0
	numCorrectNonZero <- 0
	percentCorrectNonZero <- 0
	numCorrectOverall <- 0
	percentCorrectOverall <- 0
	#print(head(matEstimates_case))
	for (i in 1:nrow(matEstimates_case)){
		for (j in 1:15){
			#print(matEstimates_case[3,1])
			#print(matEstimates_case[i,j])
			if (matEstimates_case[i,j]!=0){
				numCorrectNonZero <- numCorrectNonZero+1
				}
		}
		for (j in 16:p){
			#print(j)
			#print(matEstimates_case1[i,j])
			if (matEstimates_case[i,j]==0){
				numCorrectZero <- numCorrectZero+1
				}
		}
		percentCorrectNonZero <- (numCorrectNonZero/15.00)*100
		percentCorrectZero <- (numCorrectZero/85.00)*100
		numCorrectOverall <- numCorrectNonZero+numCorrectZero
		percentCorrectOverall <- ((numCorrectNonZero+numCorrectZero)/100)*100
		metrics[i,1] <- numCorrectNonZero
		metrics[i,2] <- percentCorrectNonZero
		metrics[i,3] <- numCorrectZero
		metrics[i,4] <- percentCorrectZero
		metrics[i,5] <- numCorrectOverall
		metrics[i,6] <- percentCorrectOverall
		numCorrectZero <- 0
		percentCorrectZero <- 0
		numCorrectNonZero <- 0
		percentCorrectNonZero <- 0
		numCorrectOverall <- 0
		percentCorrectOverall <- 0
	}

	#print(head(metrics))
	return(metrics)
}

#cvSGL
# alpha=0.95
metrics_case1 <- eval_performance(matEstimates_case1)
metrics_case2 <- eval_performance(matEstimates_case2)
metrics_case3 <- eval_performance(matEstimates_case3)
metrics_case4 <- eval_performance(matEstimates_case4)

# alpha=0.5
metrics_a05_case1 <- eval_performance(matEstimates_a05_case1)
metrics_a05_case2 <- eval_performance(matEstimates_a05_case2)
metrics_a05_case3 <- eval_performance(matEstimates_a05_case3)
metrics_a05_case4 <- eval_performance(matEstimates_a05_case4)

# alpha=0.25
metrics_a025_case1 <- eval_performance(matEstimates_a025_case1)
metrics_a025_case2 <- eval_performance(matEstimates_a025_case2)
metrics_a025_case3 <- eval_performance(matEstimates_a025_case3)
metrics_a025_case4 <- eval_performance(matEstimates_a025_case4)



## GL
metrics_gl_case1 <- eval_performance(matEstimates_gl_case1)
metrics_gl_case2 <- eval_performance(matEstimates_gl_case2)
metrics_gl_case3 <- eval_performance(matEstimates_gl_case3)
metrics_gl_case4 <- eval_performance(matEstimates_gl_case4)

## Lasso
metrics_lasso_case1 <- eval_performance(matEstimates_lasso_case1)
metrics_lasso_case2 <- eval_performance(matEstimates_lasso_case2)
metrics_lasso_case3 <- eval_performance(matEstimates_lasso_case3)
metrics_lasso_case4 <- eval_performance(matEstimates_lasso_case4)


head(matEstimates_lasso_case2)


# Evaluates the performance of the prediction by calculating the average number of correctly identified zero and non-zero predictors across all 1000 simulations and calculating the corresponding percentages.
# metric_average_case1 is a 6x1 vector with...
#		row 1 = Average Number of Correct Non-Zero Predictors Identified
#		row 2 = Average Percent of Correct Non-Zero Predictors Identified
#		row 3 = Average Number of Correct Zero Predictors Identified
#		row 4 = Average Percent of Correct Zero Predictors Identified
#		row 5 = Average Number of predictors that are correctly identified
#		row 6 = Average Percent of predictors that are correctly identified
metric_averages <- function(matEstimates_case){
	metric_avg <- rep(0, 6)
	averageNumCorrectNonZero <- mean(matEstimates_case[,1])
	averagePercentCorrectNonZero <- mean(matEstimates_case[,2])
	averageNumCorrectZero <- mean(matEstimates_case[,3])
	averagePercentCorrectZero <- mean(matEstimates_case[,4])
	averageNumCorrectOverall <- mean(matEstimates_case[,5])
	averagePercentCorrectOverall <- mean(matEstimates_case[,6])
	
	metric_avg[1] = averageNumCorrectNonZero
	metric_avg[2] = averagePercentCorrectNonZero
	metric_avg[3] = averageNumCorrectZero
	metric_avg[4] = averagePercentCorrectZero
	metric_avg[5] = averageNumCorrectOverall
	metric_avg[6] = averagePercentCorrectOverall
	#print(averageNumCorrectNonZero)
	#print(averagePercentCorrectNonZero)
	#print(averageNumCorrectZero)
	#print(averagePercentCorrectZero)
	#print(averageNumCorrectOverall)
	#print(averagePercentCorrectOverall)
	print(metric_avg)
	return(metric_avg)
}

# cvSGL
# alpha=0.95
metric_average_case1 <- metric_averages(metrics_case1)
metric_average_case2 <-metric_averages(metrics_case2)
metric_average_case3 <-metric_averages(metrics_case3)
metric_average_case4 <-metric_averages(metrics_case4)

# alpha=0.5
metric_average_a05_case1 <- metric_averages(metrics_a05_case1)
metric_average_a05_case2 <- metric_averages(metrics_a05_case2)
metric_average_a05_case3 <- metric_averages(metrics_a05_case3)
metric_average_a05_case4 <- metric_averages(metrics_a05_case4)

# alpha=0.25
metric_average_a025_case1 <- metric_averages(metrics_a025_case1)
metric_average_a025_case2 <- metric_averages(metrics_a025_case2)
metric_average_a025_case3 <- metric_averages(metrics_a025_case3)
metric_average_a025_case4 <- metric_averages(metrics_a025_case4)


## GL
metric_average_gl_case1 <- metric_averages(metrics_gl_case1)
metric_average_gl_case2 <- metric_averages(metrics_gl_case2)
metric_average_gl_case3 <- metric_averages(metrics_gl_case3)
metric_average_gl_case4 <- metric_averages(metrics_gl_case4)

## Lasso
metric_average_lasso_case1 <- metric_averages(metrics_lasso_case1)
metric_average_lasso_case2 <- metric_averages(metrics_lasso_case2)
metric_average_lasso_case3 <- metric_averages(metrics_lasso_case3)
metric_average_lasso_case4 <- metric_averages(metrics_lasso_case4)


# Evaluates the empirical bias of the 15 non-zero betas and 1 combined zero beta
#	beta_bias is the mean of all 1000 rounds of estimated beta values subtracted by the true beta value
# 		Take difference for each row -> Absolute value the difference -> Take Average of the differences
# 		For true zero, average over all 0's -> one variable
bias_eval <- function(matEstimates_case){
	beta_bias <- rep(0, 16)
	for (i in 1:15){
		beta_bias[i] <- mean(abs(matEstimates_case[,i] - beta[i]))
		}
	for (i in 16:100){
		beta_bias[16] = beta_bias[16] + mean(abs(matEstimates_case[,i] - 0))
		}
	beta_bias[16] = beta_bias[16]/85	
	print(beta_bias)
	return(beta_bias)
}

# cvSGL
beta_bias_case1 <- bias_eval(matEstimates_case1)
beta_bias_case2 <- bias_eval(matEstimates_case2)
beta_bias_case3 <- bias_eval(matEstimates_case3)
beta_bias_case4 <- bias_eval(matEstimates_case4)

## GL
beta_bias_gl_case1 <- bias_eval(matEstimates_gl_case1)
beta_bias_gl_case2 <- bias_eval(matEstimates_gl_case2)
beta_bias_gl_case3 <- bias_eval(matEstimates_gl_case3)
beta_bias_gl_case4 <- bias_eval(matEstimates_gl_case4)

## Lasso
beta_bias_lasso_case1 <- bias_eval(matEstimates_lasso_case1)
beta_bias_lasso_case2 <- bias_eval(matEstimates_lasso_case2)
beta_bias_lasso_case3 <- bias_eval(matEstimates_lasso_case3)
beta_bias_lasso_case4 <- bias_eval(matEstimates_lasso_case4)

# Evaluates the empirical standard deviation of the 15 non-zero betas and 1 combined zero beta
# 	beta_bias_sd is the standard deviation of all 1000 differences between the true beta value and the estimated beta value
bias_sd_eval <- function(matEstimates_case){
	beta_diff <- matrix(0, nrow=1000, ncol=15)
	beta_bias_sd <- rep(0, 16)
	zero <- c(0)
	for (i in 1:15){
		beta_diff[,i] <- matEstimates_case[,i] - beta[i]
		beta_bias_sd[i] <- sd(beta_diff[,i])
		}
	for (i in 16:100){
		zero <- zero + matEstimates_case[,i]
		}
	beta_bias_sd[16] = sd(zero)
	head(beta_diff)
	print(beta_bias_sd)
	return(beta_bias_sd)
}
# cvSGL
beta_bias_sd_case1 <- bias_sd_eval(matEstimates_case1)
beta_bias_sd_case2 <- bias_sd_eval(matEstimates_case2)
beta_bias_sd_case3 <- bias_sd_eval(matEstimates_case3)
beta_bias_sd_case4 <- bias_sd_eval(matEstimates_case4)

## GL
beta_bias_sd_gl_case1 <- bias_sd_eval(matEstimates_gl_case1)
beta_bias_sd_gl_case2 <- bias_sd_eval(matEstimates_gl_case2)
beta_bias_sd_gl_case3 <- bias_sd_eval(matEstimates_gl_case3)
beta_bias_sd_gl_case4 <- bias_sd_eval(matEstimates_gl_case4)

## Lasso
beta_bias_sd_lasso_case1 <- bias_sd_eval(matEstimates_lasso_case1)
beta_bias_sd_lasso_case2 <- bias_sd_eval(matEstimates_lasso_case2)
beta_bias_sd_lasso_case3 <- bias_sd_eval(matEstimates_lasso_case3)
beta_bias_sd_lasso_case4 <- bias_sd_eval(matEstimates_lasso_case4)





