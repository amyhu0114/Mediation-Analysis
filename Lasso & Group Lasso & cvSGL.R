install.packages("SGL")
library("SGL")
library("Matrix")
library("MASS")

set.seed(2)

########################################################
# Simulation with Multivariate Normal Distribution - cvSGL
# Author - Amy Hu
# Date - November 5 2022
########################################################
# n = Number of Observations; p = Number of Predictors; size.groups = 10 means dividing 100 predictors into 10 groups.
n = 50; p = 100; size.groups = 10  
index <- ceiling(1:p / size.groups)


# Covariance Matrix
#-----Case 1 - No correlation between predictors-------#
#10 by 10 block with 1 on diagonal and 0 off diagonals -> Identity Matrix
blocki_10by10_case1 <- diag(10)

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


simulation <- function(blocki_10by10, model_type){
	# cov is the block digonl matrix with 10 10x10 blocks on the main diagonal
	cov <- bdiag(replicate(10, blocki_10by10, simplify=FALSE)) 
	cov <- as.matrix(cov)
	cov[cov==NA] <- 0

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
			cvFit = cvSGL(data, index, type="linear")
			#names(cvFit)
			# lldiff - An nlam vector of cross validated negative log likelihoods (squared error loss in the linear case, along the regularization path)
			# This gives the average error for each lambda choice -> average cross validation score based on the Beta choices
			# Lower cross-validation score = Less Error <- Desirable
			cvs <- cvFit$lldiff

			id.optimal <- which.min(cvs)
			optimalLambdas[i] <- id.optimal
			#id.optimal
			#cvFit$lambdas
			cvFit$lambdas[id.optimal]
			cvFit$fit$beta[,id.optimal]
			#length(cvFit$fit$beta[,id.optimal])
			matEstimates_case[i,] <- cvFit$fit$beta[,id.optimal]
		}
		if(model_type=='SGL'){
			cvFit = SGL(data, index, type="linear")
			# Figure out how to select optimal lambda
		}
		if(model_type=='Lasso'){
			cvFit = SGL(data, index, type="linear", alpha=1)
			# Figure out how to select optimal lambda
		}
		


	}
	print(head(matEstimates_case))
	return(matEstimates_case)
}

# cvSGL
matEstimates_case1 <- simulation(blocki_10by10_case1, 'cvSGL')
matEstimates_case2 <- simulation(blocki_10by10_case2, 'cvSGL')
matEstimates_case3 <- simulation(blocki_10by10_case3, 'cvSGL')
matEstimates_case4 <- simulation(blocki_10by10_case4, 'cvSGL')

print(head(matEstimates_case1))

## SGL - Not yet working -> How to find optimal lambda and alpha?
#matEstimates_sgl_case1 <- simulation(blocki_10by10_case1, 'SGL')
#matEstimates_sgl_case2 <- simulation(blocki_10by10_case2, 'SGL')
#matEstimates_sgl_case3 <- simulation(blocki_10by10_case3, 'SGL')
#matEstimates_sgl_case4 <- simulation(blocki_10by10_case4, 'SGL')

## Lasso - Not yet working -> How to find optimal lambda and alpha?
#matEstimates_lasso_case1 <- simulation(blocki_10by10_case1, 'Lasso')
#matEstimates_lasso_case2 <- simulation(blocki_10by10_case2, 'Lasso')
#matEstimates_lasso_case3 <- simulation(blocki_10by10_case3, 'Lasso')
#matEstimates_lasso_case4 <- simulation(blocki_10by10_case4, 'Lasso')



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
	for (i in 1:nrow(matEstimates_case1)){
		for (j in 1:15){
			#print(matEstimates_case1[i,j])
			if (matEstimates_case[i,j]!=0){
				numCorrectNonZero <- numCorrectNonZero+1
				}
		}
		for (j in 16:p){
			#print(matEstimates_case1[i,j])
			if (abs(matEstimates_case1[i,j])<=0.01){
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

	print(head(metrics))
	return(metrics)
}

#cvSGL
metrics_case1 <- eval_performance(matEstimates_case1)
metrics_case2 <- eval_performance(matEstimates_case2)
metrics_case3 <- eval_performance(matEstimates_case3)
metrics_case4 <- eval_performance(matEstimates_case4)

## SGL - Not yet working
#metrics_sgl_case1 <- eval_performance(matEstimates_sgl_case1)
#metrics_sgl_case2 <- eval_performance(matEstimates_sgl_case2)
#metrics_sgl_case3 <- eval_performance(matEstimates_sgl_case3)
#metrics_sgl_case4 <- eval_performance(matEstimates_sgl_case4)

## Lasso - Not yet working
#metrics_lasso_case1 <- eval_performance(matEstimates_lasso_case1)
#metrics_lasso_case2 <- eval_performance(matEstimates_lasso_case2)
#metrics_lasso_case3 <- eval_performance(matEstimates_lasso_case3)
#metrics_lasso_case4 <- eval_performance(matEstimates_lasso_case4)




metric_averages <- function(matEstimates_case){
	averageNumCorrectNonZero <- mean(matEstimates_case[,1])
	averagePercentCorrectNonZero <- mean(matEstimates_case[,2])
	averageNumCorrectZero <- mean(matEstimates_case[,3])
	averagePercentCorrectZero <- mean(matEstimates_case[,4])
	averageNumCorrectOverall <- mean(matEstimates_case[,5])
	averagePercentCorrectOverall <- mean(matEstimates_case[,6])
	
	print(averageNumCorrectNonZero)
	print(averagePercentCorrectNonZero)
	print(averageNumCorrectZero)
	print(averagePercentCorrectZero)
	print(averageNumCorrectOverall)
	print(averagePercentCorrectOverall)
}

# cvSGL
metric_averages(metrics_case1)
metric_averages(metrics_case2)
metric_averages(metrics_case3)
metric_averages(metrics_case4)

## SGL - Not yet working
#metric_averages(metrics_sgl_case1)
#metric_averages(metrics_sgl_case2)
#metric_averages(metrics_sgl_case3)
#metric_averages(metrics_sgl_case4)

## Lasso - Not yet working
#metric_averages(metrics_lasso_case1)
#metric_averages(metrics_lasso_case2)
#metric_averages(metrics_lasso_case3)
#metric_averages(metrics_lasso_case4)



