install.packages("SGL")
library("SGL")
library("Matrix")
library("MASS")

set.seed(2)

## SGL - This function runs SGL nfold+1 times; the initial run is to find the lambda sequence, subsequent runs are used to compute the cross-validated error rate and its standard deviation

# n = Number of Observations; p = Number of Predictors; size.groups = 10 means dividing 100 predictors into 10 groups.
n = 50; p = 100; size.groups = 10  
index <- ceiling(1:p / size.groups)

# nxp design matrix that is generated from normal distribution -> rnorm(x*p, mean=0, sd=#)
X = matrix(rnorm(n * p), ncol = p, nrow = n) 
X 

# Assume only the first 5 predictors actually have non-zero Beta coefficient
# beta = (-2,-1,0,1,2) -> True beta values??
beta = (-2:2)
# y = Response vector with the first 5 non-zero beta coefficients and error
# Matrix multiplication of the 50x5 matrix (X) and the length 5 vector (beta) 
y = X[,1:5] %*% beta + 0.1*rnorm(n)   
y

# data with X matrix and Y vector to enter into the SGL function
# For type="linear", the input data should be a list with $x$ an input matrix of dimension n-obs by p-vars, and $y$ a length $n$ response vector. For type="logit", the input data should be a list with $x$, an input matrix, as before, and $y$ a length $n$ binary response vector. For type="cox", the input data should be a list with x as before, time, an n-vector corresponding to failure/censor times, and status, an n-vector indicating failure (1) or censoring (0).
data = list(x = X, y = y)
fit = SGL(data, index, type = "linear")

fit
summary(fit)
names(fit)

# Beta = A p by nlam matrix, giving the penalized MLEs for the nlam different models, where the index corresponds to the penalty parameter lambda
fit$beta
# Lambda = Sequence of penalty parameter (lambda) values used
# If a sequence of lambda values for fitting is not inputted by users, SGL self-selects the lambda values at the initial run
fit$lambda




##cvSGL
cvFit = cvSGL(data, index, type="linear")
names(cvFit)
# lldiff - An nlam vector of cross validated negative log likelihoods (squared error loss in the linear case, along the regularization path)
# This gives the average error for each lambda choice -> average cross validation score based on the Beta choices
# Lower cross-validation score = Less Error <- Desirable
cvs <- cvFit$lldiff
cvs

id.optimal <- which.min(cvs)
id.optimal
cvFit$lambdas
cvFit$lambdas[id.optimal]
cvFit$fit$beta[,id.optimal]




################################################
Simulation with Multivariate Normal Distribution
################################################
# n = Number of Observations; p = Number of Predictors; size.groups = 10 means dividing 100 predictors into 10 groups.
n = 50; p = 100; size.groups = 10  
index <- ceiling(1:p / size.groups)


# Covariance Matrix
#10 by 10 block with 1 on diagonal and 0.5 off diagonals
blocki_10by10 <- diag(10)
blocki_10by10[blocki_10by10 == 0] <- 0.5
blocki_10by10

# block digonl matrix with 10 10x10 block of 1 and 0.5 and NA (later changed to 0) off diagonal
cov <- bdiag(replicate(10, blocki_10by10, simplify=FALSE)) 
cov <- as.matrix(cov)
cov[cov==NA] <- 0
is.na(cov)
cov

# nxp design matrix that is generated from normal distribution -> rnorm(x*p, mean=0, sd=#)
set.seed(5)
X = matrix(mvrnorm(n, mu=rep(0,100), cov), ncol = p, nrow = n, byrow=TRUE) 
X

# Assume only the first 5 predictors actually have non-zero Beta coefficient
beta = c(1:10, 1:5)

matEstimates <- matrix(NA, nrow=100, ncol=100)
optimalLambdas <- rep(NA, 100)
for (i in 1:100){
# y = Response vector with the first 5 non-zero beta coefficients and error
set.seed(i)
y = X[,1:15] %*% beta + 0.1*rnorm(n)   
#y

# data with X matrix and Y vector to enter into the SGL function
data = list(x = X, y = y)



cvFit = cvSGL(data, index, type="linear")
names(cvFit)
# lldiff - An nlam vector of cross validated negative log likelihoods (squared error loss in the linear case, along the regularization path)
# This gives the average error for each lambda choice -> average cross validation score based on the Beta choices
# Lower cross-validation score = Less Error <- Desirable
cvs <- cvFit$lldiff
#cvs

id.optimal <- which.min(cvs)
optimalLambdas[i] <- id.optimal
id.optimal
#cvFit$lambdas
cvFit$lambdas[id.optimal]
cvFit$fit$beta[,id.optimal]
#length(cvFit$fit$beta[,id.optimal])
matEstimates[i,] <- cvFit$fit$beta[,id.optimal]
}
matEstimates
dim(matEstimates)
table(optimalLambdas)
optimalLambdas

#Penalization unbiased
#Refit by least square regression

############################
## Case 1: No Correlation ##
############################
# n = Number of Observations; p = Number of Predictors; size.groups = 10 means dividing 100 predictors into 10 groups.
n = 50; p = 100; size.groups = 10  
index <- ceiling(1:p / size.groups)


# Covariance Matrix
#100 by 100 identity matrix
cov <- diag(100)


# nxp design matrix that is generated from normal distribution -> rnorm(x*p, mean=0, sd=#)
set.seed(2)
X = matrix(mvrnorm(n, mu=rep(0,1000), cov), ncol = p, nrow = n, byrow=TRUE) 

# Assume only the first 5 predictors actually have non-zero Beta coefficient
beta = c(1:10, 1:5)

matEstimates <- matrix(NA, nrow=100, ncol=100)
optimalLambdas <- rep(NA, 100)
for (i in 1:100){
  # y = Response vector with the first 5 non-zero beta coefficients and error
  set.seed(i)
  y = X[,1:15] %*% beta + 0.1*rnorm(n)   
  #y
  
  # data with X matrix and Y vector to enter into the SGL function
  data = list(x = X, y = y)
  
  
  
  cvFit = cvSGL(data, index, type="linear")
  names(cvFit)
  # lldiff - An nlam vector of cross validated negative log likelihoods (squared error loss in the linear case, along the regularization path)
  # This gives the average error for each lambda choice -> average cross validation score based on the Beta choices
  # Lower cross-validation score = Less Error <- Desirable
  cvs <- cvFit$lldiff
  #cvs
  
  id.optimal <- which.min(cvs)
  optimalLambdas[i] <- id.optimal
  id.optimal
  #cvFit$lambdas
  cvFit$lambdas[id.optimal]
  cvFit$fit$beta[,id.optimal]
  #length(cvFit$fit$beta[,id.optimal])
  matEstimates[i,] <- cvFit$fit$beta[,id.optimal]
}
matEstimates
dim(matEstimates)
table(optimalLambdas)
optimalLambdas

##############################
## Case 2: Weak Correlation ##
##############################
# n = Number of Observations; p = Number of Predictors; size.groups = 10 means dividing 100 predictors into 10 groups.
n = 50; p = 100; size.groups = 10  
index <- ceiling(1:p / size.groups)


# Covariance Matrix
#10 by 10 block with 1 on diagonal and 0.2 off diagonals
blocki_10by10 <- diag(10)
blocki_10by10[blocki_10by10 == 0] <- 0.2
blocki_10by10

# block digonl matrix with 10 10x10 block of 1 and 0.5 and NA (later changed to 0) off diagonal
cov <- bdiag(replicate(10, blocki_10by10, simplify=FALSE)) 
cov <- as.matrix(cov)
cov[cov==NA] <- 0
is.na(cov)
cov

# nxp design matrix that is generated from normal distribution -> rnorm(x*p, mean=0, sd=#)
set.seed(5)
X = matrix(mvrnorm(n, mu=rep(0,1000), cov), ncol = p, nrow = n, byrow=TRUE) 
X

# Assume only the first 5 predictors actually have non-zero Beta coefficient
beta = c(1:10, 1:5)

matEstimates <- matrix(NA, nrow=100, ncol=100)
optimalLambdas <- rep(NA, 100)
for (i in 1:100){
  # y = Response vector with the first 5 non-zero beta coefficients and error
  set.seed(i)
  y = X[,1:15] %*% beta + 0.1*rnorm(n)   
  #y
  
  # data with X matrix and Y vector to enter into the SGL function
  data = list(x = X, y = y)
  
  
  
  cvFit = cvSGL(data, index, type="linear")
  names(cvFit)
  # lldiff - An nlam vector of cross validated negative log likelihoods (squared error loss in the linear case, along the regularization path)
  # This gives the average error for each lambda choice -> average cross validation score based on the Beta choices
  # Lower cross-validation score = Less Error <- Desirable
  cvs <- cvFit$lldiff
  #cvs
  
  id.optimal <- which.min(cvs)
  optimalLambdas[i] <- id.optimal
  id.optimal
  #cvFit$lambdas
  cvFit$lambdas[id.optimal]
  cvFit$fit$beta[,id.optimal]
  #length(cvFit$fit$beta[,id.optimal])
  matEstimates[i,] <- cvFit$fit$beta[,id.optimal]
}
matEstimates
dim(matEstimates)
table(optimalLambdas)
optimalLambdas

##################################
## Case 3: Moderate Correlation ##
##################################
# n = Number of Observations; p = Number of Predictors; size.groups = 10 means dividing 100 predictors into 10 groups.
n = 50; p = 100; size.groups = 10  
index <- ceiling(1:p / size.groups)


# Covariance Matrix
#10 by 10 block with 1 on diagonal and 0.5 off diagonals
blocki_10by10 <- diag(10)
blocki_10by10[blocki_10by10 == 0] <- 0.5
blocki_10by10

# block digonl matrix with 10 10x10 block of 1 and 0.5 and NA (later changed to 0) off diagonal
cov <- bdiag(replicate(10, blocki_10by10, simplify=FALSE)) 
cov <- as.matrix(cov)
cov[cov==NA] <- 0
is.na(cov)
cov

# nxp design matrix that is generated from normal distribution -> rnorm(x*p, mean=0, sd=#)
set.seed(5)
X = matrix(mvrnorm(n, mu=rep(0,1000), cov), ncol = p, nrow = n, byrow=TRUE) 
X

# Assume only the first 5 predictors actually have non-zero Beta coefficient
beta = c(1:10, 1:5)

matEstimates <- matrix(NA, nrow=100, ncol=100)
optimalLambdas <- rep(NA, 100)
for (i in 1:100){
  # y = Response vector with the first 5 non-zero beta coefficients and error
  set.seed(i)
  y = X[,1:15] %*% beta + 0.1*rnorm(n)   
  #y
  
  # data with X matrix and Y vector to enter into the SGL function
  data = list(x = X, y = y)
  
  
  
  cvFit = cvSGL(data, index, type="linear")
  names(cvFit)
  # lldiff - An nlam vector of cross validated negative log likelihoods (squared error loss in the linear case, along the regularization path)
  # This gives the average error for each lambda choice -> average cross validation score based on the Beta choices
  # Lower cross-validation score = Less Error <- Desirable
  cvs <- cvFit$lldiff
  #cvs
  
  id.optimal <- which.min(cvs)
  optimalLambdas[i] <- id.optimal
  id.optimal
  #cvFit$lambdas
  cvFit$lambdas[id.optimal]
  cvFit$fit$beta[,id.optimal]
  #length(cvFit$fit$beta[,id.optimal])
  matEstimates[i,] <- cvFit$fit$beta[,id.optimal]
}
matEstimates
dim(matEstimates)
table(optimalLambdas)
optimalLambdas
################################
## Case 4: Strong Correlation ##
################################
# n = Number of Observations; p = Number of Predictors; size.groups = 10 means dividing 100 predictors into 10 groups.
n = 50; p = 100; size.groups = 10  
index <- ceiling(1:p / size.groups)


# Covariance Matrix
#10 by 10 block with 1 on diagonal and 0.8 off diagonals
blocki_10by10 <- diag(10)
blocki_10by10[blocki_10by10 == 0] <- 0.8
blocki_10by10

# block digonl matrix with 10 10x10 block of 1 and 0.5 and NA (later changed to 0) off diagonal
cov <- bdiag(replicate(10, blocki_10by10, simplify=FALSE)) 
cov <- as.matrix(cov)
cov[cov==NA] <- 0
is.na(cov)
cov

# nxp design matrix that is generated from normal distribution -> rnorm(x*p, mean=0, sd=#)
set.seed(5)
X = matrix(mvrnorm(n, mu=rep(0,1000), cov), ncol = p, nrow = n, byrow=TRUE) 
X

# Assume only the first 5 predictors actually have non-zero Beta coefficient
beta = c(1:10, 1:5)

matEstimates <- matrix(NA, nrow=100, ncol=100)
optimalLambdas <- rep(NA, 100)
for (i in 1:100){
  # y = Response vector with the first 5 non-zero beta coefficients and error
  set.seed(i)
  y = X[,1:15] %*% beta + 0.1*rnorm(n)   
  #y
  
  # data with X matrix and Y vector to enter into the SGL function
  data = list(x = X, y = y)
  
  
  
  cvFit = cvSGL(data, index, type="linear")
  names(cvFit)
  # lldiff - An nlam vector of cross validated negative log likelihoods (squared error loss in the linear case, along the regularization path)
  # This gives the average error for each lambda choice -> average cross validation score based on the Beta choices
  # Lower cross-validation score = Less Error <- Desirable
  cvs <- cvFit$lldiff
  #cvs
  
  id.optimal <- which.min(cvs)
  optimalLambdas[i] <- id.optimal
  id.optimal
  #cvFit$lambdas
  cvFit$lambdas[id.optimal]
  cvFit$fit$beta[,id.optimal]
  #length(cvFit$fit$beta[,id.optimal])
  matEstimates[i,] <- cvFit$fit$beta[,id.optimal]
}
matEstimates
dim(matEstimates)
table(optimalLambdas)
optimalLambdas
