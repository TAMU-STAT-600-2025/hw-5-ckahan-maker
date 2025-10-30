
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)
library(testthat)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################
test_that("soft_c matches soft for several scalar inputs", {
  expect_equal(soft_c(4, 1),  soft(4, 1))
  expect_equal(soft_c(4, 5),  soft(4, 5))
  expect_equal(soft_c(-4, 3), soft(-4, 3))
})
# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################
test_that("lasso_c matches lasso for multiple inputs", {
  lambda <- 2
  
  # Example #1:
  set.seed(123)
  X1 <- matrix(rnorm(500), 100, 5)
  Y1 <- rep(0, 100)
  beta1 <- sample(1:10, 5)
  expect_equal(lasso(X1, Y1, beta1, lambda),
               lasso_c(X1, Y1, beta1, lambda),
               tolerance = 1e-8)
  
  # Example #2:
  set.seed(456)
  X2 <- matrix(rnorm(100), 50, 2)
  Y2 <- rep(50, 50)
  beta2 <- sample(-3:3, 2)
  expect_equal(lasso(X2, Y2, beta2, lambda),
               lasso_c(X2, Y2, beta2, lambda),
               tolerance = 1e-8)
  
  # Example #3:
  set.seed(789)
  X3 <- matrix(rnorm(1000), 10, 100)
  Y3 <- rnorm(10)
  beta3 <- -49:50  
  expect_equal(lasso(X3, Y3, beta3, lambda),
               lasso_c(X3, Y3, beta3, lambda),
               tolerance = 1e-8)
})

# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################

# Tests on riboflavin data
##########################
# require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
# data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression
# 
# # Make sure riboflavin$x is treated as matrix later in the code for faster computations
# class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]
# 
# # Standardize the data
# out <- standardizeXY(riboflavin$x, riboflavin$y)
# 
# # This is just to create lambda_seq, can be done faster, but this is simpler
# outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)
# 
# # The code below should assess your speed improvement on riboflavin data
# microbenchmark(
#   fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
#   fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
#   times = 10
# )
