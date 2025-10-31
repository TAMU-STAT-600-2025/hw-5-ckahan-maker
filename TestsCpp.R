
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)
library(testthat)
library(microbenchmark)

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
test_that("fitLASSOstandardized_c agrees with fitLASSOstandardized across random inputs", {
  lambda <- 0.1

  # Example 1:
  set.seed(101)
  Xtilde1 <- matrix(rnorm(40), nrow = 8, ncol = 5)
  Ytilde1 <-Xtilde1[, 2] + rnorm(8, sd = 10)
  centeredXY <- standardizeXY(Xtilde1, Ytilde1)
  cXtilde1 <- centeredXY$Xtilde
  cYtilde1 <- centeredXY$Ytilde
  beta_start1 <- rep(1, 5)
  
  expect_equal(
    as.vector(fitLASSOstandardized_c( cXtilde1, cYtilde1, lambda, beta_start1)),
    fitLASSOstandardized(cXtilde1, cYtilde1, lambda, beta_start1)$beta,
    tolerance = 1e-8
  )
  
  # Example 2:
  set.seed(202)
  Xtilde2 <- matrix(rnorm(72), nrow = 12, ncol = 6)
  Ytilde2 <- 2*Xtilde1[, 1] - 0.5*Xtilde1[, 4] + rnorm(12, sd = 5)
  centeredXY <- standardizeXY(Xtilde2, Ytilde2)
  cXtilde2 <- centeredXY$Xtilde
  cYtilde2 <- centeredXY$Ytilde
  beta_start2 <- rep(0, 6)
  
  expect_equal(
    as.vector(fitLASSOstandardized_c(cXtilde2, cYtilde2, lambda, beta_start2)),
    fitLASSOstandardized(cXtilde2, cYtilde2, lambda, beta_start2)$beta,
    tolerance = 1e-8
  )
  
  # Example 3:
  set.seed(303)
  Xtilde3 <- matrix(rnorm(90), nrow = 15, ncol = 6)
  Ytilde3 <- rowSums(Xtilde3) + rnorm(15, sd = 1)
  centeredXY <- standardizeXY(Xtilde3, Ytilde3)
  cXtilde3 <- centeredXY$Xtilde
  cYtilde3 <- centeredXY$Ytilde
  beta_start3 <- rep(-1, 6)
  
  expect_equal(
    as.vector(fitLASSOstandardized_c(cXtilde3, cYtilde3, lambda, beta_start3)),
    fitLASSOstandardized(cXtilde3, cYtilde3, lambda, beta_start3)$beta,
    tolerance = 1e-8
  )
})

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################
Xtilde <- matrix(rnorm(72), nrow = 12, ncol = 6)
Ytilde <- 2 * Xtilde[, 1] - 0.5 * Xtilde[, 4] + rnorm(12, sd = 5)
centeredXY <- standardizeXY(Xtilde, Ytilde)
cXtilde <- centeredXY$Xtilde
cYtilde <- centeredXY$Ytilde
beta_start <- rep(0, 6)
lambda <- 1

# benchmark
microbenchmark(
  R_version = fitLASSOstandardized(cXtilde, cYtilde, lambda, beta_start),
  Cpp_version = fitLASSOstandardized_c(cXtilde, cYtilde, lambda, beta_start),
  times = 100
)
# R version: 487 microseconds median
# C Version: 7.60 microseconds median


# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################
test_that("fitLASSOstandardized_seq matches fitLASSOstandardized_seq_c", {
  # Example #1
  Xtilde1 <- matrix(rnorm(40), nrow = 8, ncol = 5)
  Ytilde1 <- rnorm(8)
  centeredXY <- standardizeXY(Xtilde1, Ytilde1)
  cXtilde1 <- centeredXY$Xtilde
  cYtilde1 <- centeredXY$Ytilde
  lambda_seq1 <- c(1, 0.1, 0.01, 0.001)
  
  beta_mat_r1 <- fitLASSOstandardized_seq(cXtilde1, cYtilde1, lambda_seq = lambda_seq1, eps = 0.001)$beta_mat
  beta_mat_c1 <- fitLASSOstandardized_seq_c(cXtilde1, cYtilde1, lambda_seq1, eps = 0.001)
  
  expect_equal(beta_mat_r1, beta_mat_c1, tolerance = 1e-8)
  
  # Example 2
  Xtilde2 <- matrix(rnorm(200), nrow = 20, ncol = 10)
  Xtilde2[, 5:10] <- Xtilde2[, 1:6] + rnorm(120, sd = 0.1)  
  Ytilde2 <- Xtilde2[, 1] - 0.5 * Xtilde2[, 2] + rnorm(20)
  centeredXY <- standardizeXY(Xtilde2, Ytilde2)
  cXtilde2 <- centeredXY$Xtilde
  cYtilde2 <- centeredXY$Ytilde
  lambda_seq2 <- c(1, 0.5, 0.2, 0.1, 0.05)
  
  beta_mat_r2 <- fitLASSOstandardized_seq(cXtilde2, cYtilde2, lambda_seq = lambda_seq2, eps = 1e-4)$beta_mat
  beta_mat_c2 <- fitLASSOstandardized_seq_c(cXtilde2, cYtilde2, lambda_seq2, 1e-4)
  
  expect_equal(beta_mat_r2, beta_mat_c2, tolerance = 1e-8)
})
# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################
set.seed(123)
Xtilde_bench <- matrix(rnorm(5000), nrow = 100, ncol = 50)
Xtilde_bench[, 26:50] <- Xtilde_bench[, 1:25] + rnorm(2500, sd = 0.05)
Ytilde_bench <- Xtilde_bench[, 1] - 0.3 * Xtilde_bench[, 45] + 0.5 * Xtilde_bench[, 27] + rnorm(100)

# Standardize
centeredXY <- standardizeXY(Xtilde_bench, Ytilde_bench)
cXtilde_bench <- centeredXY$Xtilde
cYtilde_bench <- centeredXY$Ytilde

# Lambda sequence
lambda_seq_bench <- c(2, 0.5, 0.3, 0.15, 0.01)

# Benchmark R vs C++ 
bench_result <- microbenchmark(
  R_version = fitLASSOstandardized_seq(cXtilde_bench, cYtilde_bench,
                                       lambda_seq = lambda_seq_bench,
                                       eps = 1e-4)$beta_mat,
  C_version = fitLASSOstandardized_seq_c(cXtilde_bench, cYtilde_bench,
                                         lambda_seq_bench, 1e-4),
  times = 50
)
# R-version: 6460 microseconds median
# C++-version: 144 microseconds median

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
