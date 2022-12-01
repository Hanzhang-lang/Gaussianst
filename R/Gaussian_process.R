Rcpp::sourceCpp("src/Gaussian_arma.cpp")

#'
#'@title One dimensional gaussian process regression
#'@description naive implementation for gaussian process regression
#' @param X input matrix
#' @param y input value
#' @param Xtest predict matrix
#' @param alpha error parameter
#' @param kernel kernel used in gaussian process, by default we use exponential cov
#' @param theta parameters for kernel
#' @param plot_f plot f_prior or not
#'@importFrom graphics abline lines
#'@importFrom stats na.omit rnorm
#' @return a list containing posterior mean, std, posterior sample and marginal likelihood
#' @export
#'
#' @examples
#' x <- seq(-20,20,0.1)
#' y <- sin(x)/x + rnorm(401,sd=0.03)
#' res = Gaussian_reg_fit_pre(x, y , x, alpha=1e-5, theta=list(1, 1/(3.88)^2), plot_f=TRUE)
#' plot(x, y, type ="l")
#' lines(x, res$mu, col="red")
#'
Gaussian_reg_fit_pre <- function(X , y, Xtest, alpha=1e-10, kernel=exponential_cov, theta, plot_f=FALSE){
# preprocessing: na.action == na.omit()
df= na.omit(data.frame(y, X))
y = df[,1]
X = as.matrix(df[, -1])
n = length(X)

K_ss = kernel(Xtest, Xtest, theta)
L_ss = chol(K_ss +  alpha* diag(length(Xtest)))

f_prior = L_ss %*% rnorm(length(Xtest))
#plot f_prior
if (plot_f){
plot(Xtest, f_prior, type='l')
abline(h=0, col='gray')}

K = kernel(X, X, theta)
L = t(chol(K + alpha * diag(n)))

K_s = kernel(X, Xtest, theta)
Lk = solve(L, K_s)
v = solve(L, y)

# compute marginal likelihood(simple output)
marginal = -0.5 * sum(v^2) - sum(log(diag(L))) - n/2 * log(2 * pi)

mu = t(Lk) %*% v
s2 = diag(K_ss) - apply(Lk^2, 2, sum)
stdv = sqrt(s2)
L_new = chol(K_ss + 1e-10 * diag(length(Xtest)) - t(Lk) %*% Lk)
f_post = mu + t(L_new) %*% rnorm(length(Xtest))
return (list(
  "mu" = mu,
  "f_post" = f_post,
  "stdv" = stdv,
  "marginal" = marginal
))
}


#作图,包括拟合前的折线图，拟合后的图线，和不确定性度量
# Gaussian_plot <- function(X, y, Xtest, res){
#   mu = res$mu
# plot(X, y, type='l')
# lines(Xtest, mu,col='red')
# }

#' Plot using result from Gaussian process
#'
#' @param Xtest input test matrix
#' @param Xtrain input training matrix X
#' @param ytrain input training vector y
#'
#' @return plot func, no return variable
#' @import ggplot2
#' @export
#'
#' @examples
#'Xtest = seq(-5, 5, by=10/49)
#'theta = list(1, 10)
#'Xtrain = c(-4, -3, -2, -1, 1)
#'y = sin(Xtrain)
#'Gaussian_plot_pos(Xtest, Xtrain, y)
Gaussian_plot_pos <- function(Xtest,Xtrain, ytrain){
  res = Gaussian_reg_fit_pre(Xtrain, ytrain, Xtest, theta = list(1, 10))
  mu = res$mu
  stdv = res$stdv
  Xtrain = c(Xtrain, rep(Xtrain[1], length(Xtest)-length(Xtrain)))
  ytrain = c(ytrain, rep(ytrain[1], length(Xtest)-length(ytrain)))
data1 = data.frame('test'=Xtest, 'post' = mu)
p <- ggplot(data = data1, mapping = aes(x=test, y=post))
p +ylim(-3, 3)+ geom_line(size=0.8) + labs(x="X_test", y='post',title='samples and uncertainty measure from GPR ') +
  geom_ribbon(aes(ymin = mu-2*stdv, ymax=mu+2*stdv),color='grey',fill='grey',alpha=0.2)+ geom_point(aes(Xtrain, ytrain),shape='square', size=4)
# points(Xtrain, ytrain,col='red',pch=16)
}











#' main function running for spatiotemporal prediction
#'
#' @param Xtrain_s spatial train matrix
#' @param Xtrain_t temporal train matrix
#' @param ytrain ytrain matrix
#' @param Xtest_s spatial test matrix
#' @param Xtest_t temporal test matrix
#' @param kernel kernel, default: Gaussian kernel
#' @param parameters param for kernel,spatial and temporal
#'
#' @return prediction list containing mean and variance
#' @export
#'
run <- function(Xtrain_s, Xtrain_t, ytrain, Xtest_s, Xtest_t, kernel, parameters){
  #using the kernel width
  #ytrain should be the response matrix, not vec
  s_s = parameters$kw_s

K_spatial_tra_tra <- kernel(Xtrain_s, Xtrain_s, s_s)
  K_spatial_test_tra <- kernel(Xtest_s, Xtrain_s, s_s)
  K_spatial_test_test <- kernel(Xtest_s, Xtest_s, s_s)

  s_t = parameters$kw_t
  K_temporal_tra_tra <- kernel(Xtrain_t, Xtrain_t, s_t)
  K_temporal_test_tra <- kernel(Xtest_t, Xtrain_t, s_t)
  K_temporal_test_test <- kernel(Xtest_t, Xtest_t, s_t)

  Y_tra <- ytrain



  #training
  state <- spatiotemporal_gpr_train(K_spatial_tra_tra, K_temporal_tra_tra, Y_tra, parameters)

#prediction
prediction <- spatiotemporal_gpr_test(K_spatial_test_tra, K_spatial_test_test, K_temporal_test_tra, K_temporal_test_test, state)


# print(prediction$Y_mean)
return (prediction)}

#' training function for spatiotemporal gpr
#'
#' @param K_1_tra_tra kernel for train
#' @param K_2_tra_tra kernel for train
#' @param Y_tra train y
#' @param parameters params
#'
#' @return state information
spatiotemporal_gpr_train <- function(K_1_tra_tra, K_2_tra_tra, Y_tra, parameters) {
  svd_1 <- svd(K_1_tra_tra)
  svd_2 <- svd(K_2_tra_tra)

  state <- list(U_1 = svd_1$u, d_1 = svd_1$d, U_2 = svd_2$u, d_2 = svd_2$d, Y_tra = Y_tra, parameters = parameters)
}
#' fit function for spatiotemporal gpr
#'
#' @param K_1_test_tra kernel for train_test
#' @param K_1_test_test kernel for test
#' @param K_2_test_tra kernel for train_test
#' @param K_2_test_test kernel for test
#' @param state state return from train func
#'
#' @return prediction
spatiotemporal_gpr_test <- function(K_1_test_tra, K_1_test_test, K_2_test_tra, K_2_test_test, state) {
  N_1_tra <- ncol(K_1_test_tra)
  N_1_test <- nrow(K_1_test_tra)
  N_2_tra <- ncol(K_2_test_tra)
  N_2_test <- nrow(K_2_test_tra)

  diag_term <- (matrix(matrix(state$d_1, N_2_tra, N_1_tra, byrow = TRUE), N_1_tra * N_2_tra, 1) * matrix(matrix(state$d_2, N_2_tra, N_1_tra, byrow = FALSE), N_1_tra * N_2_tra, 1) + state$parameters$sigma^2)^-1

  Y_mean <- t(matrix(matrix((K_2_test_tra %*% state$U_2) %*%
                              matrix(diag_term * matrix(t(state$U_2) %*% t(state$Y_tra) %*% state$U_1, nrow = N_1_tra * N_2_tra, ncol = 1), nrow = N_2_tra, ncol = N_1_tra)
                            %*% t(K_1_test_tra %*% state$U_1), nrow = N_1_test * N_2_test, ncol = 1), N_2_test, N_1_test))

  Y_variance <- matrix(0, nrow = nrow(Y_mean), ncol = ncol(Y_mean))
  for (row in 1:nrow(Y_variance)) {
    for (column in 1:ncol(Y_variance)) {
      temp <- kronecker(K_1_test_tra[row,] %*% state$U_1, K_2_test_tra[column,] %*% state$U_2)
      Y_variance[row, column] <- K_1_test_test[row, row] * K_2_test_test[column, column] - temp %*% (diag_term * t(temp))
    }
  }

  prediction <- list(Y_mean = Y_mean, Y_variance = Y_variance)
}

# X_train_s = matrix(seq(20), nrow=10, ncol=2)
# X_train_t = matrix(rnorm(24), nrow=8, ncol=3)
# ytrain = matrix(1, nrow=10, ncol=8)
# X_test_s = matrix(c(1, 2), nrow=1, ncol=2)
# X_test_t = matrix(rnorm(3), nrow=1, ncol=3)
#
# run(X_train_s, X_train_t, ytrain, X_test_s, X_test_t, Gaussian_kernel, theta=5)
