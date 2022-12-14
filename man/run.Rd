% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gaussian_process.R
\name{run}
\alias{run}
\title{run(Gaussianst)}
\usage{
run(
  Xtrain_s,
  Xtrain_t,
  ytrain,
  Xtest_s,
  Xtest_t,
  kernel = Gaussian_kernel,
  parameters
)
}
\arguments{
\item{Xtrain_s}{spatial train matrix}

\item{Xtrain_t}{temporal train matrix}

\item{ytrain}{ytrain matrix}

\item{Xtest_s}{spatial test matrix}

\item{Xtest_t}{temporal test matrix}

\item{kernel}{kernel, \strong{default: Gaussian kernel}}

\item{parameters}{params for spatial kernel, temporal kernel and the sigma. \strong{Suggest:}The kernel width parameters are chosen as the means of pairwise Euclidean distances between training instances for all kernels.}
}
\value{
prediction list containing posterior mean and variance
}
\description{
Main function running for spatiotemporal prediction
}
\details{
Implementation for spatiotemporal prediction, including  uncertainty measuring
At the same time, the package can be used to make predictions under three different scenatios:
\enumerate{\item temporal prediction: predicting case counts for future time periods, leading to predicting outcome for each location in the future
\item spatial prediction: predicting case counts for unseen locations, leading to predicting outcome spread within the same time frame in other locations
\item  spatiotemporal prediction: predicting case counts for unseen location and future time period pairs, leading to predicting outcome spread to new locations in the future}
}
\examples{
X_train_s = matrix(seq(20), nrow=10, ncol=2)
X_train_t = matrix(rnorm(24), nrow=8, ncol=3)
ytrain = matrix(1, nrow=10, ncol=8)
X_test_s = matrix(c(1, 2), nrow=1, ncol=2)
X_test_t = matrix(rnorm(3), nrow=1, ncol=3)
parameters <- list()
parameters$sigma <- 1e-3
parameters$kw_s <- 1
parameters$kw_t <- 1
run(X_train_s, X_train_t, ytrain, X_test_s, X_test_t, Gaussian_kernel, parameters)

}
\references{
Ak, Çiğdem, et al. "Spatiotemporal prediction of infectious diseases using structured Gaussian processes with application to Crimean–Congo hemorrhagic fever." PLoS neglected tropical diseases 12.8 (2018): e0006737.\url{https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0006737}
}
\seealso{
\code{\link{spatiotemporal_gpr_train}}

\code{\link{spatiotemporal_gpr_test}}
}
