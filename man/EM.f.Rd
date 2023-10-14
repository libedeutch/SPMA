\name{EM.f}
\alias{EM.f}
\title{EM.f}
\description{
Fit Cox PH with lasso and scad penalty for interval-censored data}
\usage{
EM.f(indata, lam0, beta0, lasso.lam, ncov, npieces, cutpoints, penalty.function,
penalty.factor = NULL, nopenalty.index = NULL, thresh = 1e-6, maxit = 200)
}
\arguments{
\item{indata}{a data matrix with id, event time endpoints, covariate matrix, survival status.}
\item{lam0}{bashline hazard initialization.}
\item{beta0}{initialized coefficients.}
\item{lasso.lam}{penalty parameter \eqn{\lambda}.}
\item{ncov}{the number of covariates.}
\item{npieces}{the number of pieces in the piecewise constant exponential model.}
\item{cutpoints}{intervals of each piece.}
\item{penalty.function}{penalty, \code{scad} or \code{lasso}.}
}
\details{
Piecewise constant exponential model is employed to approximate the underlying bashline hazard function.
}

\references{
  Wu, Y. & Cook, R.J. (2015) \emph{Penalized regression for interval‐censored times of disease progression: Selection of HLA markers in psoriatic arthritis}. Biometrics, 71(3): 782-791.
}




