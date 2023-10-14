\name{penalty.factor.f}
\alias{penalty.factor.f}
\title{penalty.factor.f}
\description{
Calculate penalty function of lasso and scad }
\usage{
penalty.factor.f(beta0 = NULL, lasso.lam = NULL, ncov,
npieces, penalty.function, nopenalty.index = NULL)
}
\arguments{
\item{beta0}{initialized coefficients.}
\item{lasso.lam}{penalty parameter \eqn{\lambda}.}
\item{ncov}{the number of covariates.}
\item{npieces}{the number of pieces in the piecewise constant exponential model.}
\item{panelty.function}{penalty, \code{scad} or \code{lasso}.}
}
\details{
Lasso: \cr
SCAD: \cr
}

\references{
  Wu, Y. & Cook, R.J. (2015) \emph{Penalized regression for interval‐censored times of disease progression: Selection of HLA markers in psoriatic arthritis}. Biometrics, 71(3): 782-791.
}




