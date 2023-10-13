\name{createdata.f}
\alias{createdata.f}
\title{createdata.f}
\description{
Creat pseudo data set }
\usage{
createdata.f(indata, lam, beta, ncov, cutpoints)
}
\arguments{
\item{indata}{a data matrix with id, event time endpoints, covariate matrix, survival status.}
\item{beta}{estimated coefficients.}
\item{lam}{penalty parameter \eqn{\lambda}.}
\item{ncov}{the number of covariates.}
\item{cutpoints}{intervals of each piece.}
}
\details{

}

\references{
  Wu, Y. & Cook, R.J. (2015) \emph{Penalized regression for interval‐censored times of disease progression: Selection of HLA markers in psoriatic arthritis}. Biometrics, 71(3): 782-791.
}




