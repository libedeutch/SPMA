\name{dist2beta0}
\alias{dist2beta0}
\title{dist2beta0}
\description{
Calculate the distance between true beta and estimated beta
}
\usage{
dist2beta0(beta0,betahat,p,M)
}
\arguments{
\item{beta0}{true coefficients in generated data.}
\item{betahat}{estimated coefficients.}
\item{p}{the number of covariates in the simulation data.}
\item{M}{the number of generated studies}
}
\examples{
beta0 = c(1:30)
betahat = rnorm(30,2,5)
dist2beta(beta0,betahat,10,3)
}

