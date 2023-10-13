\name{generate_data}
\alias{generate_data}
\title{generate_data}
\description{
Generate simulation survival data for meta-analysis.
}
\usage{
generate_data(N ,p , correlate , rightCensor ,intervalCensor, M =3, q ,seed)
}
\arguments{
\item{M}{the number of studies or groups in the simulated data}
\item{N}{
the number of subjects (patients) in each group.
}
\item{p}{the number of covariates for the simulation data}
\item{correlate}{logical. If \code{TRUE}, the covariance matrix of generated covariates X has a autocorrelated structure; If \code{FALSE}, the covariance is block-structured. }
\item{rightCensor}{logical. If \code{TRUE}, the event time is generated right-censoring scheme.}
\item{IntervalCensor}{logical. If \code{TRUE}, the event time is generated from  interval-censoring scheme.}
\item{q}{the proportion of true coefficients to be non-zero }
\item{seed}{a random seed}
}
\examples{
seed = sample(1,1:10e8)
generate_data(N=c(500,500,500) ,p =10, correlate = TRUE , rightCensor ,intervalCensor, M =3, q=0.2 ,seed)
}
