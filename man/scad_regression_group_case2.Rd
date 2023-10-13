\name{scad_regression_group_case2}
\alias{scad_regression_group_case2}
\title{scad_regression_group_case2}
\description{
Fit Group SCAD for interval-censored data}
\usage{
scad_regression_group_case2(initial_val,n,J,p,data,m,lambda,epsilon,maxiter)
}
\arguments{
\item{initial_val}{a vector of initial \eqn{\beta}.}
\item{n}{the total number of subjects in m groups.}
\item{J}{number of .}
\item{p}{the total number of coefficients.}
\item{data}{a list of observation time intervals and covariate matrix.}
\item{m}{the number of studies.}
\item{lambda}{a list of candicate penalty parameter \eqn{\lambda}.}
\item{epsilon}{the tolerance parameter to stop the algorithm.}
\item{maxiter}{the number of max round of iteration.}
}
\details{
Group SCAD \cr
\deqn{Q(\bm{\beta}) = L(\bm{\beta}) + n\sum_{k=1}^Kp_\lambda(\lVert \beta_k^*\rVert_2)}
}

\references{

}



