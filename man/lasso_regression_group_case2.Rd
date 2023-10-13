\name{lasso_regression_group_case2}
\alias{lasso_regression_group_case2}
\title{lasso_regression_group_case2}
\description{
Fit Group Lasso for interval-censored data}
\usage{
lasso_regression_group_case2(initial_val,n,J,p,data,m,lambda,epsilon,maxiter)
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
Group lasso \cr
\deqn{\frac{1}{2}\lVert\bm{Y} - \sum_{j=1}^pX_j\beta_j \rVert^2 + \lambda \sum_{j=1}^p (\beta_j^\top K_j \beta_j)^{1/2}}
}

\references{
  Yuan, M. & Lin, Y. (2006) \emph{Model selection and estimation in regression with grouped variables}. Journal of the Royal Statistical Society Series B: Statistical Methodology, 68(1): 49-67}



