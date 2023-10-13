\name{l_n_case2}
\alias{l_n_case2}
\title{l_n_case2}
\description{
Calculate loglikelihood of interval-censored data
}
\usage{
l_n_case2(phi,data,beta,m,n)
}
\arguments{
\item{phi}{a term bounded by a positive number M.}
\item{data}{a list of observation time intervals and covariate matrix.}
\item{beta}{a vector of coefficients.}
\item{m}{degree of freedom of the Berstein polynomials.}
\item{n}{the number of subjects.}
}
\details{
The likelihood of interval-censored data: \cr
\deqn{L_n(\beta, \Lambda_0) = \prod_{i=1}^n S(L_i) - S(R_i) = \prod_{i=1}^n\left[ \exp\left(-\Lambda_0(L_i)e^{\beta^\top X}\right) - \exp\left(-\Lambda_0(R_i)e^{\beta^\top X}\right)\right]}
The loglikelihood of interval-censored data: \cr
\deqn{\ell(\beta,\Lambda_0) = \sum_{i=1}^n\log\left( \exp\left(-\Lambda_0(L_i)e^{\beta^\top X}\right) - \exp\left(-\Lambda_0(R_i)e^{\beta^\top X}\right) \right)}
}

\references{
  Wu, Y., Zhao, H., & Sun, J. (2023) \emph{Group variable selection for the Cox model with interval‐censored failure time data}. Biometrics.
}


