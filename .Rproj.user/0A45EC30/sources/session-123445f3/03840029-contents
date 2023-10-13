\name{ddl_n_case2}
\alias{ddl_n_case2}
\title{ddl_n_case2}
\description{
Calculate the second order derivative of loglikelihood of interval-censored data
}
\usage{
ddl_n_case2(phi,data,beta,m,n)
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
The second order derivative of loglikelihood of interval-censored data:
Let \eqn{e^L = \exp\left(-\Lambda_0(L_i)e^{\beta^\top X}\right)}, \eqn{e^R = \exp\left(-\Lambda_0(R_i)e^{\beta^\top X}\right)}, \eqn{L = -\Lambda_0(L_i)e^{\beta^\top X}}, \eqn{R = -\Lambda_0(R_i)e^{\beta^\top X}}. \cr
The second-order derivative of the log-likelihood function is:

\deqn{\ddot{\ell}(\beta, \Lambda_0) = \frac{(e^LLX LX + e^LXLX- e^RRXRX - e^RXRX) (e^L-e^R)  - (e^LLX - e^RRX)^2}{(e^L - e^R)^2}}

}

\references{
  Wu, Y., Zhao, H., & Sun, J. (2023) \emph{Group variable selection for the Cox model with interval‐censored failure time data}. Biometrics.
}


