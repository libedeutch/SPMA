\name{logL.f}
\alias{logL.f}
\title{logL.f}
\description{
loglikelihood of interval-censored data with piecewise constant exponential model}
\usage{
logL.f(T, Z, lam, beta, ncov, cutpoints)
}
\arguments{
\item{T}{a matrix of left and right endpoints of interval.}
\item{Z}{covariate matrix.}
\item{lam}{a vector of baseline hazard.}
\item{beta}{estimated coefficients.}
\item{ncov}{the number of covariates.}
\item{cutpoints}{intervals of each piece.}
}
\details{
Let \eqn{b_0 = 0, 0 <b_1<\cdots, b_{K-1} < b_K} be the \eqn{K} break-points, the baseline hazard function is \eqn{h(t;\alpha) = \lambda_k = \exp(\alpha_k)} for \eqn{t\in \mathcal{B}_k = [b_{k-1},b_k), k = 1, \cdots, K}. Thus, \eqn{h(t) = \sum_{k=1}^K\lambda_kI(t\in\mathcal{B}_k)\exp(\beta^\top X)}, \eqn{H(t) = \int_0^t h(u) = \sum_{k=1}^K \lambda_k\phi_k \exp(\beta^\top X), } where \eqn{\phi_k = \max\{0, \min(b_k - b_{k-1}, t_i -b_{k-1})\}}. \cr
The likelihood function is:

\deqn{\bar{F}(t) =  \exp(-H(t))
}
}
\references{
  Wu, Y. & Cook, R.J. (2015) \emph{Penalized regression for interval‐censored times of disease progression: Selection of HLA markers in psoriatic arthritis}. Biometrics, 71(3): 782-791.
}




