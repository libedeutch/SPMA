\name{bayes_bic_spma}
\alias{bayes_bic_spma}
\title{bayes_bic_spma}
\description{
Calculate modified BIC for Lasso, SCAD, MCP, Group Lasso, Group SCAD, Group MCP
}
\usage{
bayes_bic_spma(set.n,thetahat,sample1,sample2,sample3,p,M,intervalCensor)
}
\arguments{
\item{set.n}{sample size vector for three studies.}
\item{thetahat}{estimated coefficients.}
\item{sample1, sample2, sample3}{Covariate matrix of three studies.}
\item{p}{the number of covariates.}
\item{M}{he number of studies, default 3.}
\item{intervalCensor}{logical. If \code{TRUE}, evaluate the log likelihood of interval-censored data.}
}
\details{
For interval-censored data, \cr
\deqn{\text{MBIC} =  -2 \ell(\beta,\Lambda_0) +  {\color{black}{ \log n/n }} \sum_{m=1}^M \sum_{j=1}^pI(\hat{\beta}_j^m \neq 0 )} \cr
For right-censored data and data without censoring, \cr
\deqn{\text{MBIC} = -2\sum_{m=1}^M \sum_{i=1}^{n_m} \delta_i^m
 \left[ \bm{x}_{i}^{m\top}\bm{\hat{\beta}}^m - \log \left(\sum_{j: T_j\geq T_i} ^{n_m}\exp(\bm{x}_{j}^{m\top}\bm{\hat{\beta}}^m)\right)\right] +  {\color{black}{ \log n/n }} \sum_{m=1}^M \sum_{j=1}^pI(\hat{\beta}_j^m \neq 0 )}
}
\examples{
seed = sample(1:10e8,1)
simulData = generate_data(N=c(500,500,500), p =10, correlate = TRUE ,
rightCensor = FALSE ,intervalCensor = FALSE, M =3, q=0.2 ,seed)
modl = trainSPMA(simulData)
# Unpenalized estimator
thetaF = modl$unpel
#take first model as example
thetahat = modl$spma[[1]]$alphaM[,,10]
bayes_bic_spma(thetaF,thetahat,
modl$covF1,modl$covF2,modl$covF3,p = 10, N =c(500,500,500) )
}

