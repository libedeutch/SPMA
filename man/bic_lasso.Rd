\name{bic_lasso}
\alias{bic_lasso}
\title{bic_lasso}
\description{
Calculate modified BIC for Lasso
}
\usage{
bic_lasso(ni,thetahat,sample,p)
}
\arguments{
\item{ni}{sample size of the study.}
\item{thetahat}{estimated coefficients.}
\item{sample}{covariate matrix of the study.}
\item{p}{number of covariates.}
}
\examples{
seed = sample(1:10e8,1)
N=c(500,500,500)
simulData = generate_data(N, p =10, correlate = TRUE,
rightCensor=TRUE ,intervalCensor=FALSE, M =3, q=0.2 ,seed)
modl = trainSPMA(simulData)
# take lasso for study one as an example
thetahat = modl$lasso1$beta[,1]
bic_lasso(N[1],thetahat, simulData$sample1,p =10)
}

