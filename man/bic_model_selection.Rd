\name{bic_model_selection}
\alias{bic_model_selection}
\title{bic_model_selection}
\description{
Selection best model of a list of candicate models based on modified BIC
}
\usage{
bic_model_selection(trainSPMA.object)
}
\arguments{
\item{trainSPMA.object}{a trainSPMA object.}
}
\details{
For interval-censored data, modified BIC are calculated after training of each model; for right-censored data and data without censoring, modified BIC are calculated for each model corresponding to a penalty parameter \eqn{\lambda}.
}
\examples{
seed = sample(1:10e8,1)
N=c(500,500,500)
simulData = generate_data(N,p =10, correlate = TRUE,
rightCensor = TRUE,intervalCensor = FALSE, M =3, q=0.2 ,seed)
modl = trainSPMA(simulData)
bic_model_selection(modl)
}
