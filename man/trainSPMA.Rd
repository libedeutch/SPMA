\name{trainSPMA}
\alias{trainSPMA}
\title{trainSPMA}
\description{
Implement unpenalized Cox PH model, Cox PH model with penalities Lasso, SCAD, MCP, Group Lasso, Group SCAD, Group MCP, SPMA.
}
\usage{
trainSPMA(simulData, iter,m)
}
\arguments{
\item{simulData}{a generae_data object.}
\item{iter}{
the max iteration for ADMM algorithm of SPMA.
}
\item{m}{the degree of Berstein polynomial (default 5) for Baseline hazard model approximation for interval-censored data.}
}
\examples{
seed = sample(1,1:10e8)
simulData = (N=c(500,500,500) ,p =10, correlate = TRUE , rightCensor ,intervalCensor, M =3, q=0.2 ,seed)
trainSPMA(simulData)
}

