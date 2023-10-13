\name{rmse}
\alias{rmse}
\title{rmse}
\description{
Calculate error metrics for eight methods, including ME, MRME, TP, FP
}
\usage{
rmse(trainSPMA.object, best.idx, spma.bic, Data)
}
\arguments{
\item{trainSPMA.object}{a trainSPMA object.}
\item{best.idx}{index of best model for Lasso, SCAD, MCP, Group Lasso, Group SCAD and Group MCP.}
\item{spma.bic}{the index of the select model for SPMA.}
\item{Data}{a generate_data object.}
}
\examples{
seed = sample(1,1:10e8)
simulData = (N=c(500,500,500) ,p =10, correlate = TRUE , rightCensor ,intervalCensor, M =3, q=0.2 ,seed)
modl = trainSPMA(simulData)
spma.bic = select_model_spma(modl)$idx.spma
best.idx = bic_model_selection(modl)
rmse(modl, best.idx,spma.bic, simulData)
}

