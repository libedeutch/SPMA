\name{printResult}
\alias{printResult}
\title{printResult}
\description{
Generate the averaged error metrics across replciates
}
\usage{
printResult(rre, intervalCensor)
}
\arguments{
\item{rre}{a list of trainSPMA object.}
\item{intervalCensor}{logical. If \code{TRUE}, error metrics are calculuted for interval-censored data.}

}
\examples{
rre = list()
seed = runif(10,0,10e8)
for(i in 1:10){
  cc = generate_data(p=8,intervalCensor = FALSE,rightCensor = TRUE,
  seed = seed[i], N = c(1000,1000,1000))
  dd = trainSPMA(cc,m = 5)
  best.idx = bic_model_selection(dd)
  spma.bic = select_model_spma(dd)
  rmse0 = rmse(trainSPMA.object=dd, best.idx = best.idx, spma.bic = spma.bic, Data = cc )
  rre[[i]] = rmse0
}
printResult(rre, intervalCensor = TRUE)
}

