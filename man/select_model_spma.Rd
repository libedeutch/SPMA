\name{select_model_spma}
\alias{select_model_spma}
\title{select_model_spma}
\description{
Select the best model for SPMA based on modified BIC
}
\usage{
select_model_spma(trainSPMA.object)
}
\arguments{
\item{trainSPMA.object}{a trainSPMA object.}
}
\details{
\deqn{\text{MBIC} =\sum_{m=1}^M n_m (\bm{\hat{\beta}}^m - \bm{\hat{\beta}}^m_F)^\top \hat{\bm{I}}_F^m (\bm{\hat{\beta}}^m - \bm{\hat{\beta}}^m_F) +  {\color{black}{ \log n/n }}\sum_{m=1}^M \sum_{j=1}^pI(\hat{\beta}_j^m \neq 0 ) } \cr
}
\examples{
seed = sample(1:10e8,1)
simulData = generate_data(N=c(500,500,500) ,p =10, correlate = TRUE , rightCensor = TRUE,
intervalCensor = FALSE, M =3, q=0.2 ,seed)
modl = trainSPMA(simulData)
spma.bic = select_model_spma(modl)
}

