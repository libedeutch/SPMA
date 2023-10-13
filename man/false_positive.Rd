\name{flase_positive}
\alias{false_positive}
\title{false_positive}
\description{
Calculate false positive rate from confusion matrix
}
\usage{
false_positive(confuse)
}
\arguments{
\item{confuse}{a confuse matrix for two classes.}
}
\details{
\tabular{cccc}{
 \tab \tab  True \cr
 \tab \tab P \tab N \cr
 Prediction \tab P \tab TP \tab FP\cr
 \tab N  \tab FN \tab TN\cr
}
\deqn{\text{FPR} = \frac{FP}{FP+TN}}
}


