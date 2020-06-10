\name{aucC-package}
\alias{aucC-package}
\alias{aucC}
\docType{package}
\title{
  Quick AUC Computation
}
\description{
  Computes the AUC using a simple for-loop implemented in C++. As of 2020, this implementation was 10x faster than pROC::auc.
}
\details{
  Use the aucC_wrap function. See the docs for aucC_wrap.
}
\author{
Jeremiah Jones.

Maintainer: Jeremiah Jones <jeremiah_jones@urmc.rochester.edu>
}

\keyword{ package }
\seealso{
  \code{\link{aucC_wrap}}
}
\examples{
  \dontrun{
     y=rbinom(100, size=1, prob=1/2)
     label = factor(y, levels=c(0,1), labels=c("No", "Yes"))
     scores = runif(100) # higher scores predict label == "Yes"
     aucC::aucC_wrap(label, scores)
  }
}
