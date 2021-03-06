% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/print.R
\name{print.enpls.ad}
\alias{print.enpls.ad}
\title{Print enpls.ad Object}
\usage{
\method{print}{enpls.ad}(x, ...)
}
\arguments{
\item{x}{An object of class \code{enpls.ad}.}

\item{...}{Other parameters to be passed on to \code{print}.}
}
\description{
This function prints enpls.ad object.
}
\details{
This function prints enpls.ad object.
}
\examples{
data(logS)
x = logS$x
y = logS$y
x.test = logS$x.test1
y.test = logS$y.test1

set.seed(42)
ad = enpls.ad(x, y, x.test, y.test, MCtimes = 10)
print(ad)
}
\author{
Min-feng Zhu <\email{wind2zhu@163.com}>,
         Nan Xiao <\email{road2stat@gmail.com}>
}
\seealso{
See \code{\link{enpls.ad}} for applicability domain with ensemble PLS.
}

