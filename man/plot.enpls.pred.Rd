% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.enpls.pred}
\alias{plot.enpls.pred}
\title{Plot enpls.pred Object}
\usage{
\method{plot}{enpls.pred}(x, y.test, main = NULL, ...)
}
\arguments{
\item{x}{An object of class \code{enpls.pred}.}

\item{y.test}{Response vector for test.}

\item{main}{plot title}

\item{...}{Other graphical parameters to be passed on to \code{ggplot}.}
}
\description{
This function plots enpls.pred object.
}
\details{
This function plots enpls.pred object.
}
\examples{
data(logS)
x = logS$x
y = logS$y
x.test = logS$x.test1
y.test = logS$y.test1

set.seed(42)
enpls.fit = enpls.en(x, y, MCtimes = 10)
y.pred = predict(enpls.fit, newx = x.test)
plot(y.pred, y.test)
}
\author{
Min-feng Zhu <\email{wind2zhu@163.com}>,
        Nan Xiao <\email{road2stat@gmail.com}>
}
\seealso{
See \code{\link{predict}} for  make predictions on new data by fitted enpls.en object.
}

