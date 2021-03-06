% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/cv.enpls.R
\name{cv.enpls}
\alias{cv.enpls}
\title{Cross Validation for Ensemble Partial Least Squares Regression}
\usage{
cv.enpls(x, y, nfolds = 5L, verbose = TRUE, ...)
}
\arguments{
\item{x}{predictor matrix}

\item{y}{response vector}

\item{nfolds}{number of folds - default is \code{5}.}

\item{verbose}{shall we print the cross validation process}

\item{...}{other arguments that can be passed to \code{enpls.en}}
}
\value{
A list containing four components:
\itemize{
\item \code{ypred} - a matrix containing two columns: real y and predicted y
\item \code{residual} - cross validation result (y.pred - y.real)
\item \code{RMSE} - RMSE
\item \code{R2} - R2
}
}
\description{
This function performs k-fold cross validation for
ensemble partial least squares regression.
}
\details{
This function performs k-fold cross validation for
ensemble partial least squares regression.
}
\examples{
data(alkanes)
x = alkanes$x
y = alkanes$y

set.seed(42)
cv.enpls.fit = cv.enpls(x, y, MCtimes = 10)
print(cv.enpls.fit)
plot(cv.enpls.fit)
}
\author{
Min-feng Zhu <\email{wind2zhu@163.com}>,
        Nan Xiao <\email{road2stat@gmail.com}>
}
\references{
Dongsheng Cao, Yizeng Liang, Qingsong Xu, Yifeng Yun, and Hongdong Li.
"Toward better QSAR/QSPR modeling: simultaneous outlier detection and
variable selection using distribution of model features."
\emph{Journal of computer-aided molecular design} 25, no. 1 (2011): 67--80.
}
\seealso{
See \code{\link{enpls.en}} for ensemble PLS regression.
}

