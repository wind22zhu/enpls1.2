% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/enpls.od.R
\name{enpls.od}
\alias{enpls.od}
\title{Ensemble Partial Least Squares for Outlier Detection}
\usage{
enpls.od(x, y, maxcomp = NULL, MCtimes = 500L, method = c("mc",
  "bootstrap"), ratio = 0.8, parallel = 1L)
}
\arguments{
\item{x}{predictor matrix}

\item{y}{response vector}

\item{maxcomp}{Maximum number of components included within the models,
if not specified, default is the variable (column) numbers in x.}

\item{MCtimes}{times of Monte-Carlo}

\item{method}{\code{"mc"} or \code{"bootstrap"}. Default is \code{"mc"}.}

\item{ratio}{sample ratio used when \code{method = "mc"}}

\item{parallel}{Integer. Number of parallel processes to use.
Default is \code{1}, which means run serially.}
}
\value{
A list containing four components:
\itemize{
\item \code{error.mean} - error mean for all samples (absolute value)
\item \code{error.median} - error median for all samples
\item \code{error.sd} - error sd for all samples
\item \code{predict.error.matrix} - the original prediction error matrix
}
}
\description{
This function performs outlier detection with ensemble partial least squares.
}
\details{
This function performs outlier detection with ensemble partial least squares.
}
\examples{
data(alkanes)
x = alkanes$x
y = alkanes$y

set.seed(42)
od = enpls.od(x, y, MCtimes = 100)
print(od)
plot(od)
plot(od, criterion = 'sd')
}
\author{
Min-feng Zhu <\email{wind2zhu@163.com}>,
        Nan Xiao <\email{road2stat@gmail.com}>
}
\references{
DongSheng Cao, Yizeng Liang, Qingsong Xu, Hongdong Li, and Xian Chen.
"A new strategy of outlier detection for QSAR/QSPR."
\emph{Journal of computational chemistry} 31, no. 3 (2010): 592--602.

Dongsheng Cao, Yizeng Liang, Qingsong Xu, Yifeng Yun, and Hongdong Li.
"Toward better QSAR/QSPR modeling: simultaneous outlier detection and
variable selection using distribution of model features."
\emph{Journal of computer-aided molecular design} 25, no. 1 (2011): 67--80.
}
\seealso{
See \code{\link{enpls.fs}} for feature selection with ensemble PLS.
See \code{\link{enpls.en}} for ensemble PLS regression.
See \code{\link{enpls.ad}} for applicability domain with ensemble PLS
}

