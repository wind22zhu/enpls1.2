#' Cross Validation for  Partial Least Squares Regression
#' 
#' This function performs k-fold cross validation for 
#' partial least squares regression.
#'
#' This function performs k-fold cross validation for 
#' partial least squares regression.
#'
#' @param x predictor matrix
#' @param y response vector
#' @param random randomization of y
#' @param nfolds number of folds - default is \code{5}.
#' @param maxcomp Maximum number of components included within the models, 
#' if not specified, default is the variable (column) numbers in x.
#' @param use.scale scale 
#' @param verbose shall we print the cross validation process
#' @param ... other arguments that can be passed to \code{enpls.en}
#'
#' @return A list containing four components:
#' \itemize{
#' \item \code{y} - response vector
#' \item \code{ypred} - predicted y
#' \item \code{residual} - cross validation result (y.pred - y.real)
#' \item \code{RMSE} - RMSE
#' \item \code{R2} - R2
#' }
#'
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>,
#'         Nan Xiao <\email{road2stat@@gmail.com}>
#'         
#' @seealso See \code{\link{cv.enpls}} for Cross Validation for ensemble PLS regression.
#'
#' @export cv.pls
#'
#'
#' @examples
#' data(alkanes)
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' cv.enpls.fit = cv.enpls(x, y, MCtimes = 10)
#' print(cv.enpls.fit)
#' plot(cv.enpls.fit)

cv.pls <- function(x, y, random=TRUE,
                   nfolds = 5, maxcomp = NULL, use.scale=TRUE, 
                   verbose = TRUE, ...) {
  
  if (is.null(maxcomp)) maxcomp = ncol(x) 
  mx = dim(x)[1]
  if (random)
  {
    In = sample(1:mx, mx, replace=FALSE)
    x = x[In, ]
    y = y[In]
  }
  ypred = y
  index = rep(1:nfolds, nrow(x))
  ind = index[1:nrow(x)]
  for (k in 1:nfolds) {
    if (verbose) cat('Beginning fold', k, '\n')
    xcal = x[ind != k, ] 
    ycal = y[ind != k]
    xtest = x[ind == k, ]
    ytest = y[ind==k]
    if (use.scale) {
      sd0 = which(apply(xcal, 2L, sd) == 0)
      if (length(sd0) > 0){
        xcal = xcal[, -sd0]
        xtest = xtest[, -sd0]
      }
      xcal = scale(xcal, center = TRUE, scale = TRUE)
      xcen = attributes(xcal)$'scaled:center'
      xsca = attributes(xcal)$'scaled:scale'
      
      xtest = scale(xtest, xcen, xsca)
      ycal = scale(ycal, center=TRUE, scale=F)
      ycen = attributes(ycal)$'scaled:center'
      
    }
    mvrout = plsr(ycal~., ncomp = min(c(maxcomp, dim(xcal))), 
                 data = data.frame(xcal, ycal), scale=FALSE,
                 method="simpls");
    ypred[ind==k] = predict(mvrout, comps = 1:min(c(maxcomp, dim(xcal))), xtest) + ycen  
  }; cat("\n")
  
  RMSE = sqrt(t(y - ypred)%*%(y - ypred)/mx)
  R2 = 1 - t(y - ypred)%*%(y - ypred)/(t(y - mean(y))%*%(y - mean(y))) 
  residual = y - ypred
  list(y = y, ypred = ypred, residual = residual, 
       RMSE = RMSE, R2 = R2)
}