#' Ensemble Partial Least Squares Regression
#'
#' This function performs ensemble partial least squares regression.
#'
#' This function performs ensemble partial least squares regression.
#'
#' @param x predictor matrix
#' @param y response vector
#' @param maxcomp Maximum number of components included within the models, 
#' if not specified, default is the variable (column) numbers in x.
#' @param MCtimes times of Monte-Carlo
#' @param method \code{"mc"} or \code{"bootstrap"}. Default is \code{"mc"}.
#' @param ratio sample ratio used when \code{method = "mc"}
#' @param parallel Integer. Number of parallel processes to use. 
#' Default is \code{1}, which means run serially.
#'
#' @return A list containing all PLS model objects.
#'
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>,
#'         Nan Xiao <\email{road2stat@@gmail.com}>
#'         
#'
#' @seealso See \code{\link{enpls.fs}} for feature selection with ensemble PLS. 
#' See \code{\link{enpls.od}} for outlier detection with ensemble PLS. See 
#' \code{\link{enpls.ad}} for applicability domain with ensemble PLS.
#'
#' @export enpls.en
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#'
#' @references
#' Dongsheng Cao, Yizeng Liang, Qingsong Xu, Yifeng Yun, and Hongdong Li. 
#' "Toward better QSAR/QSPR modeling: simultaneous outlier detection and 
#' variable selection using distribution of model features." 
#' \emph{Journal of computer-aided molecular design} 25, no. 1 (2011): 67--80.
#'
#' @examples
#' data(alkanes)
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' enpls.fit = enpls.en(x, y, MCtimes = 10)
#' print(enpls.fit)
#' predict(enpls.fit, newx = x)

enpls.en = function(x, y, 
                    maxcomp = NULL, 
                    MCtimes = 500L, 
                    method = c('mc', 'bootstrap'), ratio = 0.8, 
                    parallel = 1L) {

  if (is.null(maxcomp)) maxcomp = ncol(x)

  method = match.arg(method)

  x.row = nrow(x)
  samp.idx = vector('list', MCtimes)

  if (method == 'mc') {
    for (i in 1L:MCtimes) samp.idx[[i]] = sample(1L:x.row, floor(x.row * ratio))
  }

  if (method == 'bootstrap') {
    for (i in 1L:MCtimes) samp.idx[[i]] = sample(1L:x.row, x.row, replace = TRUE)
  }

  if (parallel < 1.5) {

    modellist = vector('list', MCtimes)
    for (i in 1L:MCtimes) {
      plsdf.x = x[samp.idx[[i]], ]
      plsdf.y = y[samp.idx[[i]]]
      modellist[[i]] = suppressWarnings(enpls.en.core(plsdf.x, plsdf.y, maxcomp))
    }

  } else {

    registerDoParallel(parallel)
    modellist = foreach(i = 1L:MCtimes) %dopar% {
      x = x[samp.idx[[i]], ]
      y = y[samp.idx[[i]]]
      enpls.en.core(x, y, maxcomp)
    }

  }

  class(modellist) = 'enpls.en'
  return(modellist)

}

#' core function for enpls.en
#'
#' select the best ncomp with cross-validation and
#' use it to fit the complete training set again.
#' scale = TRUE
#'
#' @return the coefficients
#' 
#' @keywords internal

enpls.en.core = function(plsdf.x, plsdf.y, maxcomp) {
 
    num = apply(plsdf.x, 2L, function(x) sum(x == 0))
    cv.ratio = num/dim(plsdf.x)[1]
    weight = which(cv.ratio > 0.7)
    if(length(weight) > 0) {
      plsdf.x = plsdf.x[, -weight]
    }
    if(maxcomp > ncol(plsdf.x)) maxcomp = ncol(plsdf.x)
    
 
  plsr.cvfit = plsr(plsdf.y ~ ., data = data.frame(plsdf.x, plsdf.y),  
                    ncomp  = maxcomp, 
                    scale  = TRUE, 
                    method = 'simpls', 
                    validation = 'CV', segments = 5L)

  # choose best component number using adjusted CV
  cv.bestcomp = which.min(RMSEP(plsr.cvfit)[['val']][2L, 1L, -1L])

  plsr.fit = plsr(plsdf.y ~ ., data = data.frame(plsdf.x, plsdf.y), 
                  ncomp  = cv.bestcomp, 
                  scale  = TRUE, 
                  method = 'simpls', 
                  validation = 'none')

  enpls.core.fit = list(plsr.fit, cv.bestcomp)  # save cv.bestcomp for predict.enpls
  return(enpls.core.fit)

}
