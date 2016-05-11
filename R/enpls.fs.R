#' Ensemble Partial Least Squares for Feature Selection
#'
#' This function performs feature selection with ensemble partial least squares.
#'
#' This function performs feature selection with ensemble partial least squares.
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
#' @return A list containing two components:
#' \itemize{
#' \item \code{variable.importance} - a vector of variable importance
#' \item \code{coefficient.matrix} - original coefficient matrix
#' }
#'
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>,
#'         Nan Xiao <\email{road2stat@@gmail.com}>
#'         
#' @seealso See \code{\link{enpls.od}} for outlier detection with ensemble PLS. 
#' See \code{\link{enpls.en}} for ensemble PLS regression.
#' \code{\link{enpls.ad}} for applicability domain with ensemble PLS.
#'
#' @export enpls.fs
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
#' varimp = enpls.fs(x, y, MCtimes = 100)
#' print(varimp)
#' plot(varimp)

enpls.fs = function(x, y, 
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

    coeflist = vector('list', MCtimes)
    for (i in 1L:MCtimes) {
      plsdf.x = x[samp.idx[[i]], ]
      plsdf.y = y[samp.idx[[i]]]
      coeflist[[i]] = suppressWarnings(enpls.fs.core(plsdf.x, plsdf.y, maxcomp))
    }

  } else {

    registerDoParallel(parallel)
    coeflist = foreach(i = 1L:MCtimes) %dopar% {
      x = x[samp.idx[[i]], ]
      y = y[samp.idx[[i]]]
      enpls.fs.core(x, y, maxcomp)
    }

  }

  coefmat = do.call(rbind, coeflist)
   
  varimp = abs(colMeans(coefmat))/apply(coefmat, 2L, sd)
  
  varimp[which(is.na(varimp))] = 0

  object = list('variable.importance' = varimp, 
                'coefficient.matrix'  = coefmat)
  class(object) = 'enpls.fs'
  return(object)

}

#' core function for enpls.fs
#'
#' select the best ncomp with cross-validation and
#' use it to fit the complete training set again.
#' scale = FALSE
#'
#' @return fitted coefficients
#' 
#' @keywords internal

enpls.fs.core = function(plsdf.x, plsdf.y, maxcomp = NULL) {
  
  plsr.coef = matrix(0, 1, dim(plsdf.x)[2])
  colnames(plsr.coef) = colnames(plsdf.x)
  num = apply(plsdf.x, 2L, function(x) sum(x == 0))
  cv.ratio = num/dim(plsdf.x)[1]
  weight = which(cv.ratio > 0.8)
  if(length(weight) > 0) {
    plsdf.x = plsdf.x[, -weight]
    maxcomp = ncol(plsdf.x)
  }
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
  
  
  coe = drop(coef(plsr.fit))
  for(i in 1:dim(plsr.coef)[2]) {
    idx = which(colnames(plsr.coef)[i] == names(coe))
    if(length(idx)!=0) plsr.coef[1,i] = coe[idx]
  }
  
  plsr.coef = drop(plsr.coef)
  return(plsr.coef)

}
