#' Ensemble Partial Least Squares for Applicability Domain
#'
#' This function performs applicability domain with ensemble partial least squares.
#'
#' This function performs applicability domain with ensemble partial least squares.
#'
#' @param x predictor matrix.
#' @param y response vector.
#' @param x.test predictor matrix for test.
#' @param y.test response vector for test.
#' @param maxcomp Maximum number of components included within the models, 
#' if not specified, default is the variable (column) numbers in x.
#' @param MCtimes times of Monte-Carlo.
#' @param method \code{"mc"} or \code{"bootstrap"} or \code{"jagging"}. Default is \code{"mc"}.
#' @param verbose shall we print the MCtimes process.
#' @param ratio sample ratio used when \code{method = "mc"} and \code{method = "jagging"}.
#' @param parallel Integer. Number of parallel processes to use. Default is \code{1}, which means run serially.
#'
#' @return A list containing four components:
#' \itemize{
#' \item \code{STD.cv} - STD value for training set 
#' \item \code{STD.te} - STD value for test set
#' \item \code{error.cv} - absolute prediction error of training set
#' \item \code{error.te} - absolute prediction error of test set
#' }
#'
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>,
#'         Nan Xiao <\email{road2stat@@gmail.com}>
#'         
#' @seealso See \code{\link{enpls.fs}} for feature selection with ensemble PLS. 
#' See \code{\link{enpls.en}} for ensemble PLS regression.
#' See \code{\link{enpls.od}} for Outlier Detection with ensemble PLS
#' 
#'
#' @export enpls.ad
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#'
#' @references
#' Kaneko H, Funatsu K. 
#' "Applicability Domain Based on Ensemble Learning in Classification and Regression Analyses." 
#' \emph{Journal of chemical information and modeling} 54, no. 9 (2014): 2469-2482.
#'
#' @examples
#' data(logS)
#' x = logS$x
#' y = logS$y
#' x.test1 = logS$x.test1
#' y.test1 = logS$y.test1
#' 
#' set.seed(42)
#' ad_test1 = enpls.ad(x, y, x.test1, y.test1, MCtimes = 10)
#' print(ad_test1)
#' plot(ad_test1)
#' 
#' x.test2 = logS$x.test2
#' y.test2 = logS$y.test2
#' 
#' set.seed(42)
#' ad_test2 = enpls.ad(x, y, x.test2, y.test2, MCtimes = 10)
#' print(ad_test2)
#' plot(ad_test2)
#' 

enpls.ad = function(x, y, 
                    x.test = NULL, 
                    y.test = NULL, 
                    maxcomp = NULL, 
                    MCtimes = 500L, 
                    verbose = FALSE, 
                    method = c('mc', 'bootstrap', 'jagging'), ratio = 0.8,  
                    parallel = 1L) {
  
  if (is.null(x.test)) x.test = x
  
  
  method = match.arg(method)
  
  x.row = nrow(x)
  x.col = ncol(x)
  xte.row = nrow(x.test)
  samp.idx = vector('list', MCtimes)
  
  if (method == 'mc') {
    for (i in 1L:MCtimes) {
      samp.idx[[i]] = sample(1L:x.row, floor(x.row * ratio))
    }
  }
  
  if (method == 'bootstrap') {
    for (i in 1L:MCtimes) {
      samp.idx[[i]] = sample(1L:x.row, x.row, replace = TRUE)
    }
  }
  
  if (method == 'jagging') {
    for (i in 1L:MCtimes) {
      samp.idx[[i]] = sample(1L:x.col, floor(x.col * ratio))
    }
  }
  
  #plsdf = as.data.frame(cbind(x, y))
  if(method == 'mc' | method == 'bootstrap'){
     
    if (is.null(maxcomp)) maxcomp = ncol(x)
    if (parallel < 1.5) {
      predcvlist = vector('list', MCtimes)
      predtelist = vector('list', MCtimes)
      for (i in 1L:MCtimes) {
        if (verbose) cat('Beginning MCtimes', i, '\n')
        plsdf.x = x[samp.idx[[i]], ]
        plsdf.y = y[samp.idx[[i]]]
        predcvlist[[i]] = suppressWarnings(enpls.ad.cv(plsdf.x, plsdf.y, maxcomp = maxcomp))
        predtelist[[i]] = enpls.ad.core(plsdf.x, plsdf.y, x.test, maxcomp = maxcomp)
      }
      
    } else {
    
    registerDoParallel(parallel)
    predcvlist = foreach(i = 1L:MCtimes) %dopar% {
      plsdf.x = x[samp.idx[[i]], ]
      plsdf.y = y[samp.idx[[i]]]
      enpls.ad.cv(plsdf.x, plsdf.y, maxcomp = maxcomp)
     }
    
    predtelist = foreach(i = 1L:MCtimes) %dopar% {
      plsdf.x = x[samp.idx[[i]], ]
      plsdf.y = y[samp.idx[[i]]]
      enpls.ad.core(plsdf.x, plsdf.y, x.test, maxcomp = maxcomp)
    }
   }
   
   predcvmat = matrix(NA, ncol = x.row, nrow = MCtimes)
   for (i in 1L:MCtimes) {
     for (j in 1L:length(samp.idx[[i]])) {
       predcvmat[i, samp.idx[[i]][j]] = predcvlist[[i]][j]
     }
    }
   
   predtemat = matrix(NA, ncol = xte.row, nrow = MCtimes)
   for (i in 1L:MCtimes) {
       predtemat[i, ] = predtelist[[i]]
     }
  }
  

  if(method == 'jagging'){
   
   if (is.null(maxcomp)) maxcomp = length(samp.idx[[1]])
   if (parallel < 1.5) {
    predcvlist = vector('list', MCtimes)
    predtelist = vector('list', MCtimes)
    for (i in 1L:MCtimes) {
      if (verbose) cat('Beginning MCtimes', i, '\n')
      plsdf.x = x[, samp.idx[[i]]]
      plsdf.y = y
      predcvlist[[i]] = suppressWarnings(enpls.ad.cv(plsdf.x, plsdf.y, maxcomp = maxcomp))
      predtelist[[i]] = enpls.ad.core(plsdf.x, plsdf.y, x.test, maxcomp = maxcomp)
    }
    
   }else {
     
     registerDoParallel(parallel)
     predcvlist = foreach(i = 1L:MCtimes) %dopar% {
       plsdf.x = x[, samp.idx[[i]]]
       plsdf.y = y
       enpls.ad.cv(plsdf.x, plsdf.y, maxcomp = maxcomp)
     }
     
     predtelist = foreach(i = 1L:MCtimes) %dopar% {
       plsdf.x = x[, samp.idx[[i]]]
       plsdf.y = y
       enpls.ad.core(plsdf.x, plsdf.y, x.test, maxcomp = maxcomp)
     }
   }
   
   predcvmat = matrix(NA, ncol = x.row, nrow = MCtimes)
   for (i in 1L:MCtimes) {
       predcvmat[i, ] = predcvlist[[i]]
   }
   
   predtemat = matrix(NA, ncol = xte.row, nrow = MCtimes)
   for (i in 1L:MCtimes) {
     predtemat[i, ] = predtelist[[i]]
   }
   
  }
  
  
  
  ypredcv.mean = apply(predcvmat, 2L, mean, na.rm = TRUE)
  error.cv = abs(y - ypredcv.mean)
  
  STD.cv = c()
  for(i in 1L:x.row){
    
    subcv.ypred = na.omit(predcvmat[,i])
    STD.cv[i] = sqrt(sum((subcv.ypred - ypredcv.mean[i])^2)/(length(subcv.ypred) - 1))
  }

  
  ypredte.mean = apply(predtemat, 2L, mean, na.rm = TRUE)
  if (is.null(y.test)) error.te = NULL
  else error.te = abs(y.test - ypredte.mean)
  
  STD.te = c()
  for(i in 1L:xte.row){
    STD.te[i] = sqrt(sum((predtemat[,i] - ypredte.mean[i])^2)/(MCtimes - 1))
  }
  
  object = list('STD.cv' = STD.cv, 
                'STD.te' = STD.te, 
                'error.cv' = error.cv, 
                'error.te' = error.te)
  
  
  class(object) = 'enpls.ad'
  return(object)
  
}

#' cv function for enpls.od
#'
#' This function performs k-fold cross validation for 
#' ensemble partial least squares regression.
#'
#' @return 5-fold cross validation result predicted y 
#' 
#' @keywords internal

enpls.ad.cv = function(plsdf.x, plsdf.y, 
                       nfolds = 5L, 
                       maxcomp= NULL) {
  
  plsdf.x.row = nrow(plsdf.x)
  index = rep_len(1L:nfolds, plsdf.x.row)
  
  ypred = plsdf.y
  
  for (i in 1L:nfolds) {
    xtrain = plsdf.x[index != i, ]
    ytrain = plsdf.y[index != i]
    xtest  = plsdf.x[index == i, ]
    ytest  = plsdf.y[index == i]
      num = apply(xtrain, 2L, function(x) sum(x == 0))
      nan.ratio = num/dim(xtrain)[1]
      weight = which(nan.ratio > 0.7)
      if(length(weight > 0)) {
      xtrain = xtrain[, -weight]
      xtest = xtest[, -weight]
      if(maxcomp > ncol(xtrain)) maxcomp = ncol(xtrain)
      }
    plsr.cvfit = plsr(ytrain ~ ., data = data.frame(xtrain, ytrain), 
                      ncomp  = maxcomp, 
                      scale  = TRUE, 
                      method = 'simpls', 
                      validation = 'CV', segments = 5L)
    
    # choose best component number using adjusted CV
    cv.bestcomp = which.min(RMSEP(plsr.cvfit)[['val']][2L, 1L, -1L])
    
    plsr.fit = plsr(ytrain ~ ., data = data.frame(xtrain, ytrain), 
                    ncomp  = cv.bestcomp, 
                    scale  = TRUE, 
                    method = 'simpls', 
                    validation = 'none')
    ypredvec  = predict(plsr.fit, comps = 1:cv.bestcomp, xtest)
    ypred[index == i] = ypredvec
  }
  
  return(ypred)
  
}

#' core function for enpls.od
#'
#' select the best ncomp with cross-validation and
#' use it to fit the complete training set again, 
#' then predict on the test set. scale = TRUE
#'
#' @return predicted y for test set
#' 
#' @keywords internal
 

enpls.ad.core = function(plsdf.x, plsdf.y, x.test, maxcomp = NULL) {
  
  num = apply(plsdf.x, 2L, function(x) sum(x == 0))
  ratio = num/dim(plsdf.x)[1]
  weight = which(ratio > 0.7)
  if(length(weight > 0)) {
    plsdf.x = plsdf.x[, -weight]
    x.test = x.test[, -weight]
    if(maxcomp > ncol(plsdf.x)) maxcomp = ncol(plsdf.x)
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
  ypredvec  = predict(plsr.fit, comps = 1:cv.bestcomp, x.test)
  return(ypredvec)
}