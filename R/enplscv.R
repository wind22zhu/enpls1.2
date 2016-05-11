#' Ensemble Partial Least Squares Cross-Valdidation for feature selection
#'
#' This function shows the cross-validated prediction performance of
#' models with sequentially reduced number of predictors (ranked by
#' variable importance) via a nested cross-validation procedure.
#'
#' This function shows the cross-validated prediction performance of
#' models with sequentially reduced number of predictors (ranked by
#' variable importance) via a nested cross-validation procedure.
#'
#' @param trainx matrix or data frame containing columns of predictor
#' @param trainy vector of response, must have length equal to the number
#' of rows in \code{trainx}.
#' @param cv.fold number of folds in the cross-validation
#' @param scale if \code{"log"}, reduce a fixed proportion (\code{step})
#' of variables at each step, otherwise reduce \code{step} variables at a
#' time.
#' @param step if \code{log=TRUE}, the fraction of variables to remove at
#' each step, else remove this many variables at a time
#' @param recursive whether variable importance is (re-)assessed at each
#' step of variable reduction 
#' @param MCtimes times of Monte-Carlo
#' @param ... other arguments passed on to \code{enpls.en}
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{n.var} - vector of number of variables used at each step.
#' \item \code{error.cv} - corresponding vector of error rates or MSEs at each step.
#' \item \code{predicted} - list of \code{n.var} components, each containing
#' the predicted values from the cross-validation.
#' \item \code{res} - list of \code{n.var} components, each containing
#' the sum of coefficient values from the cross-validation.
#' \item \code{imp} - list of \code{n.var} components, each containing
#' the sum of coefficient values from the cross-validation.
#' }
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>,
#'         Nan Xiao <\email{road2stat@@gmail.com}>
#'  
#' @seealso See \code{\link{enpls.en}} for feature selection with ensemble PLS. 
#' See \code{\link{enpls.en}} for ensemble PLS regression              
#'
#' @references Svetnik, V., Liaw, A., Tong, C. and Wang, T., "Application of Breiman's
#' Random Forest to Modeling Structure-Activity Relationships of
#' Pharmaceutical Molecules", MCS 2004, Roli, F. and Windeatt, T. (Eds.)
#' pp. 334-343.
#'
#' @examples
#' data(logS)
#' x = logS$x
#' y = logS$y
#' 
#' result = enplscv(x, y, recursive = TRUE)
#' with(result, plot(n.var, error.cv, type="b", lwd=2))
#'

enplscv <- function(trainx, trainy, cv.fold=5, scale="log", step=0.5,
                  recursive = FALSE, MCtimes = 10, ...) {
  classRF <- is.factor(trainy)
  n <- nrow(trainx)
  p <- ncol(trainx)
  if (scale == "log") {
    k <- floor(log(p, base=1/step))
    n.var <- round(p * step^(0:(k-1)))
    same <- diff(n.var) == 0
    if (any(same)) n.var <- n.var[-which(same)]
    if (1 %in% n.var) n.var <- n.var[-which(n.var ==1 )]
  } else {
    n.var <- seq(from=p, to=1, by=step)
    if (1 %in% n.var) n.var <- n.var[-which(n.var ==1 )]
  }
  k <- length(n.var)
  cv.pred <- vector(k, mode="list")
  imp.var.num <- vector(k, mode="list")
  imp.var.idx <- vector(k, mode="list")
  for (i in 1:k) cv.pred[[i]] <- trainy
  for (i in 1:k) imp.var.num[[i]] <- matrix(NA, cv.fold, n.var[i])
  for (i in 1:k) imp.var.idx[[i]] <- matrix(NA, cv.fold, n.var[i])
  ## Generate the indices of the splits
  ## For regression, bin the response into 5 bins and stratify.
    f <- factor(rep(1:5, length=length(trainy))[order(order(trainy))])
  
  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <-  sample(rep(1:cv.fold, length=nlvl[i]))
  }
  
  for (i in 1:cv.fold) {
    ## cat(".")
    all.enpls <- enpls.en(trainx[idx != i, , drop=FALSE],
                           trainy[idx != i], MCtimes = MCtimes)
    
    cv.pred[[1]][idx == i] <- predict(all.enpls, newx = trainx[idx == i, , drop=FALSE])$ypred
    all.imp <- enpls.fs(trainx[idx != i, , drop=FALSE],trainy[idx != i], MCtimes = 10)
    
    impvar <- (1:p)[order(all.imp$variable.importance, decreasing=TRUE)]
    imp.var.num[[1]][i,] <- as.numeric(all.imp$variable.importance)
    imp.var.idx[[1]][i,] <- impvar
    impidx <- impvar
    for (j in 2:k) {
      imp.idx <- impvar[1:n.var[j]]
      sub.enpls <- enpls.en(trainx[idx != i, imp.idx, drop=FALSE],
                            trainy[idx != i], MCtimes = MCtimes)
      cv.pred[[j]][idx == i] <-  predict(sub.enpls,newx = trainx[idx == i,imp.idx , drop=FALSE])$ypred
      sub.imp <- enpls.fs(trainx[idx != i,imp.idx , drop=FALSE],trainy[idx != i], MCtimes = 10)
      
      imp.var.num[[j]][i,] <- as.numeric(sub.imp$variable.importance)
      
      imp <- (1:length(imp.idx))[order(sub.imp$variable.importance, decreasing=TRUE)]
 
      impvar.idx <- c()
      for (l in 1L:n.var[j]){
        impvar.idx[l] = impidx[imp[l]]
      } 
      impidx <- impvar.idx
  ## For recursive selection, use importance measures from the sub-model.
      
    if (recursive) {
      impvar <- impidx
    }
      imp.var.idx[[j]][i,] <- impvar.idx
      
    }
    NULL
  }
  
  imp.var = vector(k, mode="list")
  for (i in 1:k) imp.var[[i]] <- matrix(0, cv.fold, p)
  imp.var[[1]] = imp.var.num[[1]]
  for (i in 2L:k) {
    for (j in 1L:cv.fold) {
      for (l in 1L:n.var[i]){
        imp.var[[i]][j, imp.var.idx[[i]][j,l]] = imp.var.num[[i]][j,l]
      }
    }
  }
  ## cat("\n")
  res <- vector(k, mode="list")
  for(i in 1:k) {
    res[[i]] <- apply(imp.var[[i]], 2L, sum)
  }
  
  
  if(classRF) {
    error.cv <- sapply(cv.pred, function(x) mean(trainy != x))
  } else {
    error.cv <- sapply(cv.pred, function(x) mean((trainy - x)^2))
  }
  names(error.cv) <- names(cv.pred) <- n.var
  list(n.var=n.var, error.cv=error.cv, predicted=cv.pred, 
       res=res)
}
