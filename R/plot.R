#' Plot cv.enpls Object
#'
#' This function plots cv.enpls object.
#'
#' This function plots cv.enpls object.
#'
#' @param x An object of class \code{cv.enpls}.
#' @param main plot title
#' @param ... Other graphical parameters to be passed on to \code{plot}.
#'
#' @author  Min-feng Zhu <\email{wind2zhu@@163.com}>,
#'          Nan Xiao <\email{road2stat@@gmail.com}>
#'        
#'
#' @seealso See \code{\link{cv.enpls}} for performing ensemble PLS regression.
#'
#' @method plot cv.enpls
#'
#' @export
#'
#' @examples
#' data(alkanes)
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' cv.enpls.fit = cv.enpls(x, y, MCtimes = 10)
#' plot(cv.enpls.fit)

plot.cv.enpls = function(x, main = NULL, ...) {

  if (!inherits(x, 'cv.enpls'))
    stop('This function only works for objects of class "cv.enpls"')

  y.data = x$ypred
  y.real = y.data[, 'y.real']
  y.pred = y.data[, 'y.pred']

  plot(y.real, y.pred, 
       xlim = range(y.real), ylim = range(y.real), 
       xlab = 'Real Response', ylab = 'Predicted Response', 
       main = main, ...)
  abline(a = 0L, b = 1L)

}

#' Plot enpls.fs Object
#'
#' This function plots enpls.fs object.
#'
#' This function plots enpls.fs object.
#'
#' @param x An object of class \code{enpls.fs}.
#' @param sort Should the variables be sorted in decreasing order of importance?
#' @param nvar How many variables to show? Ignored if \code{sort = FALSE}.
#' @param main plot title
#' @param ... Other graphical parameters to be passed on to \code{dotchart}.
#'
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>,
#'         Nan Xiao <\email{road2stat@@gmail.com}>
#'         
#' @seealso See \code{\link{enpls.fs}} for feature selection with ensemble PLS.
#'
#' @method plot enpls.fs
#'
#' @export
#'
#' @examples
#' data(alkanes)
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' varimp = enpls.fs(x, y, MCtimes = 100)
#' plot(varimp)
#' plot(varimp, nvar = 10L)

plot.enpls.fs = function(x, 
                         sort = TRUE, nvar = NULL, 
                         main = NULL, ...) {

  if (!inherits(x, 'enpls.fs'))
    stop('This function only works for objects of class "enpls.fs"')

  varimp = x$variable.importance
  if (is.null(nvar)) nvar = length(varimp)

  if (sort == TRUE) {
    dotchart(sort(varimp, TRUE)[nvar:1], main = main, ...)
  } else {
    dotchart(rev(varimp), main = main, ...)
  }

}

#' Plot enpls.od Object
#'
#' This function plots enpls.od object.
#'
#' This function plots enpls.od object.
#'
#' @param x An object of class \code{enpls.od}.
#' @param criterion Criterion of being outlier, 
#' could be \code{'quantile'} or \code{'sd'}.
#' @param prob the quantile
#' @param sdtimes the times of sd
#' @param main plot title
#' @param ... Other graphical parameters to be passed on to \code{plot}.
#'
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>,
#'         Nan Xiao <\email{road2stat@@gmail.com}>
#'         
#' @seealso See \code{\link{enpls.od}} for outlier detection with ensemble PLS.
#'
#' @method plot enpls.od
#'
#' @export
#'
#' @examples
#' data(alkanes)
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' od = enpls.od(x, y, MCtimes = 100)
#' plot(od, criterion = 'quantile')
#' plot(od, criterion = 'sd')

plot.enpls.od = function(x, 
                         criterion = c('quantile', 'sd'), 
                         prob = 0.05, sdtimes = 3L, 
                         main = NULL, ...) {

  if (!inherits(x, 'enpls.od'))
    stop('This function only works for objects of class "enpls.od"')

  criterion = match.arg(criterion)

  error.mean = x$error.mean
  error.sd = x$error.sd

  if (criterion == 'quantile') {
    vpos = quantile(error.mean, 1 - prob)
    hpos = quantile(error.sd, 1 - prob)
  } else {
    vpos = mean(error.mean) + (sdtimes * sd(error.mean))
    hpos = mean(error.sd) + (sdtimes * sd(error.sd))
  }

  yout = intersect(which(error.mean >= vpos), which(error.sd <= hpos))
  Xout = intersect(which(error.mean <= vpos), which(error.sd >= hpos))
  abnormal = intersect(which(error.mean >= vpos), which(error.sd >= hpos))

  plot(error.mean, error.sd, 
       xlab = 'Error Mean', ylab = 'Error SD', main = main, ...)

  abline(h = hpos, col = 'gray', lty = 2)
  abline(v = vpos, col = 'gray', lty = 2)

  if (length(yout) != 0L) text(error.mean[yout], error.sd[yout], 
                               labels = as.character(yout), 
                               col = 'red', cex = 0.7, pos = 3)
  if (length(Xout) != 0L) text(error.mean[Xout], error.sd[Xout], 
                               labels = as.character(Xout), 
                               col = 'blue', cex = 0.7, pos = 1)
  if (length(abnormal) != 0L) text(error.mean[abnormal], error.sd[abnormal], 
                                   labels = as.character(abnormal), 
                                   col = 'purple', cex = 0.7, pos = 1)

}


#' Plot enpls.ad Object
#'
#' This function plots enpls.ad object.
#'
#' This function plots enpls.ad object.
#'
#' @param x An object of class \code{enpls.ad}.
#' @param main plot title
#' @param ... Other graphical parameters to be passed on to \code{ggplot}.
#'
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>,
#'         Nan Xiao <\email{road2stat@@gmail.com}>
#'         
#' @seealso See \code{\link{enpls.ad}} for applicability domain with ensemble PLS.
#'
#' @method plot enpls.ad
#'
#' @export
#'
#' @examples
#' data(logS)
#' x = logS$x
#' y = logS$y
#' x.test = logS$x.test1
#' y.test = logS$y.test1
#' 
#' set.seed(42)
#' ad = enpls.ad(x, y, x.test, y.test, MCtimes = 10)
#' plot(ad)
#' 

plot.enpls.ad = function(x, main = NULL, ...) {
  if (!inherits(x, 'enpls.ad'))
    stop('This function only works for objects of class "enpls.ad"')
    
    STD = c(x$STD.cv, x$STD.te)
    Error = c(x$error.cv, x$error.te)
    Lable = c(rep('cv', length(x$STD.cv)), rep('test',length(x$STD.te)))
    ad.gg = data.frame(STD = STD, Error = Error, Lable = Lable)
    ggplot(ad.gg, aes(x = STD, y = Error, colour = Lable, shape = Lable)) + geom_point() +
    scale_color_brewer(palette = 'Set2') + 
    theme_bw() + theme(legend.position = 'none') +
    xlab('STD') + ylab('Absolute Prediction Error') +  labs(title = main)
}

#' Plot enpls.pred Object
#'
#' This function plots enpls.pred object.
#'
#' This function plots enpls.pred object.
#'
#' @param x An object of class \code{enpls.pred}.
#' @param y.test Response vector for test.
#' @param main plot title
#' @param ... Other graphical parameters to be passed on to \code{ggplot}.
#'
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>,
#'         Nan Xiao <\email{road2stat@@gmail.com}>
#'         
#' @seealso See \code{\link{predict}} for  make predictions on new data by fitted enpls.en object.
#' @method plot enpls.pred
#'
#' @export
#'
#' @examples
#' data(logS)
#' x = logS$x
#' y = logS$y
#' x.test = logS$x.test1
#' y.test = logS$y.test1
#' 
#' set.seed(42)
#' enpls.fit = enpls.en(x, y, MCtimes = 10)
#' y.pred = predict(enpls.fit, newx = x.test)
#' plot(y.pred, y.test)
#' 

plot.enpls.pred = function(x, y.test, main = NULL, ...) {
   
  if (!inherits(x, 'enpls.pred'))
    stop('This function only works for objects of class "enpls.pred"')
  
  predmat = x$predmat
  ypred = x$ypred
  mx = nrow(x$predmat)
  yte = y.test
  se = apply(predmat, 1L, function(x){sqrt(sum((x - mean(x))^2)/100)})
  Ind = c(1 : mx)
  Ind = factor(Ind)
  if(max(y.test) > max(ypred)) 
    { ymax = max(y.test) + 1 
  } else { ymax = max(ypred) + 1 }
  if(min(y.test) < min(ypred)) 
  { ymin = min(y.test) - 1
  }else { ymin = min(ypred) - 1 }
  dm = data.frame(yte = yte, ypred = ypred, se = se,Ind = Ind)
  P = ggplot(dm, aes(x = yte, y = ypred, colour = Ind, group = Ind)) + 
    geom_errorbar(aes(ymin = ypred - se, ymax = ypred + se),
                  width = .2, size = 0.25) +
    geom_point(size = 1.5) 
  P + xlim(ymin, ymax) + ylim(ymin, ymax) +
    geom_abline(intercept = 0, slope = 1, colour = "#CCCCCC") +
    xlab("Experimental") + ylab("Predicted") +
    theme_bw() + 
    theme(legend.position="none")+
    labs(title = main)
}