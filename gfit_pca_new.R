gfit_pca_new <- function (target, sumstat, index, cprob = 0.1, xlim = NULL, ylim = NULL, 
          ...) 
{
  loc2plot = function(x, y, cprob, ...) {
    fit = locfit(~lp(x, y, nn = 0.2), maxk = 1000, mint = 100, 
                 maxit = 100)
    lev = sort(fitted(fit))[floor(cprob * length(x))]
    plot(fit, lev = lev, m = 100, drawlabels = FALSE, ...)
    return(list(fit = fit, lev = lev))
  }
  if (is.vector(target)) {
    target = t(as.data.frame(target))
    if (is.data.frame(sumstat)) {
      colnames(target) = names(sumstat)
    }
    if (is.matrix(sumstat)) {
      colnames(target) = colnames(sumstat)
    }
  }
  res.prcomp = prcomp(sumstat, scale = T, center = T)
  nmod = length(table(index))
  theindex = names(table(index))
  if (is.null(xlim)) {
    xlim = ylim
  }
  if (is.null(ylim)) {
    ylim = xlim
  }
  if (!is.null(xlim)) {
    plot(0, type = "n", xlim = xlim, ylim = ylim, xlab = "PC1", 
         ylab = "PC2")
  }
  for (i in 1:nmod) {
    ind = index == theindex[i]
    if ((i == 1) & (is.null(xlim))) {
      add = FALSE
    }
    else {
      add = TRUE
    }
    loc2plot(res.prcomp$x[ind, 1], res.prcomp$x[ind, 2], 
             cprob, col = i, lty = 1, lwd = 2, add = add, ...)
  }
  points(predict(res.prcomp, target)[1], predict(res.prcomp, 
                                                 target)[2], col = 1, cex = 2, pch = 3, lwd = 2)
  legend("topright", legend = theindex, cex = 1.5, col = c(1:nmod), 
         lty = 0, pch = 15)
}
