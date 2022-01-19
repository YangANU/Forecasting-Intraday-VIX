ftsm_robrsvd <- function(y, ngrid = max(500, ncol(y$y)), mean = TRUE, level = FALSE, weight = FALSE)
{
  if (length(colnames(y$y)) > 0) {
    y$time = ts(as.numeric(colnames(y$y)), start = head(as.numeric(colnames(y$y)), 
                                                        1), end = tail(as.numeric(colnames(y$y)), 1), frequency = 1/diff(as.numeric(colnames(y$y)))[1])
  }
  else {
    y$time = ts(1:ncol(y$y), start = 1, end = ncol(y$y))
  }
  n <- ncol(y$y)
  m <- length(y$x)
  yy <- matrix(NA, nrow = ngrid, ncol = n)
  for (i in 1:n) {
    miss <- is.na(y$y[, i])
    yy[, i] <- spline(y$x[!miss], y$y[!miss, i], n = ngrid)$y
  }
  xx <- seq(min(y$x[!miss]), max(y$x[!miss]), l = ngrid)
  mean.se = approx(xx, sqrt(apply(yy, 1, var)/n), xout = y$x)$y
  ###
  ytsp <- tsp(y$time)
  y.robrsvd = fdrsvd(y$x, y$y, ngrid = ngrid, mean = mean, level = level)
  order = y.robrsvd$order
  coeff <- ts(y.robrsvd$coeff, start = ytsp[1], frequency = ytsp[3])
  ###
  basis <- y.robrsvd$basis
  fits <- fts(1:length(y$x), basis %*% t(coeff), start = ytsp[1], 
              frequency = ytsp[3], xname = y$xname, yname = paste("Fitted", 
                                                                  y$yname))
  fits$x = y$x
  if (mean) {
    colnames(basis) = c("mean", paste("phi", 1:order, 
                                        sep = ""))
    colnames(coeff) = c("mean", paste("beta", 1:order, 
                                        sep = ""))
  }
  else {
    colnames(basis) = c(paste("phi", 1:order, sep = ""))
    colnames(coeff) = c(paste("beta", 1:order, sep = ""))
  }
  res <- fts(y$x, y$y - fits$y, start = ytsp[1], frequency = ytsp[3], 
             xname = y$xname, yname = paste("Residuals", y$yname))

  out <- list(x1 = as.numeric(colnames(y$y)), y1 = as.numeric(rownames(y$y)), 
              y = fts(y$x, y$y, xname = y$xname, yname = y$yname), 
              basis = basis, coeff = coeff, fitted = fits, residuals = res, 
              v = ts(y.robrsvd$v, start = ytsp[1], 
              frequency = ytsp[3]), basis2 = y.robrsvd$basis2, 
              coeff2 = y.robrsvd$coeff2, mean.se = mean.se, call = match.call())
  
  return(structure(out, class = c("ftsm", "fm")))
}