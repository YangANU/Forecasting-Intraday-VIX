L1median2 = ftsa:::L1median2

fdrsvd = function (x, y, ngrid = 500, mean = mean, level = level, ...) 
{
  xnam = x
  coefnam = colnames(y)
  x = 1:length(x)
  n <- ncol(y)
  m <- length(x)
  if (ngrid < n) 
    stop("Grid should be larger than number of observations per time period.")
  if (m != nrow(y)) 
    stop("x and y of incompatible dimension")
  yy <- matrix(NA, nrow = ngrid, ncol = n)
  for (i in 1:n) {
    miss <- is.na(y[, i])
    yy[, i] <- spline(x[!miss], y[!miss, i], n = ngrid)$y
  }
  xx <- seq(min(x[!miss]), max(x[!miss]), l = ngrid)
  delta <- xx[2] - xx[1]
  if (mean) {
      ax <- L1median2(t(yy), method = "hoss")
      yy <- sweep(yy, 1, ax)
      axse <- approx(xx, sqrt(apply(yy, 1, var)/n), xout = x)$y
  }
  else axse <- NULL
  if (level) {
    bx <- colMeans(yy, na.rm = TRUE)
    yy <- sweep(yy, 2, bx)
  }
  if (level) {
    coeff <- as.matrix(bx)
    basis <- as.matrix(rep(1, m))
    colnames(coeff) <- colnames(basis) <- "level"
  }
  else coeff <- basis <- NULL
  if (mean) {
    coeff <- cbind(rep(1, n), coeff)
    basis <- cbind(approx(xx, ax, xout = x)$y, basis)
    colnames(coeff)[1] <- colnames(basis)[1] <- "mean"
  }
  ####
  robusteig <- rapca_robrsvd(yy, mean = FALSE)
  B <- robusteig$coeff
  Phi <- robusteig$basis
  yyhat <- Phi %*% t(B)
  v <- colSums((yy - yyhat)^2) * delta
  order = robusteig$order
  ###
  Phinorm = matrix(NA, length(x[!miss]), order)
  Phinormngrid = matrix(NA, ngrid, order)
  for (i in 1:order) {
    Phinorm[, i] = approx(xx, Phi[, i], xout = x[!miss])$y/delta/(sqrt(sum((approx(xx, 
                                                                                   Phi[, i], xout = x[!miss])$y/delta)^2)))
    Phinormngrid[, i] = approx(x[!miss], Phinorm[, i], xout = xx)$y
  }
  B <- t(yy) %*% Phinormngrid
  v <- colSums((yy - Phinormngrid %*% t(B))^2) * delta
  colnames(B) <- paste("beta", 1:order, sep = "")
  coeffdummy = B * delta
  colmeanrm = matrix(colMeans(coeffdummy), dim(B)[2], 1)
  coeff <- cbind(coeff, sweep(coeffdummy, 2, colmeanrm))
  m <- ncol(basis)
  if (mean == TRUE) {
    basis = basis[!miss] + Phinorm %*% colmeanrm
    colnames(basis) = "mean"
  }
  for (i in 1:order) {
    basis <- cbind(basis, Phinorm[, i])
    if (sum(basis[, i + m]) < 0) {
      basis[, i + m] <- -basis[, i + m]
      coeff[, i + m] <- -coeff[, i + m]
    }
  }
  colnames(basis)[m + (1:order)] <- paste("phi", 1:order, sep = "")

  yy <- yy - Phinormngrid %*% t(B)
  s <- try(La.svd(t(yy)), silent = TRUE)
  if (class(s) == "try-error") {
    s <- svd(t(yy), LINPACK = TRUE)
    s$vt <- t(s$v)
  }
  Phi2 <- as.matrix(t(s$vt)[, s$d > 1e-06])
  m <- ncol(Phi2)
  basis2 <- coeff2 <- NULL
  Phinorm2 = matrix(NA, length(x[!miss]), m)
  Phinorm2ngrid = matrix(NA, ngrid, m)
  if (m > 0) {
    for (i in 1:m) {
      Phinorm2[, i] = approx(xx, Phi2[, i], xout = x[!miss])$y/delta/(sqrt(sum((approx(xx, 
                                                                                       Phi2[, i], xout = x[!miss])$y/delta)^2)))
      Phinorm2ngrid[, i] = approx(x[!miss], Phinorm2[, 
                                                     i], xout = xx)$y
    }
    B2 <- t(yy) %*% Phinorm2ngrid
    colnames(B2) <- paste("beta", order + (1:ncol(B2)), sep = "")
    coeff2dummy <- B2 * delta
    colmeanrm2 = matrix(colMeans(coeff2dummy), dim(B2)[2], 
                        1)
    coeff2 = sweep(coeff2dummy, 2, colmeanrm2)
    for (i in 1:m) {
      basis2 <- cbind(basis2, Phinorm2[, i])
      if (sum(basis2[, i]) < 0) {
        basis2[, i] <- -basis2[, i]
        coeff2[, i] <- -coeff2[, i]
      }
    }
    colnames(basis2) <- paste("phi", order + (1:m), sep = "")
  }
  basiscomb = matrix(, nrow(y), ncol(basis))
  basiscomb[!miss, ] = basis
  basis2comb = matrix(, nrow(y), ncol(basis2))
  basis2comb[!miss, ] = basis2
  return(list(naindex = as.numeric(which(miss)), basis = basiscomb, order = order,
              coeff = coeff, v = v, basis2 = basis2comb, coeff2 = coeff2, mean.se = axse))
}