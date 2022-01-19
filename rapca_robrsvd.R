rapca_robrsvd = function (x, FUN = Qn, mean = TRUE) 
{
  X <- t(x)
  n <- nrow(X)
  p <- ncol(X)
  if (mean) {
    med <- colMeans(X)
    xx <- sweep(X, 2, med)
  }
  else xx <- X
  ###
  dum = basis_rsvd(xx)
  tmp <- matrix(dum$right_vector, ncol = dum$N_rsvd)
  if (mean) {
    med <- c(med + tmp[, 1])
    xx <- sweep(X, 2, med)
    basis <- cbind(med, tmp[, (1:dum$N_rsvd)])
    coef <- cbind(rep(1, n), xx %*% basis[, -1])
  }
  else {
    basis <- tmp
    coef <- xx %*% basis
  }
  return(list(basis = basis, coeff = coef, X = xx, order = dum$N_rsvd ))
}