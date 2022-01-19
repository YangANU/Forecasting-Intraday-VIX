###############################################################################################
# basis selection using robust regularized singular value decomposition by Zhang et al. (2013)
###############################################################################################

basis_rsvd = function(y, prop = 0.9, ngrid = max(500, ncol(y)), iugcv = FALSE, ivgcv = FALSE)
{
  my = colMeans(y, na.rm = TRUE)
  decenter_y = sweep(y, 2, my)
  dummy = RobRSVD_ginv(decenter_y, irobust = TRUE, iugcv = iugcv, ivgcv = ivgcv)
  s = dummy$s
  v = dummy$v
  u = dummy$u
  res = decenter_y - dummy$s * dummy$u %*% t(dummy$v)
  var_sum = (dummy$s)^2
  order = 1
  
  repeat
  {
    my = colMeans(res, na.rm = TRUE)
    decenter_y = sweep(res, 2, my)
    dummy = RobRSVD_ginv(decenter_y, irobust = TRUE, iugcv = iugcv, ivgcv = ivgcv)
    s = c(s, dummy$s)
    v = cbind(v, dummy$v)
    u = cbind(u, dummy$u)
    res = decenter_y - dummy$s * dummy$u %*% t(dummy$v)
    var_sum = var_sum + (dummy$s)^2
    order = order + 1
    if((dummy$s)^2/var_sum < (1-prop)/prop)
    {
      break
    }
  }
  
  res_second_last = sweep(decenter_y, 2, -my)
  return(list(N_rsvd = order-1, singular = s[1:(order-1)]  , left_vector = u[,1:(order-1)], right_vector = v[, 1:(order-1)], residual = res_second_last))
}


test = fts(1:100, data)


basis_rsvd(test)


temp = data - data.robrsvd$s * data.robrsvd$u %*% t(data.robrsvd$v)
