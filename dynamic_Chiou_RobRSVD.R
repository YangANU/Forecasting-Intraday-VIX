dynamic_Chiou_robrsvd = function(dat, newdata, holdoutdata, order_k_percent = 0.9, order_m_percent = 0.9, 
                                 robust_lambda = 2.33, bootrep = 100, pointfore)
{
  if(ncol(newdata) == 1)
  {
     stop("Can't perform RobRSVD")
     # data_first = t(scale(dat[1,], scale = FALSE))
     # data_last  = t(scale(t(dat[(ncol(newdata) + 1):nrow(dat),]), center = TRUE, scale = FALSE))
  }
  if(ncol(holdoutdata) == 1)
  {
     stop("Can't perform RobRSVD")
     # data_first = t(scale(t(dat[1:ncol(newdata),]), center = TRUE, scale = FALSE))
     # data_last  = t(scale(dat[nrow(dat),], center = TRUE, scale = FALSE))
  }
  if(ncol(newdata) > 1 & ncol(holdoutdata) > 1)
  {
      data_first = t(scale(t(dat[1:ncol(newdata),]), center = TRUE, scale = FALSE))
      data_last  = t(scale(t(dat[(ncol(newdata)+1):nrow(dat),]), center = TRUE, scale = FALSE))
  }
  cross_cov  = data_last %*% t(data_first)
  
  ################
  # RobRSVD (1st)
  ################
  
  dum_1st = basis_rsvd(crossprod(t(data_first)), prop = order_k_percent)
  phi_k = matrix(dum_1st$right_vector)
  xi_k  = t(data_first) %*% phi_k
  reconstruct_k = phi_k %*% t(xi_k)
  resi_k = data_first - reconstruct_k
  lambda_k = as.numeric(dum_1st$singular)
  
  if(pointfore == FALSE)
  {
    n_xi = nrow(xi_k)
    p_xi = ncol(xi_k)
    out_xi <- meboot(ts(as.numeric(xi_k),start=1,end=n_xi*p_xi), reps = bootrep)
    
    n_resi = nrow(t(resi_k))  
    p_resi = ncol(t(resi_k))
    out_resi <- meboot(ts(as.numeric(t(resi_k)),start=1,end=n_resi*p_resi), reps = bootrep)
    
    data_first_boot = array(, dim = c(nrow(data_first), ncol(data_first), bootrep))
    phi_k_boot = array(, dim = c(nrow(phi_k), ncol(phi_k), bootrep))
    lambda_k_boot = matrix(, length(lambda_k), bootrep)
    for(ik in 1:bootrep)
    {
      xi_k_boot = matrix(out_xi$ensemble[,ik], n_xi, p_xi)
      resi_k_boot = t(matrix(out_resi$ensemble[,ik], n_resi, p_resi))
      if(!any(is.finite(resi_k_boot)))
      {
        data_first_boot[,,ik] = phi_k %*% t(xi_k_boot)
      }  
      else
      {
        data_first_boot[,,ik] = phi_k %*% t(xi_k_boot) + resi_k_boot        
      }
      phi_k_boot[,,ik] = matrix(basis_rsvd(crossprod(t(data_first_boot[,,ik])), prop = order_k_percent)$right_vector)
      lambda_k_boot[,ik] = as.numeric(basis_rsvd(crossprod(t(data_first_boot[,,ik])), prop = order_k_percent)$singular)
    }
  }
  
  ################
  # RobRSVD (2nd)
  ################
  
  dum_2nd = basis_rsvd(crossprod(t(data_last)), prop = order_k_percent)
  psi_m = matrix(dum_2nd$right_vector)
  eta_m = t(data_last) %*% psi_m
  reconstruct_m = psi_m %*% t(eta_m)
  resi_m = data_last - reconstruct_m
  if(pointfore == FALSE)
  {
    n_eta = nrow(eta_m)
    p_eta = ncol(eta_m)
    out_eta = meboot(ts(as.numeric(eta_m),start=1,end=n_eta*p_eta), reps = bootrep)
    
    n_resi = nrow(t(resi_m))
    p_resi = ncol(t(resi_m))
    out_resi = meboot(ts(as.numeric(t(resi_m)),start=1,end=n_resi*p_resi), reps = bootrep)
    
    data_last_boot = array(, dim = c(nrow(data_last), ncol(data_last), bootrep))
    psi_m_boot = array(, dim = c(nrow(psi_m), ncol(psi_m), bootrep))
    for(ik in 1:bootrep)
    {
      eta_m_boot = matrix(out_eta$ensemble[,ik], n_eta, p_eta)
      resi_m_boot = t(matrix(out_resi$ensemble[,ik], n_resi, p_resi))
      if(!any(is.finite(resi_m_boot)))
      {
        data_last_boot[,,ik] = psi_m %*% t(eta_m_boot)
      }
      else
      {
        data_last_boot[,,ik] = psi_m %*% t(eta_m_boot) + resi_m_boot
      }
      psi_m_boot[,,ik] = matrix(basis_rsvd(crossprod(t(data_last_boot[,,ik])), prop = order_k_percent)$right_vector)
    }
  }
  
  if(pointfore == TRUE)
  {
    order_k = length(as.numeric(dum_1st$singular))
    order_m = length(as.numeric(dum_2nd$singular))
    lambda_km = matrix(, order_k, order_m)
    for(j in 1:order_k)
    {
      for(i in 1:order_m)
      {
        lambda_km[j,i] = matrix(psi_m[,i], nrow = 1) %*% cross_cov %*% matrix(phi_k[,j], nrow = length(newdata))
      }
    }
    
    beta = array(, dim = c(length(newdata), length(holdoutdata), order_m, order_k))
    for(j in 1:order_k)
    {
      for(i in 1:order_m)
      {
        beta[,,i,j] = lambda_km[j,i]/lambda_k[j] * matrix(phi_k[,j], ncol = 1) %*% matrix(psi_m[,i], nrow = 1)
      }
    }
    
    finalbeta = apply(beta, c(1,2), sum)
    
    if(ncol(newdata) == 1)
    {
      update_forecast = rowMeans(dat[(ncol(newdata)+1):nrow(dat),]) + matrix(newdata - mean(dat[1:ncol(newdata),]), nrow=1)%*%finalbeta
    }  
    if(ncol(holdoutdata) == 1)
    {
      update_forecast = mean(dat[(ncol(newdata)+1):nrow(dat),]) + matrix(newdata -rowMeans(dat[1:ncol(newdata),]), nrow = 1)%*%finalbeta
    }
    if(ncol(newdata) > 1 & ncol(holdoutdata) > 1)
    {
      update_forecast = rowMeans(dat[(ncol(newdata)+1):nrow(dat),]) + matrix(newdata -rowMeans(dat[1:ncol(newdata),]), nrow = 1)%*%finalbeta
    }
    err = matrix(c(ftsa:::mae(update_forecast, holdoutdata), ftsa:::mse(update_forecast, holdoutdata)),nrow=1)
    colnames(err) = c("MAFE","MSFE")
    return(list(update_forecast = update_forecast, holdoutdata = holdoutdata, err = err, order_k = order_k, order_m = order_m))
  }
  else
  {
    lambda_km_boot = array(, dim = c(order_k, order_m, bootrep))
    for(ik in 1:bootrep)
    {
      for(j in 1:order_k)
      {
        for(i in 1:order_m)
        {
          cross_cov  = data_last_boot[,,ik] %*% t(data_first_boot[,,ik])
          lambda_km_boot[j,i,ik] = matrix(psi_m_boot[,i,ik], nrow = 1) %*% cross_cov %*% matrix(phi_k_boot[,j,ik], nrow = length(newdata))
        }
      }
    }
    
    beta_boot = array(, dim = c(length(newdata), length(holdoutdata), order_m, order_k, bootrep))
    for(ik in 1:bootrep)
    {
      for(j in 1:order_k)
      {
        for(i in 1:order_m)
        {
          beta_boot[,,i,j,ik] = lambda_km_boot[j,i,ik]/lambda_k_boot[j,ik] * matrix(phi_k_boot[,j,ik], ncol = 1) %*% matrix(psi_m_boot[,i,ik], nrow = 1)
        }
      }
    }
    
    update_forecast = matrix(, nrow(dat) - length(newdata), bootrep)
    for(ik in 1:bootrep)  
    {
      finalbeta_boot = apply(beta_boot[,,,,ik], c(1,2), sum)
      if(ncol(newdata) == 1)
      {
        update_forecast[,ik] = rowMeans(dat[(ncol(newdata)+1):nrow(dat),]) + matrix(newdata - mean(dat[1:ncol(newdata),]), nrow=1)%*%finalbeta_boot
      }  
      if(ncol(holdoutdata) == 1)
      {
        finalbeta_boot = as.matrix(rowSums(beta_boot[,,,,ik]))
        update_forecast[,ik] = mean(dat[(ncol(newdata)+1):nrow(dat),]) + matrix(newdata -rowMeans(dat[1:ncol(newdata),]), nrow = 1)%*%finalbeta_boot
      }
      if(ncol(newdata) > 1 & ncol(holdoutdata) > 1)
      {
        update_forecast[,ik] = rowMeans(dat[(ncol(newdata)+1):nrow(dat),]) + matrix(newdata -rowMeans(dat[1:ncol(newdata),]), nrow = 1)%*%finalbeta_boot
      }
    }
    dum = err_comp(dat, ik = nrow(dat) - length(holdoutdata), bootrep = bootrep)
    update_comb = update_forecast + dum
    update_comb_lb_ub = apply(update_comb, 1, quantile, c(0.1,0.9))
    return(update_comb_lb_ub)  
  } 
}
