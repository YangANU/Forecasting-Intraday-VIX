####################
# TS (non-updating)
####################

TS_noupdate <- function(data_return, pcdmethod = c("classical","M"), jk)
{
  pcdmethod = match.arg(pcdmethod)
  if(pcdmethod == "classical")
  {
    # creating storage for errors
    ts_point_update_mse = vector(, (length(trading_time)-2))
    ts_point_update_mae = vector(, (length(trading_time)-2))
    
    # doing functional PCA
    fts_object = fts(1:1621, data_return[,1:(jk+84)])
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = pcdmethod, lambda = 2.33)$varprop)>=0.9),1)
    for(i in 1:(length(trading_time)-2))
    {
      dum = dynupdate(data = fts_object, newdata = data_return[1:(i+1),(jk+85)], holdoutdata = data_return[(i+2):1621,(jk+85)],  	
                              method = "ts", pcdmethod = pcdmethod, order = ftsm_object_order, interval = FALSE) 
      ts_point_update_mse[i] = dum$errormse
      ts_point_update_mae[i] = dum$errormae
    }
    return(list(ts_point_update_mae = ts_point_update_mae, ts_point_update_mse = ts_point_update_mse))
  }
  if(pcdmethod == "M")
  {
    # creating storage for errors
    ts_point_update_M_mse = vector(, (length(trading_time)-2))
    ts_point_update_M_mae = vector(, (length(trading_time)-2))
    
    # doing functional PCA
    fts_object = fts(1:1621, data_return[,1:(jk+84)])
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = pcdmethod, lambda = 2.33)$varprop)>=0.9),1)
    for(i in 1:(length(trading_time)-2))
    {
      dum = dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(jk+85)],  	
                        method = "ts", pcdmethod = pcdmethod, order = ftsm_object_order, robust_lambda = 2.33, interval = FALSE) 
      ts_point_update_M_mse[i] = dum$errormse
      ts_point_update_M_mae[i] = dum$errormae
    }
    return(list(ts_point_update_M_mae = ts_point_update_M_mae, ts_point_update_M_mse = ts_point_update_M_mse))
  }
}


# pcdmethod = "classical"
htm = proc.time()

library(doParallel)
cl <- makeCluster(10) 
registerDoParallel(cl)

ts_point_update_err_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% TS_noupdate(fts_return, "classical", jk = j)

# record errors
ts_point_update_mae = ts_point_update_mse = matrix(, 1619, 40)
for (j in 1:40)
{
  ts_point_update_mae[,j] = ts_point_update_err_list[[j]]$ts_point_update_mae
  ts_point_update_mse[,j] = ts_point_update_err_list[[j]]$ts_point_update_mse
}

ts_point_update_err = cbind(round(rowMeans(ts_point_update_mse), 4), round(rowMeans(ts_point_update_mae), 4))
colnames(ts_point_update_err) = c("mse", "mae")

# pcdmethod = "M"

ts_point_update_M_err_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% TS_noupdate(fts_return, "M", jk = j)

# record errors
ts_point_update_M_mae = ts_point_update_M_mse = matrix(, 1619, 40)
for (j in 1:40)
{
  ts_point_update_M_mae[,j] = ts_point_update_M_err_list[[j]]$ts_point_update_M_mae
  ts_point_update_M_mse[,j] = ts_point_update_M_err_list[[j]]$ts_point_update_M_mse
}

ts_point_update_M_err = cbind(round(rowMeans(ts_point_update_M_mse), 4), round(rowMeans(ts_point_update_M_mae), 4))
colnames(ts_point_update_M_err) = c("mse", "mae")

proc.time() - htm

##############
# BM updating
##############

htm = proc.time()

BM_update <- function(data_return, pcdmethod = c("classical","M"), jk)
{
  pcdmethod = match.arg(pcdmethod)
  if(pcdmethod == "classical")
  {
    # creating storage for errors
    bm_point_update_mse = vector(, (length(trading_time)-2))
    bm_point_update_mae = vector(, (length(trading_time)-2))
    
    # doing functional PCA
    fts_object = fts(1:1621, data_return[,1:(jk+84)])
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "classical")$varprop)>=0.9),1)
    for(i in 1:(length(trading_time)-2))
    {
      dum = dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(j+85)],  	
                      method = "block", pcdmethod = "classical", order = ftsm_object_order, interval = FALSE) 
      bm_point_update_mse[i] = dum$errormse
      bm_point_update_mae[i] = dum$errormae
    }
    return(list(bm_point_update_mae = bm_point_update_mae, bm_point_update_mse = bm_point_update_mse))
  }
  
  if(pcdmethod == "M")
  {
    # creating storage for errors
    bm_point_update_M_mse = vector(, (length(trading_time)-2))
    bm_point_update_M_mae = vector(, (length(trading_time)-2))
    
    # doing functional PCA
    fts_object = fts(1:1621, data_return[,1:(jk+84)])
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "M", lambda=2.33)$varprop)>=0.9),1)
    for(i in 1:(length(trading_time)-2))
    {
      dum = dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(jk+85)],  	
                      method = "block", pcdmethod = "M", order = ftsm_object_order, robust_lambda = 2.33, interval = FALSE)
      bm_point_update_M_mse[i] = dum$errormse
      bm_point_update_M_mae[i] = dum$errormae
    }
    return(list(bm_point_update_M_mae = bm_point_update_M_mae, bm_point_update_M_mse = bm_point_update_M_mse))
  }
}

# pcdmethod = "classical"

library(doParallel)
cl <- makeCluster(10) 
registerDoParallel(cl)

bm_point_update_err_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% BM_update(fts_return, "classical", jk = j)

# record errors
bm_point_update_mae = bm_point_update_mse = matrix(, 1619, 40)
for (j in 1:40)
{
  bm_point_update_mae[,j] = bm_point_update_err_list[[j]]$bm_point_update_mae
  bm_point_update_mse[,j] = bm_point_update_err_list[[j]]$bm_point_update_mse
}

bm_point_update_err = cbind(round(rowMeans(bm_point_update_mse), 4), round(rowMeans(bm_point_update_mae), 4))
colnames(bm_point_update_err) = c("mse", "mae")

# pcdmethod = "M"

library(doParallel)
cl <- makeCluster(10) 
registerDoParallel(cl)

bm_point_update_M_err_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% BM_update(fts_return, "M", jk = j)

# record errors
bm_point_update_M_mae = bm_point_update_M_mse = matrix(, 1619, 40)
for (j in 1:40)
{
  bm_point_update_M_mae[,j] = bm_point_update_M_err_list[[j]]$bm_point_update_M_mae
  bm_point_update_M_mse[,j] = bm_point_update_M_err_list[[j]]$bm_point_update_M_mse
}

bm_point_update_M_err = cbind(round(rowMeans(bm_point_update_M_mse), 4), round(rowMeans(bm_point_update_M_mae), 4))
colnames(bm_point_update_M_err) = c("mse", "mae")

proc.time() - htm

###############
# PLS updating
###############

# function for selecting lambda

pls_find_lambda <- function(find_lambda, pcdmethod = c("classical","M"), 
                            error = c("mse", "mae", "score"), ik)
{
  pcdmethod = match.arg(pcdmethod)
  error = match.arg(error)
  if(error == "mse")
  {
    pls_point_update_mse_train = vector(,42)
    for(j in 1:42)
    {
      fts_object = fts(1:1621, fts_return[,1:(j+42)])
      ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = pcdmethod, lambda=2.33)$varprop)>=0.9),1)
      dum = dynupdate(data = fts_object, newdata = fts_return[1:ik,(j+43)], holdoutdata =fts_return[(ik+1):1621,(j+43)],  	
                      method = "pls", pcdmethod = pcdmethod, order = ftsm_object_order, interval = FALSE, lambda = find_lambda) 
      pls_point_update_mse_train[j] = dum$errormse
    }
    return(mean(pls_point_update_mse_train))
  }
  if(error == "mae")
  {
    pls_point_update_mae_train = vector(,42)
    for(j in 1:42)
    {
      fts_object = fts(1:1621, fts_return[,1:(j+42)])
      ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = pcdmethod, lambda=2.33)$varprop)>=0.9),1)
      dum = dynupdate(data = fts_object, newdata = fts_return[1:ik,(j+43)], holdoutdata =fts_return[(ik+1):1621,(j+43)],  	
                      method = "pls", pcdmethod = pcdmethod, order = ftsm_object_order, interval = FALSE, lambda = find_lambda) 
      pls_point_update_mae_train[j] = dum$errormae
    }
    return(median(pls_point_update_mae_train))
  }
  if(error == "score")
  {
    pls_interval_update_score_train = vector(,42)
    for(j in 1:42)
    {
      holdout = as.numeric(fts_return[(ik+1):1621,(j+43)])
      fts_object = fts(1:1621, fts_return[,1:(j+42)])
      ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = pcdmethod, lambda=2.33)$varprop)>=0.9),1)
      dum = dynupdate(data = fts_object, newdata = fts_return[1:ik,(j+43)], holdoutdata = holdout,
                      method = "pls", pcdmethod = pcdmethod, order = ftsm_object_order, interval = TRUE, lambda = find_lambda)
      dum_lb = as.numeric(dum$low$y)
      dum_ub = as.numeric(dum$up$y)
      dum_lb_ind = ifelse(holdout < dum_lb, 1, 0)
      dum_ub_ind = ifelse(holdout > dum_ub, 1, 0)
      pls_interval_update_score_train[j] = mean((dum_ub - dum_lb) + 2/0.2 * ((dum_lb - holdout)*dum_lb_ind +  (holdout - dum_ub)*dum_ub_ind))
    }
    return(mean(pls_interval_update_score_train))  
  }
}

# function for doing PLS updating

PLS_update = function(data_return, pcdmethod = c("classical","M"), jk)
{
  pcdmethod = match.arg(pcdmethod)
  if(pcdmethod == "classical")
  {
    # creating storage for errors
    pls_point_update_mse = vector(, (length(trading_time)-2))
    pls_point_update_mae = vector(, (length(trading_time)-2))
    
    # doing functional PCA
    fts_object = fts(1:1621, data_return[,1:(jk+84)])                    
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "classical")$varprop)>=0.9),1)
    for(i in 2:(length(trading_time)-1))
    {
      # mse
      dum = dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(jk+85)],
                        method = "pls", lambda = pls_find_lambda_mse[ceiling((i-1)/10)], order = ftsm_object_order, pcdmethod = "classical", interval = FALSE)
      pls_point_update_mse[i-1] = dum$errormse
      
      # mae
      dum_mae = dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(jk+85)],
                            method = "pls", lambda = pls_find_lambda_mae[ceiling((i-1)/10)], order = ftsm_object_order, pcdmethod = "classical", interval = FALSE)
      pls_point_update_mae[i-1] = dum$errormae
    }
    return(list(pls_point_update_mae = pls_point_update_mae, pls_point_update_mse = pls_point_update_mse))
  }
  
  if(pcdmethod == "M")
  {
    # creating storage for errors
    pls_point_update_M_mse = vector(, (length(trading_time)-2))
    pls_point_update_M_mae = vector(, (length(trading_time)-2))
    
    # doing functional PCA
    fts_object = fts(1:1621, fts_return[,1:(jk+84)])
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "M", lambda = 2.33)$varprop)>=0.9),1)
    for(i in 2:(length(trading_time)-1))
    {
      # mse
      dum = dynupdate(data = fts_object, newdata = fts_return[1:i,(jk+85)], holdoutdata = fts_return[(i+1):1621,(jk+85)],
                        method = "pls", pcdmethod = "M", lambda = pls_find_lambda_M_mse[ceiling((i-1)/10)], order = ftsm_object_order, interval = FALSE)
      pls_point_update_M_mse[i-1] = dum$errormse
      
      # mae
      dum_mae = dynupdate(data = fts_object, newdata = fts_return[1:i,(jk+85)], holdoutdata = fts_return[(i+1):1621,(jk+85)],
                            method = "pls",  pcdmethod = "M", lambda = pls_find_lambda_M_mae[ceiling((i-1)/10)], order = ftsm_object_order, interval = FALSE)
      pls_point_update_M_mae[i-1] = dum$errormae
    }
    return(list(pls_point_update_M_mae = pls_point_update_M_mae, pls_point_update_M_mse = pls_point_update_M_mse))
  }
}

# selecting optimal lambdas, pcdmethod = "classical"

library(doParallel)
cl <- makeCluster(10) 
registerDoParallel(cl)

####
htm = proc.time()
select_lambda_ind = seq(2,length(trading_time)-1, 10)
w_classical = foreach(i = select_lambda_ind, .packages = c("ftsa")) %dopar% optimise(pls_find_lambda, interval = c(0, 10000), pcdmethod = "classical", ik = i, error = "mse")
htm_pls_lambda_mse = proc.time() - htm

htm = proc.time()
select_lambda_ind = seq(2,length(trading_time)-1, 10)
w_classical_mae = foreach(i = select_lambda_ind, .packages = c("ftsa")) %dopar% optimise(pls_find_lambda, interval = c(0, 10000), pcdmethod = "classical", ik = i, error = "mae")
htm_pls_lambda_mae = proc.time() - htm

###

htm = proc.time()
w_classical = foreach(i = 2:(length(trading_time)-1), .packages = c("ftsa")) %dopar% optimise(pls_find_lambda, interval = c(0, 10000), pcdmethod = "classical", ik = i, error = "mse")
htm_pls_lambda_mse = proc.time() - htm

htm = proc.time()
w_classical_mae = foreach(i = 2:(length(trading_time)-1), .packages = c("ftsa")) %dopar% optimise(pls_find_lambda, interval = c(0, 10000), pcdmethod = "classical", ik = i, error = "mae")
htm_pls_lambda_mae = proc.time() - htm

pls_find_lambda_mse = pls_find_objective_mse = vector(,(length(trading_time)-2))
for(i in 1:(length(trading_time)-2))
{
  pls_find_lambda_mse[i] = w_classical[[i]]$minimum
  pls_find_objective_mse[i] = w_classical[[i]]$objective
}

pls_find_lambda_mae = pls_find_objective_mae = vector(,(length(trading_time)-2))
for(i in 1:(length(trading_time)-2))
{
  pls_find_lambda_mae[i] = w_classical_mae[[i]]$minimum
  pls_find_objective_mae[i] = w_classical_mae[[i]]$objective
}


# updating errors, pcdmethod = "classical"

htm = proc.time()

pls_point_update_err_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% PLS_update(fts_return, "classical", jk = j)

# record errors
pls_point_update_mae = pls_point_update_mse = matrix(, 1619, 40)
for (j in 1:40)
{
  pls_point_update_mae[,j] = pls_point_update_err_list[[j]]$pls_point_update_mae
  pls_point_update_mse[,j] = pls_point_update_err_list[[j]]$pls_point_update_mse
}

pls_point_update_err = cbind(round(rowMeans(pls_point_update_mse), 4), round(rowMeans(pls_point_update_mae), 4))
colnames(pls_point_update_err) = c("mse", "mae")

htm_PLS = proc.time() - htm


# selecting optimal lambdas, pcdmethod = "M"

library(doParallel)
cl <- makeCluster(10) 
registerDoParallel(cl)

htm = proc.time()
select_lambda_ind = seq(2,length(trading_time)-1, 10)
w_M = foreach(i = select_lambda_ind, .packages = c("ftsa")) %dopar% optimise(pls_find_lambda, interval = c(0, 10000), pcdmethod = "M", ik = i, error = "mse")
htm_pls_lambda_M_mse = proc.time() - htm

htm = proc.time()
select_lambda_ind = seq(2,length(trading_time)-1, 10)
w_M_mae = foreach(i = select_lambda_ind, .packages = c("ftsa")) %dopar% optimise(pls_find_lambda, interval = c(0, 10000), pcdmethod = "M", ik = i, error = "mae")
htm_pls_lambda_M_mae = proc.time() - htm

pls_find_lambda_M_mse = pls_find_objective_M_mse = vector(,length(select_lambda_ind))
for(i in 1:length(select_lambda_ind))
{
  pls_find_lambda_M_mse[i] = w_M[[i]]$minimum
  pls_find_objective_M_mse[i] = w_M[[i]]$objective
}

pls_find_lambda_M_mae = pls_find_objective_M_mae = vector(,length(select_lambda_ind))
for(i in 1:length(select_lambda_ind))
{
  pls_find_lambda_M_mae[i] = w_M_mae[[i]]$minimum
  pls_find_objective_M_mae[i] = w_M_mae[[i]]$objective
}


# updating errors, pcdmethod = "M"

htm = proc.time()
library(doParallel)
cl <- makeCluster(10) 
registerDoParallel(cl)
pls_point_update_M_err_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% PLS_update(fts_return, "M", jk = j)

# record errors
pls_point_update_M_mae = pls_point_update_M_mse = matrix(, 1619, 40)
for (j in 1:40)
{
  pls_point_update_M_mae[,j] = pls_point_update_M_err_list[[j]]$pls_point_update_M_mae
  pls_point_update_M_mse[,j] = pls_point_update_M_err_list[[j]]$pls_point_update_M_mse
}

pls_point_update_M_err = cbind(round(rowMeans(pls_point_update_M_mse), 4), round(rowMeans(pls_point_update_M_mae), 4))
colnames(pls_point_update_M_err) = c("mse", "mae")

htm_PLS_M = proc.time() - htm

##############
# RR updating
##############

# function for selecting lambda

rr_find_lambda <- function(find_lambda, pcdmethod = c("classical","M"), error = c("mse", "mae", "score"), ik)
{
    pcdmethod = match.arg(pcdmethod)
    error = match.arg(error)
    if(error == "mse")
    {
      rr_point_update_mse_train = vector("numeric",42)
      for(j in 1:42)
      {
        fts_object = fts(1:1621, fts_return[,1:(j+42)])
        ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = pcdmethod, lambda=2.33)$varprop)>=0.9),1)
        dum = dynupdate(data = fts_object, newdata = fts_return[1:ik,(j+43)], holdoutdata = fts_return[(ik+1):1621,(j+43)],  	
                        method = "ridge", pcdmethod = pcdmethod, order = ftsm_object_order, interval = FALSE, lambda = find_lambda) 
        rr_point_update_mse_train[j] = dum$errormse
      }
      return(mean(rr_point_update_mse_train))
    }
    if(error == "mae")
    {
      rr_point_update_mae_train = vector("numeric",42)
      for(j in 1:42)
      {
        fts_object = fts(1:1621, fts_return[,1:(j+42)])
        ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = pcdmethod, lambda=2.33)$varprop)>=0.9),1)
        dum = dynupdate(data = fts_object, newdata = fts_return[1:ik,(j+43)], holdoutdata = fts_return[(ik+1):1621,(j+43)],  	
                        method = "ridge", pcdmethod = pcdmethod, order = ftsm_object_order, interval = FALSE, lambda = find_lambda) 
        rr_point_update_mae_train[j] = dum$errormae
      }
      return(median(rr_point_update_mae_train))
    }
}

# function for doing RR updating

RR_update = function(data_return, pcdmethod = c("classical","M"), jk)
{
    pcdmethod = match.arg(pcdmethod)
    if(pcdmethod == "classical")
    {
      # creating storage for errors
      rr_point_update_mse = vector("numeric", (length(trading_time)-2))
      rr_point_update_mae = vector("numeric", (length(trading_time)-2))
      
      # doing functional PCA
      fts_object = fts(1:1621, data_return[,1:(jk+84)])
      ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "classical")$varprop)>=0.9),1)
      for(i in 2:(length(trading_time)-1))
      {
        # mse
        dum = dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(jk+85)],
                        method = "ridge", lambda = rr_find_lambda_mse[ceiling((i-1)/10)], order = ftsm_object_order, pcdmethod = pcdmethod, interval = FALSE)
        rr_point_update_mse[i-1] = dum$errormse
        
        # mae
        dum_mae = dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(jk+85)],
                            method = "ridge", lambda = rr_find_lambda_mae[ceiling((i-1)/10)], order = ftsm_object_order, pcdmethod = pcdmethod, interval = FALSE)
        rr_point_update_mae[i-1] = dum$errormae
      }
      return(list(rr_point_update_mae = rr_point_update_mae, rr_point_update_mse = rr_point_update_mse))
    }
    
    if(pcdmethod == "M")
    {
      # creating storage for errors
      rr_point_update_M_mse = vector("numeric", (length(trading_time)-2))
      rr_point_update_M_mae = vector("numeric", (length(trading_time)-2))
      
      # doing functional PCA
      fts_object = fts(1:1621, fts_return[,1:(jk+84)])
      ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "M", lambda = 2.33)$varprop)>=0.9),1)
      for(i in 2:(length(trading_time)-1))
      {
        # mse
        dum = dynupdate(data = fts_object, newdata = fts_return[1:i,(jk+85)], holdoutdata = fts_return[(i+1):1621,(jk+85)],
                        method = "ridge", pcdmethod = pcdmethod, lambda = rr_find_lambda_M_mse[ceiling((i-1)/10)], order = ftsm_object_order, interval = FALSE)
        rr_point_update_M_mse[i-1] = dum$errormse
        
        # mae
        dum_mae = dynupdate(data = fts_object, newdata = fts_return[1:i,(jk+85)], holdoutdata = fts_return[(i+1):1621,(jk+85)],
                            method = "ridge",  pcdmethod = pcdmethod, lambda = rr_find_lambda_M_mae[ceiling((i-1)/10)], order = ftsm_object_order, interval = FALSE)
        rr_point_update_M_mae[i-1] = dum$errormae
      }
      return(list(rr_point_update_M_mae = rr_point_update_M_mae, rr_point_update_M_mse = rr_point_update_M_mse))
    }
}

# selecting optimal lambdas, pcdmethod = "classical"


htm = proc.time()
cl <- makeCluster(10) 
registerDoParallel(cl)

select_lambda_ind = seq(2,length(trading_time)-1, 10)
rr_classical = foreach(i = select_lambda_ind, .packages = c("ftsa")) %dopar% optimise(rr_find_lambda, interval = c(0, 10000), pcdmethod = "classical", ik = i, error = "mse")
htm_rr_lambda_mse = proc.time() - htm

save(rr_classical, file = "rr_classical.RData")

htm = proc.time()
cl <- makeCluster(10) 
registerDoParallel(cl)

select_lambda_ind = seq(2,length(trading_time)-1, 10)
rr_classical_mae = foreach(i = select_lambda_ind, .packages = c("ftsa")) %dopar% optimise(rr_find_lambda, interval = c(0, 10000), pcdmethod = "classical", ik = i, error = "mae")
htm_rr_lambda_mae = proc.time() - htm

save(rr_classical_mae, file = "rr_classical_mae.RData")

rr_find_lambda_mse = rr_find_objective_mse = vector("numeric",length(select_lambda_ind))
for(i in 1:length(select_lambda_ind))
{
    rr_find_lambda_mse[i] = rr_classical[[i]]$minimum
    rr_find_objective_mse[i] = rr_classical[[i]]$objective
}

rr_find_lambda_mae = rr_find_objective_mae = vector("numeric",length(select_lambda_ind))
for(i in 1:length(select_lambda_ind))
{
    rr_find_lambda_mae[i] = rr_classical_mae[[i]]$minimum
    rr_find_objective_mae[i] = rr_classical_mae[[i]]$objective
}

# updating errors, pcdmethod = "classical"

htm = proc.time()

library(doParallel)
cl <- makeCluster(10) 
registerDoParallel(cl)

rr_point_update_err_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% RR_update(fts_return, "classical", jk = j)

# recording errors
rr_point_update_mae = rr_point_update_mse = matrix(NA, 1619, 40)
for(j in 1:40)
{
    rr_point_update_mae[,j] = rr_point_update_err_list[[j]]$rr_point_update_mae
    rr_point_update_mse[,j] = rr_point_update_err_list[[j]]$rr_point_update_mse
}

rr_point_update_err = cbind(round(rowMeans(rr_point_update_mse), 4), round(rowMeans(rr_point_update_mae), 4))
colnames(rr_point_update_err) = c("mse", "mae")
htm_RR = proc.time() - htm

# selecting optimal lambdas, pcdmethod = "M"

htm = proc.time()
cl <- makeCluster(10) 
registerDoParallel(cl)

select_lambda_ind = seq(2,length(trading_time)-1, 10)
rr_M = foreach(i = select_lambda_ind, .packages = c("ftsa")) %dopar% optimise(rr_find_lambda, interval = c(0, 10000), pcdmethod = "M", ik = i, error = "mse")
htm_rr_lambda_M_mse = proc.time() - htm

save(rr_M, file = "rr_M.RData")

htm = proc.time()
cl <- makeCluster(10) 
registerDoParallel(cl)

select_lambda_ind = seq(2,length(trading_time)-1, 10)
rr_M_mae = foreach(i = select_lambda_ind, .packages = c("ftsa")) %dopar% optimise(rr_find_lambda, interval = c(0, 10000), pcdmethod = "M", ik = i, error = "mae")
htm_rr_lambda_M_mae = proc.time() - htm

save(rr_M_mae, file = "rr_M_mae.RData")

rr_find_lambda_M_mse = rr_find_objective_M_mse = vector("numeric",length(select_lambda_ind))
for(i in 1:length(select_lambda_ind))
{
    rr_find_lambda_M_mse[i] = rr_M[[i]]$minimum
    rr_find_objective_M_mse[i] = rr_M[[i]]$objective
}

rr_find_lambda_M_mae = rr_find_objective_M_mae = vector("numeric",length(select_lambda_ind))
for(i in 1:length(select_lambda_ind))
{
    rr_find_lambda_M_mae[i] = rr_M_mae[[i]]$minimum
    rr_find_objective_M_mae[i] = rr_M_mae[[i]]$objective
}

# updating errors, pcdmethod = "M"

htm = proc.time()

library(doParallel)
cl <- makeCluster(10) 
registerDoParallel(cl)

rr_point_update_M_err_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% RR_update(fts_return, "M", jk = j)

# record errors
rr_point_update_M_mae = rr_point_update_M_mse = matrix(NA, 1619, 40)
for(j in 1:40)
{
    rr_point_update_M_mae[,j] = rr_point_update_M_err_list[[j]]$rr_point_update_M_mae
    rr_point_update_M_mse[,j] = rr_point_update_M_err_list[[j]]$rr_point_update_M_mse
}

rr_point_update_M_err = cbind(round(rowMeans(rr_point_update_M_mse), 4), round(rowMeans(rr_point_update_M_mae), 4))
colnames(rr_point_update_M_err) = c("mse", "mae")
htm_RR_M = proc.time() - htm

#################
# Chiou's method 
#################

Chiou_pca_update <- function(data_return, pcd_method = c("classical", "M"), jk, robust_lambda = 2.33)
{
    pcd_method = match.arg(pcd_method)
    if(pcd_method == "classical")
    {
        # creating storage for errors
        Chiou_point_update_mse = vector("numeric", length(trading_time))
        Chiou_point_update_mae = vector("numeric", length(trading_time))
        cor_sign_prob = vector("numeric", length(trading_time))
        # updating errors
        for(i in 6:1620)
        {
            dum = dynamic_Chiou_pca(dat = data_return[,1:(jk+84)], newdata = matrix(data_return[1:i, (jk+85)], ncol=i), 
                                holdoutdata = matrix(data_return[(i+1):1621,(jk+85)], ncol=1621-i),
                                order_k_percent = 0.9, order_m_percent = 0.9, pointfore = TRUE)
            neg_fore = sign(as.numeric(dum$update_forecast))
            neg_true = sign(as.numeric(dum$holdoutdata))
            cor_sign_prob[i-1] = length(which(neg_fore == neg_true))/length(neg_true)		                            
            Chiou_point_update_mse[i-1] = dum$err[,2]
            Chiou_point_update_mae[i-1] = dum$err[,1]
        }
        return(list(Chiou_point_update_mse = Chiou_point_update_mse, Chiou_point_update_mae = Chiou_point_update_mae, cor_sign_prob = cor_sign_prob))
    }
    if(pcd_method == "M")
    {
        # creating storage for errors
        Chiou_point_update_M_mse = vector("numeric", length(trading_time))
        Chiou_point_update_M_mae = vector("numeric", length(trading_time))
        # updating errors
        for(i in 6:1620)
        {
            dum = dynamic_Chiou_pca(dat = data_return[,1:(jk+84)], newdata = matrix(data_return[1:i, (jk+85)], ncol=i), 
                                holdoutdata = matrix(data_return[(i+1):1621,(jk+85)], ncol=1621-i),
                                order_k_percent = 0.9, order_m_percent = 0.9, pcd_method = "M", robust_lambda = robust_lambda, pointfore = TRUE)
            Chiou_point_update_M_mse[i-1] = dum$err[,2]
            Chiou_point_update_M_mae[i-1] = dum$err[,1]
        }
        return(list(Chiou_point_update_M_mse = Chiou_point_update_M_mse, Chiou_point_update_M_mae = Chiou_point_update_M_mae))
    }
}


# pcdmethod = "classical"
htm = proc.time()
library(doParallel)
cl <- makeCluster(10) 
registerDoParallel(cl)


Chiou_point_update_err_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% Chiou_pca_update(fts_return, pcd_method = "classical", jk = j)

# record errors
Chiou_point_update_mae = Chiou_point_update_mse = cor_sign_prob = matrix(, 1621, 40)
for (j in 1:40)
{
  Chiou_point_update_mae[,j] = Chiou_point_update_err_list[[j]]$Chiou_point_update_mae
  Chiou_point_update_mse[,j] = Chiou_point_update_err_list[[j]]$Chiou_point_update_mse
  cor_sign_prob[,j] = Chiou_point_update_err_list[[j]]$cor_sign_prob
}

Chiou_point_update_err = cbind(round(rowMeans(Chiou_point_update_mse), 4), round(rowMeans(Chiou_point_update_mae), 4))
colnames(Chiou_point_update_err) = c("mse", "mae")
htm_Chiou = proc.time() - htm

# pcdmethod = "M"
htm = proc.time()
library(doParallel)
cl <- makeCluster(10) 
registerDoParallel(cl)


Chiou_point_update_M_err_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% Chiou_pca_update(fts_return, "M", jk = j)

# record errors
Chiou_point_update_M_mae = Chiou_point_update_M_mse = matrix(, 1621, 40)
for (j in 1:40)
{
  Chiou_point_update_M_mae[,j] = Chiou_point_update_M_err_list[[j]]$Chiou_point_update_M_mae
  Chiou_point_update_M_mse[,j] = Chiou_point_update_M_err_list[[j]]$Chiou_point_update_M_mse
}

Chiou_point_update_M_err = cbind(round(rowMeans(Chiou_point_update_M_mse), 4), round(rowMeans(Chiou_point_update_M_mae), 4))
colnames(Chiou_point_update_M_err) = c("mse", "mae")
htm_Chiou_M = proc.time() - htm

# pcdmethod = "M", lambda = 3
htm = proc.time()
library(doParallel)
cl <- makeCluster(10) 
registerDoParallel(cl)


Chiou_point_update_M_3_err_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% Chiou_pca_update(fts_return, "M", jk = j, robust_lambda = 3)

# record errors
Chiou_point_update_M_3_mae = Chiou_point_update_M_3_mse = matrix("numeric", 1621, 40)
for (j in 1:40)
{
  Chiou_point_update_M_3_mae[,j] = Chiou_point_update_M_3_err_list[[j]]$Chiou_point_update_M_mae
  Chiou_point_update_M_3_mse[,j] = Chiou_point_update_M_3_err_list[[j]]$Chiou_point_update_M_mse
}

Chiou_point_update_M_3_err = cbind(round(rowMeans(Chiou_point_update_M_3_mse), 4), round(rowMeans(Chiou_point_update_M_3_mae), 4))
colnames(Chiou_point_update_M_3_err) = c("mse", "mae")
htm_Chiou_M_3

# pcdmethod = "M", lambda = 1.81
htm = proc.time()
library(doParallel)
cl <- makeCluster(10) 
registerDoParallel(cl)

Chiou_point_update_M_1.81_err_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% Chiou_pca_update(fts_return, "M", jk = j, robust_lambda = 1.81)

# record errors
Chiou_point_update_M_1.81_mae = Chiou_point_update_M_1.81_mse = matrix("numeric", 1621, 40)
for (j in 1:40)
{
  Chiou_point_update_M_1.81_mae[,j] = Chiou_point_update_M_1.81_err_list[[j]]$Chiou_point_update_M_mae
  Chiou_point_update_M_1.81_mse[,j] = Chiou_point_update_M_1.81_err_list[[j]]$Chiou_point_update_M_mse
}

Chiou_point_update_M_1.81_err = cbind(round(rowMeans(Chiou_point_update_M_1.81_mse), 4), round(rowMeans(Chiou_point_update_M_1.81_mae), 4))
colnames(Chiou_point_update_M_1.81_err) = c("mse", "mae")
htm_Chiou_M_1.81

# summarising point forecast results 

point_update_mse_all = cbind(ts_point_update_err[,1], bm_point_update_err[,1], pls_point_update_err[,1], rr_point_update_err[,1], Chiou_point_update_err[-(1:2),1])
point_update_mae_all = cbind(ts_point_update_err[,2], bm_point_update_err[,2], pls_point_update_err[,2], rr_point_update_err[,2], Chiou_point_update_err[-(1:2),2])
colnames(point_update_mse_all) = colnames(point_update_mae_all) = c("TS", "BM", "PLS", "RR", "FLR")

trading_time2 = seq(9.5, 16.25, by = 15/3600)[-c(1,2)]

# plots without realisations

savepdf("dynamic_update_1_continuous",width=12,height=10,toplines=0.9,pointsize=10)
plot(trading_time2, point_update_mse_all[,1], type="l", ylim = range(point_update_mse_all) , xlab="Five-minute time interval", ylab="MSFE")

lines(trading_time2, point_update_mse_all[,2], col=2)

lines(trading_time2, point_update_mse_all[,3], col=3)

lines(trading_time2, point_update_mse_all[,4], col=4)

lines(trading_time2, point_update_mse_all[,5], col=5)
legend(10,10,c("TS","BM","PLS","RR","FLR"),cex=0.6,col=1:5, lty=rep(1,5))
dev.off()

savepdf("dynamic_update_2_continuous",width=12,height=10,toplines=0.9,pointsize=10)
plot(trading_time2, point_update_mae_all[,1], type = "l", ylim = range(point_update_mae_all), xlab="Five-minute time interval", ylab="MAFE", cex=0.5)

lines(trading_time2, point_update_mae_all[,2], col=2)

lines(trading_time2, point_update_mae_all[,3], col=3)

lines(trading_time2, point_update_mae_all[,4], col=4)

lines(trading_time2, point_update_mae_all[,5], col=5)

legend(10,1.5,c("TS","BM","PLS","RR","FLR"),cex=0.6,col=1:5, lty=rep(1,5))
dev.off()

# plots with realisations

savepdf("dynamic_update_1",width=12,height=10,toplines=0.9,pointsize=10)
plot(trading_time2, point_update_mse_all[,1], type="l", ylim = range(point_update_mse_all) , xlab="Time of Day", ylab="MSFE")
points(trading_time2, point_update_mse_all[,1], col=1, cex=0.25)

lines(trading_time2, point_update_mse_all[,2], col=2)
points(trading_time2, point_update_mse_all[,2], col=2, pch=2, cex=0.25)

lines(trading_time2, point_update_mse_all[,3], col=3)
points(trading_time2, point_update_mse_all[,3], col=3, pch=3, cex=0.25)

lines(trading_time2, point_update_mse_all[,4], col=4)
points(trading_time2, point_update_mse_all[,4], col=4, pch=4, cex=0.25)

lines(trading_time2[3:1617], point_update_mse_all[3:1617,5], col=5)
points(trading_time2[3:1617], point_update_mse_all[3:1617,5], col=5, pch=5, cex=0.25)
legend("bottomleft",c("TS","BM","PLS","RR","FLR"),cex=0.6,pch =1:5,col=1:5)
dev.off()

savepdf("dynamic_update_2",width=12,height=10,toplines=0.9,pointsize=10)
plot(trading_time2, point_update_mae_all[,1], type = "l", ylim = range(point_update_mae_all), xlab="Time of Day", ylab="MAFE", cex=0.5)
points(trading_time2, point_update_mae_all[,1], col=1, cex = 0.25)

lines(trading_time2, point_update_mae_all[,2], col=2)
points(trading_time2, point_update_mae_all[,2], col=2, pch=2, cex=0.25)

lines(trading_time2, point_update_mae_all[,3], col=3)
points(trading_time2, point_update_mae_all[,3], col=3, pch=3, cex=0.25)

lines(trading_time2, point_update_mae_all[,4], col=4)
points(trading_time2, point_update_mae_all[,4], col=4, pch=4, cex=0.25)

lines(trading_time2[3:1617], point_update_mae_all[3:1617,5], col=5)
points(trading_time2[3:1617], point_update_mae_all[3:1617,5], col=5, pch=5, cex=0.25)
legend("bottomleft",c("TS","BM","PLS","RR","FLR"),cex=0.6,pch = 1:5, col = 1:5)
dev.off()







