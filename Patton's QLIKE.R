#############################################################################
# Evaluate point forecasts using the Patton's quasi-likelihood loss function
#############################################################################
library(doParallel)

# QLIKE: L(X_hat, X) = log(X_hat) + X/X_hat

# TS (non-updating)

TS_noupdate_forecasts <- function(data_return, pcdmethod = c("classical","M"), jk)
{
  pcdmethod = match.arg(pcdmethod)
  if(pcdmethod == "classical")
  {
    # creating storage for errors
    ts_noupdate_fore = list()
    
    # doing functional PCA
    fts_object = fts(1:1621, data_return[,1:(jk+84)])
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = pcdmethod, lambda = 2.33)$varprop)>=0.9),1)
    for(i in 1:(length(trading_time)-2))
    {
      ts_noupdate_fore[[i]] = dynupdate(data = fts_object, newdata = data_return[1:(i+1),(jk+85)], holdoutdata = data_return[(i+2):1621,(jk+85)],  	
                      method = "ts", pcdmethod = pcdmethod, order = ftsm_object_order, interval = FALSE, value = TRUE) 
    }
    return(list(ts_noupdate_fore = ts_noupdate_fore))
  }
}

cl <- makeCluster(30) 
registerDoParallel(cl)

ts_point_update_err_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% TS_noupdate_forecasts(fts_return, "classical", jk = j)





# Chiou's method 

Chiou_pca_forecasts <- function(data_return, pcd_method = c("classical", "M"), jk, robust_lambda = 2.33)
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
}















