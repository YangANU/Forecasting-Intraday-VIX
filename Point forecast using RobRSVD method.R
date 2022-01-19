##############################
# Chiou's method with RobRSVD 
##############################

Chiou_rsvd_update <- function(data_return, jk, robust_lambda = 2.33)
{
    
    
    # creating storage for errors
    Chiou_point_rsvd_update_mse = vector("numeric", length(trading_time))
    Chiou_point_rsvd_update_mae = vector("numeric", length(trading_time))
    # updating errors
    for(i in 10:1618)
    {
        dum = dynamic_Chiou_robrsvd(dat = data_return[,1:(jk+84)], newdata = matrix(data_return[1:i, (jk+85)], ncol=i), 
                                holdoutdata = matrix(data_return[(i+1):1621,(jk+85)], ncol=1621-i),
                                order_k_percent = 0.9, order_m_percent = 0.9, robust_lambda = robust_lambda, pointfore = TRUE)
        Chiou_point_rsvd_update_mse[i-1] = dum$err[,2]
        Chiou_point_rsvd_update_mae[i-1] = dum$err[,1]
    }
    return(list(Chiou_point_rsvd_updatemse = Chiou_point_rsvd_update_mse, Chiou_point_rsvd_update_mae = Chiou_point_rsvd_update_mae))
}


htm = proc.time()
library(doParallel)
cl <- makeCluster(30) 
registerDoParallel(cl)


Chiou_point_rsvd_update_err_list = foreach(j = 1:40, .packages = c("ftsa", "RobRSVD")) %dopar% Chiou_rsvd_update(fts_return, jk = j)

# record errors
Chiou_point_rsvd_update_mae = Chiou_point_rsvd_update_mse = matrix(, 1621, 40)
for (j in 1:40)
{
  Chiou_point_rsvd_update_mae[,j] = Chiou_point_rsvd_update_err_list[[j]]$Chiou_point_rsvd_update_mae
  Chiou_point_rsvd_update_mse[,j] = Chiou_point_rsvd_update_err_list[[j]]$Chiou_point_rsvd_update_mse
}

Chiou_point_rsvd_update_err = cbind(round(rowMeans(Chiou_point_rsvd_update_mse), 4), round(rowMeans(Chiou_point_rsvd_update_mae), 4))
colnames(Chiou_point_rsvd_update_err) = c("mse", "mae")
htm_Chiou = proc.time() - htm




