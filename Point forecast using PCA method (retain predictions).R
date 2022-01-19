####################
# TS (non-updating)
####################

TS_noupdate_forecasts <- function(data_return, pcdmethod = c("classical","M"), jk)
{
  pcdmethod = match.arg(pcdmethod)
  if(pcdmethod == "classical")
  {
    # creating storage for forecasts
    ts_point_update = list()
    
    # doing functional PCA
    fts_object = fts(1:1621, data_return[,1:(jk+84)])
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = pcdmethod, lambda = 2.33)$varprop)>=0.9),1)
    for(i in 2:(length(trading_time)-2))
    {
      ts_point_update[[i]] = dynupdate(data = fts_object, newdata = data_return[1:(i+1),(jk+85)], holdoutdata = data_return[(i+2):1621,(jk+85)],  	
                      method = "ts", pcdmethod = pcdmethod, order = ftsm_object_order, interval = FALSE, value = TRUE) 

    }
    return(ts_point_update)
  }
  if(pcdmethod == "M")
  {
    # creating storage for forecasts
    ts_point_update_M = list()
    
    # doing functional PCA
    fts_object = fts(1:1621, data_return[,1:(jk+84)])
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = pcdmethod, lambda = 2.33)$varprop)>=0.9),1)
    for(i in 1:(length(trading_time)-2))
    {
      ts_point_update_M[[i]] = dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(jk+85)],  	
                      method = "ts", pcdmethod = pcdmethod, order = ftsm_object_order, robust_lambda = 2.33, interval = FALSE, value = TRUE) 
    }
    return(ts_point_update_M)
  }
}

# pcdmethod = "classical"
htm = proc.time()

library(doParallel)
cl <- makeCluster(30) 
registerDoParallel(cl)

ts_point_update_forecasts_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% TS_noupdate_forecasts(fts_return, "classical", jk = j)

save(ts_point_update_forecasts_list, file = "ts_point_update_forecasts.RData")

# pcdmethod = "M"
ts_point_update_M_forecasts_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% TS_noupdate_forecasts(fts_return, "M", jk = j)

save(ts_point_update_M_forecasts_list, file = "ts_point_update_M_forecasts.RData")
proc.time() - htm

##############
# BM updating
##############

htm = proc.time()

BM_update_forecasts <- function(data_return, pcdmethod = c("classical","M"), jk)
{
  pcdmethod = match.arg(pcdmethod)
  if(pcdmethod == "classical")
  {
    # creating storage for forecasts
    bm_point_update = list()
    
    # doing functional PCA
    fts_object = fts(1:1621, data_return[,1:(jk+84)])
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "classical")$varprop)>=0.9),1)
    for(i in 1:(length(trading_time)-2))
    {
      bm_point_update[[i]] = dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(j+85)],  	
                      method = "block", pcdmethod = "classical", order = ftsm_object_order, interval = FALSE, value = TRUE) 
    }
    return(bm_point_update = bm_point_update)
  }
  
  if(pcdmethod == "M")
  {
    # creating storage for forecasts
    bm_point_update_M = list()
    
    # doing functional PCA
    fts_object = fts(1:1621, data_return[,1:(jk+84)])
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "M", lambda=2.33)$varprop)>=0.9),1)
    for(i in 1:(length(trading_time)-2))
    {
      bm_point_update_M[[i]] = dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(jk+85)],  	
                      method = "block", pcdmethod = "M", order = ftsm_object_order, robust_lambda = 2.33, interval = FALSE, value = TRUE)
    }
    return(bm_point_update_M = bm_point_update_M)
  }
}

# pcdmethod = "classical"

cl <- makeCluster(30) 
registerDoParallel(cl)

bm_point_update_forecasts_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% BM_update_forecasts(fts_return, "classical", jk = j)

save(bm_point_update_forecasts_list, file = "bm_point_update_forecasts.RData")
# pcdmethod = "M"

cl <- makeCluster(30) 
registerDoParallel(cl)

bm_point_update_M_forecasts_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% BM_update_forecasts(fts_return, "M", jk = j)

save(bm_point_update_M_forecasts_list, file = "bm_point_update_M_forecasts.RData")
proc.time() - htm

###############
# PLS updating
###############

# function for doing PLS updating

PLS_update_forecasts = function(data_return, pcdmethod = c("classical","M"), jk)
{
  pcdmethod = match.arg(pcdmethod)
  if(pcdmethod == "classical")
  {
    # creating storage for forecasts
    pls_point_update_mse_forecasts = list()
    pls_point_update_mae_forecasts = list()
    
    # doing functional PCA
    fts_object = fts(1:1621, data_return[,1:(jk+84)])                    
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "classical")$varprop)>=0.9),1)
    for(i in 2:(length(trading_time)-1))
    {
      pls_point_update_mse_forecasts[[i]] = dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(jk+85)],
                      method = "pls", lambda = pls_find_lambda_mse[ceiling((i-1)/10)], order = ftsm_object_order, pcdmethod = "classical", interval = FALSE, value = TRUE)
      
      pls_point_update_mae_forecasts[[i]] = dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(jk+85)],
                          method = "pls", lambda = pls_find_lambda_mae[ceiling((i-1)/10)], order = ftsm_object_order, pcdmethod = "classical", interval = FALSE, value = TRUE)

    }
    return(list(pls_point_update_mae_forecasts = pls_point_update_mae_forecasts, pls_point_update_mse_forecasts = pls_point_update_mse_forecasts))
  }
  
  if(pcdmethod == "M")
  {
    # creating storage for forecasts
    pls_point_update_M_mse_forecasts = list()
    pls_point_update_M_mae_forecasts = list()
    
    # doing functional PCA
    fts_object = fts(1:1621, fts_return[,1:(jk+84)])
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "M", lambda = 2.33)$varprop)>=0.9),1)
    for(i in 2:(length(trading_time)-1))
    {
      pls_point_update_M_mse_forecasts[[i]] = dynupdate(data = fts_object, newdata = fts_return[1:i,(jk+85)], holdoutdata = fts_return[(i+1):1621,(jk+85)],
                      method = "pls", pcdmethod = "M", lambda = pls_find_lambda_M_mse[ceiling((i-1)/10)], order = ftsm_object_order, interval = FALSE, value = TRUE)
      
      pls_point_update_M_mae_forecasts[[i]] = dynupdate(data = fts_object, newdata = fts_return[1:i,(jk+85)], holdoutdata = fts_return[(i+1):1621,(jk+85)],
                          method = "pls",  pcdmethod = "M", lambda = pls_find_lambda_M_mae[ceiling((i-1)/10)], order = ftsm_object_order, interval = FALSE, value = TRUE)
      
    }
    return(list(pls_point_update_M_mae_forecasts = pls_point_update_M_mae_forecasts, pls_point_update_M_mse_forecasts = pls_point_update_M_mse_forecasts))
  }
}

# updating errors, pcdmethod = "classical"

htm = proc.time()

pls_point_update_forecasts_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% PLS_update_forecasts(fts_return, "classical", jk = j)

save(pls_point_update_forecasts_list, file = "pls_point_update_forecasts.RData")

htm_PLS = proc.time() - htm

# updating errors, pcdmethod = "M"

htm = proc.time()

cl <- makeCluster(30) 
registerDoParallel(cl)
pls_point_update_M_forecasts_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% PLS_update_forecasts(fts_return, "M", jk = j)

save(pls_point_update_M_forecasts_list, file = "pls_point_update_M_forecasts.RData")

htm_PLS_M = proc.time() - htm

##############
# RR updating
##############


# function for doing RR updating

RR_update_forecasts = function(data_return, pcdmethod = c("classical","M"), jk)
{
  pcdmethod = match.arg(pcdmethod)
  if(pcdmethod == "classical")
  {
    # creating storage for errors
    rr_point_update_mse_forecasts = list()
    rr_point_update_mae_forecasts = list()
    
    # doing functional PCA
    fts_object = fts(1:1621, data_return[,1:(jk+84)])
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "classical")$varprop)>=0.9),1)
    for(i in 2:(length(trading_time)-1))
    {
      # mse
      rr_point_update_mse_forecasts[[i]]= dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(jk+85)],
                      method = "ridge", lambda = rr_find_lambda_mse[ceiling((i-1)/10)], order = ftsm_object_order, pcdmethod = pcdmethod, interval = FALSE, value = TRUE)
      
      # mae
      rr_point_update_mae_forecasts[[i]] = dynupdate(data = fts_object, newdata = data_return[1:i,(jk+85)], holdoutdata = data_return[(i+1):1621,(jk+85)],
                          method = "ridge", lambda = rr_find_lambda_mae[ceiling((i-1)/10)], order = ftsm_object_order, pcdmethod = pcdmethod, interval = FALSE, value = TRUE)
    }
    return(list(rr_point_update_mae_forecasts = rr_point_update_mae_forecasts, rr_point_update_mse_forecasts = rr_point_update_mse_forecasts))
  }
  
  if(pcdmethod == "M")
  {
    # creating storage for errors
    rr_point_update_M_mse_forecasts = list()
    rr_point_update_M_mae_forecasts = list()
    
    # doing functional PCA
    fts_object = fts(1:1621, fts_return[,1:(jk+84)])
    ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "M", lambda = 2.33)$varprop)>=0.9),1)
    for(i in 2:(length(trading_time)-1))
    {
      # mse
      rr_point_update_M_mse_forecasts[i] = dynupdate(data = fts_object, newdata = fts_return[1:i,(jk+85)], holdoutdata = fts_return[(i+1):1621,(jk+85)],
                      method = "ridge", pcdmethod = pcdmethod, lambda = rr_find_lambda_M_mse[ceiling((i-1)/10)], order = ftsm_object_order, interval = FALSE, value = TRUE)
      
      
      # mae
      rr_point_update_M_mae_forecasts[i] = dynupdate(data = fts_object, newdata = fts_return[1:i,(jk+85)], holdoutdata = fts_return[(i+1):1621,(jk+85)],
                          method = "ridge",  pcdmethod = pcdmethod, lambda = rr_find_lambda_M_mae[ceiling((i-1)/10)], order = ftsm_object_order, interval = FALSE, value = TRUE)
      
    }
    return(list(rr_point_update_M_mae_forecasts = rr_point_update_M_mae_forecasts, rr_point_update_M_mse_forecasts = rr_point_update_M_mse_forecasts))
  }
}


# updating errors, pcdmethod = "classical"

htm = proc.time()

cl <- makeCluster(30) 
registerDoParallel(cl)

rr_point_update_forecasts_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% RR_update_forecasts(fts_return, "classical", jk = j)

save(rr_point_update_forecasts_list, file = "rr_point_update_forecasts.RData")

htm_RR = proc.time() - htm

# updating errors, pcdmethod = "M"

htm = proc.time()

cl <- makeCluster(30) 
registerDoParallel(cl)

rr_point_update_M_forecasts_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% RR_update_forecasts(fts_return, "M", jk = j)

save(rr_point_update_M_forecasts_list, file = "rr_point_update_M_forecasts.RData")

htm_RR_M = proc.time() - htm

#################
# Chiou's method 
#################

Chiou_pca_update_forecasts <- function(data_return, pcd_method = c("classical", "M"), jk, robust_lambda = 2.33)
{
  pcd_method = match.arg(pcd_method)
  if(pcd_method == "classical")
  {
    # creating storage for errors
    Chiou_pca_mmeu = vector("numeric", length(trading_time))
    Chiou_pca_mmeo = vector("numeric", length(trading_time))
    Chiou_pca_mcpdc = vector("numeric", length(trading_time))
    # updating errors
    for(i in 6:1620)
    {
      dum = dynamic_Chiou_pca(dat = data_return[,1:(jk+84)], newdata = matrix(data_return[1:i, (jk+85)], ncol=i), 
                              holdoutdata = matrix(data_return[(i+1):1621,(jk+85)], ncol=1621-i),
                              order_k_percent = 0.9, order_m_percent = 0.9, pointfore = TRUE)
      Chiou_pca_mmeu[i] = mmeu(dum$holdoutdata, dum$update_forecast)
      Chiou_pca_mmeo[i] = mmeo(dum$holdoutdata, dum$update_forecast)
      Chiou_pca_mcpdc[i] = mcpdc(dum$holdoutdata, dum$update_forecast)
    }
    return(list(Chiou_pca_mmeu = Chiou_pca_mmeu, Chiou_pca_mmeo = Chiou_pca_mmeo, Chiou_pca_mcpdc= Chiou_pca_mcpdc))
  }
  if(pcd_method == "M")
  {
    # creating storage for errors
    Chiou_pca_M_mmeu = vector("numeric", length(trading_time))
    Chiou_pca_M_mmeo = vector("numeric", length(trading_time))
    Chiou_pca_M_mcpdc = vector("numeric", length(trading_time))
    # updating errors
    for(i in 6:1620)
    {
      dum = dynamic_Chiou_pca(dat = data_return[,1:(jk+84)], newdata = matrix(data_return[1:i, (jk+85)], ncol=i), 
                              holdoutdata = matrix(data_return[(i+1):1621,(jk+85)], ncol=1621-i),
                              order_k_percent = 0.9, order_m_percent = 0.9, pcd_method = "M", robust_lambda = robust_lambda, pointfore = TRUE)
      Chiou_pca_M_mmeu[i] = mmeu(dum$holdoutdata, dum$update_forecast)
      Chiou_pca_M_mmeo[i] = mmeo(dum$holdoutdata, dum$update_forecast)
      Chiou_pca_M_mcpdc[i] = mcpdc(dum$holdoutdata, dum$update_forecast)
    }
    return(list(Chiou_pca_M_mmeu = Chiou_pca_M_mmeu, Chiou_pca_M_mmeo = Chiou_pca_M_mmeo, Chiou_pca_M_mcpdc= Chiou_pca_M_mcpdc))
  }
}


# pcdmethod = "classical"
htm = proc.time()
cl <- makeCluster(30) 
registerDoParallel(cl)

Chiou_point_update_forecasts_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% Chiou_pca_update_forecasts(fts_return, pcd_method = "classical", jk = j)

# record results
Chiou_pca_mmeu = Chiou_pca_mmeo = Chiou_pca_mcpdc = matrix(, 1621, 40)
for (j in 1:40)
{
  Chiou_pca_mmeu[,j] = Chiou_point_update_forecasts_list[[j]]$Chiou_pca_mmeu
  Chiou_pca_mmeo[,j] = Chiou_point_update_forecasts_list[[j]]$Chiou_pca_mmeo
  Chiou_pca_mcpdc[,j] = Chiou_point_update_forecasts_list[[j]]$Chiou_pca_mcpdc
}

Chiou_point_update_forecasts = cbind(round(rowMeans(Chiou_pca_mmeu), 4), round(rowMeans(Chiou_pca_mmeo), 4), round(rowMeans(Chiou_pca_mcpdc), 4))
colnames(Chiou_point_update_forecasts) = c("mmeu", "mmeo", "mcpdc")

save(Chiou_point_update_forecasts, file = "Chiou_point_update_forecasts.RData")

htm_Chiou = proc.time() - htm

# pcdmethod = "M"
htm = proc.time()

cl <- makeCluster(30) 
registerDoParallel(cl)

Chiou_point_update_M_forecasts_list = foreach(j = 1:40, .packages = c("ftsa")) %dopar% Chiou_pca_update_forecasts(fts_return, "M", jk = j)

# record results
Chiou_pca_M_mmeu = Chiou_pca_M_mmeo = Chiou_pca_M_mcpdc = matrix(, 1621, 40)
for (j in 1:40)
{
  Chiou_pca_M_mmeu[,j] = Chiou_point_update_M_forecasts_list[[j]]$Chiou_pca_M_mmeu
  Chiou_pca_M_mmeo[,j] = Chiou_point_update_M_forecasts_list[[j]]$Chiou_pca_M_mmeo
  Chiou_pca_M_mcpdc[,j] = Chiou_point_update_M_forecasts_list[[j]]$Chiou_pca_M_mcpdc
}

Chiou_point_update_M_forecasts = cbind(round(rowMeans(Chiou_pca_M_mmeu), 4), round(rowMeans(Chiou_pca_M_mmeo), 4), round(rowMeans(Chiou_pca_M_mcpdc), 4))
colnames(Chiou_point_update_M_forecasts) = c("mmeu", "mmeo", "mcpdc")

save(Chiou_point_update_M_forecasts, file = "Chiou_point_update_M_forecasts.RData")
htm_Chiou_M = proc.time() - htm




#############################
# original plotting commands
#############################

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







