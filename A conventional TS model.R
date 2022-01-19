#######################################################################
# Prediction using the conentional ARIMA model suggested by Reviewer 1
#######################################################################
library(ftsa)
library(forecast)

# (1621 series with maximum of 125 observarions each) Making forecasts
ts.fore = matrix(NA, nrow = 1621, ncol = 40)

for (i in 1:40)
{
 for(j in 1:1621)
 {
   fit = auto.arima(fts_return[j,1:(84+i)])
   ts.fore[j,i] =  forecast(fit, h = 1)$mean
 }
}

# Computing errors
mae  = ftsa:::mae
mse  = ftsa:::mse

ts_fore_mae = ts_fore_mse = matrix(NA, 1620, 40)
for(j in 1:40)
{
  for(i in 1:1620)
  {
    ts_fore_mae[i,j] = mae(fts_return[(i+1):1621,(85+j)], ts.fore[(i+1):1621,j])
    ts_fore_mse[i,j] = mse(fts_return[(i+1):1621,(85+j)], ts.fore[(i+1):1621,j])
  }
}

ts_fore_err= cbind(round(rowMeans(ts_fore_mse), 4), round(rowMeans(ts_fore_mae), 4))
colnames(ts_fore_err) = c("mse", "mae")

savepdf("dynamic_update_1_withts",width=12,height=10,toplines=0.9,pointsize=10)
plot(trading_time2, point_update_mse_all[,1], type="l", ylim = c(3, 35) , xlab="Time of Day", ylab="MSFE")
points(trading_time2, point_update_mse_all[,1], col=1, cex=0.25)

lines(trading_time2, point_update_mse_all[,2], col=2)
points(trading_time2, point_update_mse_all[,2], col=2, pch=2, cex=0.25)

lines(trading_time2, point_update_mse_all[,3], col=3)
points(trading_time2, point_update_mse_all[,3], col=3, pch=3, cex=0.25)

lines(trading_time2, point_update_mse_all[,4], col=4)
points(trading_time2, point_update_mse_all[,4], col=4, pch=4, cex=0.25)

lines(trading_time2[3:1617], point_update_mse_all[3:1617,5], col=5)
points(trading_time2[3:1617], point_update_mse_all[3:1617,5], col=5, pch=5, cex=0.25)

lines(trading_time2, ts_fore_err[2:1620,1], col=6)
points(trading_time2, ts_fore_err[2:1620,1], col=6, pch =6, cex=0.25)

legend("bottomleft",c("TS","BM","PLS","RR","FLR", "ARIMA"),cex=0.6,col=1:6, pch = 1:6,lty=rep(1,6))
dev.off()


savepdf("dynamic_update_2_withts",width=12,height=10,toplines=0.9,pointsize=10)
plot(trading_time2, point_update_mae_all[,1], type = "l", ylim = c(1.2 ,4), xlab="Time of Day", ylab="MAFE", cex=0.5)

points(trading_time2, point_update_mae_all[,1], col=1, cex = 0.25)

lines(trading_time2, point_update_mae_all[,2], col=2)
points(trading_time2, point_update_mae_all[,2], col=2, pch=2, cex=0.25)

lines(trading_time2, point_update_mae_all[,3], col=3)
points(trading_time2, point_update_mae_all[,3], col=3, pch=3, cex=0.25)

lines(trading_time2, point_update_mae_all[,4], col=4)
points(trading_time2, point_update_mae_all[,4], col=4, pch=4, cex=0.25)

lines(trading_time2[3:1617], point_update_mae_all[3:1617,5], col=5)
points(trading_time2[3:1617], point_update_mae_all[3:1617,5], col=5, pch=5, cex=0.25)

lines(trading_time2, ts_fore_err[2:1620,2], col=6)
points(trading_time2, ts_fore_err[2:1620,2], col=6, pch =6, cex=0.25)

legend("bottomleft",c("TS","BM","PLS","RR","FLR", "ARIMA"),cex=0.6,col=1:6, pch = 1:6, lty=rep(1,6))
dev.off()