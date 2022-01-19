#################################
# summarising interval forecasts
#################################
library(meboot)

dim(fts_return)

fts_object = fts(trading_time, fts_return[,1:124])
ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "classical")$varprop)>=0.9),1)
ts_PI = forecast(ftsm(fts_object, order = ftsm_object_order), h = 1, method = "arima", pimethod = "nonparametric", B = 1000)

ik = 840
bm_PI = dynupdate(data = fts_object, newdata = fts_return[1:ik,125], order = ftsm_object_order,
                  holdoutdata = fts_return[(ik+1):1621,125], method = "block", value = TRUE, 
                  interval = TRUE, pimethod = "nonparametric")

pls_PI = dynupdate(data = fts_object, newdata = fts_return[1:ik,125], order = ftsm_object_order,
                   holdoutdata = fts_return[(ik+1):1621,125], method = "pls", value = TRUE, 
                   lambda = pls_find_lambda_mse[ceiling((ik-1)/10)],
                   interval = TRUE, pimethod = "nonparametric")

ts_PI_lb_ub = apply(ts_PI$bootsamp[,,1], 1, quantile, c(0.1,0.9))
bm_PI_lb_ub = apply(bm_PI$boot_samp, 1, quantile, c(0.1,0.9))
pls_PI_lb_ub = t(cbind(pls_PI$low$y, pls_PI$up$y))

Chiou_interval = dynamic_Chiou_pca(dat = fts_object$y, newdata = matrix(fts_return[1:ik, 125], ncol = ik),
                               holdoutdata = matrix(fts_return[(ik+1):1621,125], ncol = 1621-ik), 
                               pcd_method = "classical", pointfore = FALSE, bootrep = 1000)



savepdf("PI_update_111",width=12,height=10,toplines=0.9)
plot(trading_time[841:1621], ts_PI_lb_ub[1,841:1621], type = "l", ylim = range(ts_PI_lb_ub), col = 2, lty = 1, xlab="Trading time", ylab="Cumulative intraday return")
lines(trading_time[841:1621], ts_PI_lb_ub[2,841:1621], col = 2, lty = 1)
lines(trading_time[841:1621], bm_PI_lb_ub[1,1:781], col = 4, lty = 2)
lines(trading_time[841:1621], bm_PI_lb_ub[2,1:781], col = 4, lty = 2)
lines(trading_time[841:1621], pls_PI_lb_ub[1,1:781], col = "purple", lty = 4)
lines(trading_time[841:1621], pls_PI_lb_ub[2,1:781], col = "purple", lty = 4)
lines(trading_time[841:1621], Chiou_interval[1,1:781], col = "cyan", lty = 5)
lines(trading_time[841:1621], Chiou_interval[2,1:781], col = "cyan", lty = 5)
points(trading_time[841:1621], fts_return[841:1621,125], col = "black", pch = 1, cex = 0.5)
legend(13,-2, c("Observations", "TS prediction interval", "BM prediction interval", 
                    "PLS prediction interval", "FLR prediction interval"),
       col = c(1, 2, 4, "purple","cyan"), lty = c(NA, 1, 2, 4, 5), cex = 0.5, pch = c(1, NA, NA, NA, NA))
dev.off()


#