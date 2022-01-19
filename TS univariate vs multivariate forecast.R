farforecast <- function (object, h = 10, var_type = "const", level = 80, PI = FALSE) 
{
  x = as.numeric(rownames(object$y$y))
  order = ncol(object$basis) - 1
  ar_fit = try(ar(object$coeff[, 2:(order + 1)]), silent = TRUE)
  if (class(ar_fit) == "try-error" | ar_fit$order == 0) {
    var_order = 2
  }
  else {
    var_order = min(5, ar_fit$order)
  }
  if (requireNamespace("vars", quietly = TRUE)) {
    var_pred = predict(vars::VAR(object$coeff[, 2:(order + 1)], p = var_order, type = var_type), n.ahead = h, 
                       ci = level/100)
  }
  else {
    stop("Please install vars")
  }
  qconf <- qnorm(0.5 + level/200)
  meanfcast <- varfcast <- matrix(NA, nrow = h, ncol = order)
  for (i in 1:order) {
    var_fit_pred = var_pred$fcst[[i]]
    meanfcast[, i] = var_fit_pred[, 1]
    varfcast[, i] = ((var_fit_pred[, 3] - var_fit_pred[, 2])/(2 * qconf))^2
  }
  point_fore = object$basis[, 2:(order + 1)] %*% t(meanfcast) + object$basis[, 1]
  colnames(point_fore) = 1:h
  point_fore_fts = fts(x, point_fore, yname = "Forecasts")
  if (PI == TRUE) 
  {
    n.curve = ncol(object$y$y)
    L = max(round(n.curve/5), order)
    insample_fore = matrix(, nrow(object$y$y), (ncol(object$y$y) - L))
    for (i in 1:(ncol(object$y$y) - L)) {
      dum = ftsm(fts(object$y$x, object$y$y[, 1:(L + i - 1)]), order = order)
      var_pred = predict(vars::VAR(dum$coeff[, 2:(order + 1)], p = var_order, type = var_type), n.ahead = 1)
      meanfcast = matrix(NA, nrow = 1, ncol = order)
      for (j in 1:order) {
        var_fit_pred = var_pred$fcst[[j]]
        meanfcast[, j] = var_fit_pred[, 1]
      }
      insample_fore[, i] = dum$basis[, 2:(order + 1)] %*% t(meanfcast) + dum$basis[, 1]
    }
    insample_test = object$y$y[, (L + 1):ncol(object$y$y)]
    resi = insample_test - insample_fore
    lb_resi = apply(resi, 1, quantile, (100 - level)/200, na.rm = TRUE)
    ub_resi = apply(resi, 1, quantile, 1 - (100 - level)/200, na.rm = TRUE)
    lb = point_fore + lb_resi
    ub = point_fore + ub_resi
    colnames(lb) = colnames(ub) = 1:h
    PI_lb = fts(x, lb, yname = "Lower bound")
    PI_ub = fts(x, ub, yname = "Upper bound")
    return(list(point_fore = point_fore_fts, PI_lb = PI_lb, PI_ub = PI_ub))
  }
  else {
    return(point_fore_fts)
  }
}


###################
# Aue et al (2015)
###################

## var_type = "const"

ts_point_fore_var = ts_lb_fore_var = ts_ub_fore_var = matrix(,1621,40)
for(j in 1:40)
{
  fts_object = fts(1:1621, fts_return[,1:(84+j)])
  ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "classical")$varprop)>=0.9),1)
  dum = farforecast(ftsm(fts_object, order = ftsm_object_order), h = 1, var_type = "const", PI = TRUE)
  ts_point_fore_var[,j] = dum$point_fore$y
  ts_lb_fore_var[,j] = dum$PI_lb$y
  ts_ub_fore_var[,j] = dum$PI_ub$y
}

# MAFE and MSFE

ts_point_fore_var_mse = round(colMeans((ts_point_fore_var - fts_return[,86:125])^2), 4)
ts_point_fore_var_mae = round(colMeans(abs(ts_point_fore_var - fts_return[,86:125])), 4)

# interval score

ts_lb_ind_var = ifelse(fts_return[,86:125] < ts_lb_fore_var, 1, 0)
ts_ub_ind_var = ifelse(fts_return[,86:125] > ts_ub_fore_var, 1, 0)
ts_score_fore_var = round(colMeans((ts_ub_fore_var - ts_lb_fore_var) + 
                                     2/0.2 * ((ts_lb_fore_var - fts_return[,86:125]) * ts_lb_ind_var + 
                                                (fts_return[,86:125] - ts_ub_fore_var) * ts_ub_ind_var)), 4)

## var_type = "trend"

ts_point_fore_var_trend = ts_lb_fore_var_trend = ts_ub_fore_var_trend = matrix(,1621,40)
for(j in 1:40)
{
  fts_object = fts(1:1621, fts_return[,1:(84+j)])
  ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "classical")$varprop)>=0.9),1)
  dum = farforecast(ftsm(fts_object, order = ftsm_object_order), h = 1, var_type = "trend", PI = TRUE)
  ts_point_fore_var_trend[,j] = dum$point_fore$y
  ts_lb_fore_var_trend[,j] = dum$PI_lb$y
  ts_ub_fore_var_trend[,j] = dum$PI_ub$y
}

# MAFE and MSFE

ts_point_fore_var_trend_mse = round(colMeans((ts_point_fore_var_trend - fts_return[,86:125])^2), 4)
ts_point_fore_var_trend_mae = round(colMeans(abs(ts_point_fore_var_trend - fts_return[,86:125])), 4)


## var_type = "both"

ts_point_fore_var_both = ts_lb_fore_var_both = ts_ub_fore_var_both = matrix(,1621,40)
for(j in 1:40)
{
  fts_object = fts(1:1621, fts_return[,1:(84+j)])
  ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "classical")$varprop)>=0.9),1)
  dum = farforecast(ftsm(fts_object, order = ftsm_object_order), h = 1, var_type = "both", PI = TRUE)
  ts_point_fore_var_both[,j] = dum$point_fore$y
  ts_lb_fore_var_both[,j] = dum$PI_lb$y
  ts_ub_fore_var_both[,j] = dum$PI_ub$y
}

# MAFE and MSFE

ts_point_fore_var_both_mse = round(colMeans((ts_point_fore_var_both - fts_return[,86:125])^2), 4)
ts_point_fore_var_both_mae = round(colMeans(abs(ts_point_fore_var_both - fts_return[,86:125])), 4)


## var_type = "none"

ts_point_fore_var_none = ts_lb_fore_var_none = ts_ub_fore_var_none = matrix(,1621,40)
for(j in 1:40)
{
  fts_object = fts(1:1621, fts_return[,1:(84+j)])
  ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "classical")$varprop)>=0.9),1)
  dum = farforecast(ftsm(fts_object, order = ftsm_object_order), h = 1, var_type = "none", PI = TRUE)
  ts_point_fore_var_none[,j] = dum$point_fore$y
  ts_lb_fore_var_none[,j] = dum$PI_lb$y
  ts_ub_fore_var_none[,j] = dum$PI_ub$y
}

# MAFE and MSFE

ts_point_fore_var_none_mse = round(colMeans((ts_point_fore_var_none - fts_return[,86:125])^2), 4)
ts_point_fore_var_none_mae = round(colMeans(abs(ts_point_fore_var_none - fts_return[,86:125])), 4)

ts_point_fore_var_all_mse = cbind(ts_point_fore_var_mse, ts_point_fore_var_trend_mse, ts_point_fore_var_both_mse, ts_point_fore_var_none_mse)
ts_point_fore_var_all_mae = cbind(ts_point_fore_var_mae, ts_point_fore_var_trend_mae, ts_point_fore_var_both_mae, ts_point_fore_var_none_mae)
colnames(ts_point_fore_var_all_mse) = colnames(ts_point_fore_var_all_mae) = c("const","trend","both","none")


###########################
# Hyndman and Ullah (2007)
###########################

ts_point_fore_uni = ts_lb_fore_uni = ts_ub_fore_uni = matrix(,1621,40)
for(j in 1:40)
{
  fts_object = fts(1:1621, fts_return[,1:(84+j)])
  ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "classical")$varprop)>=0.9),1)
  dum = forecast(ftsm(fts_object, order = ftsm_object_order), h = 1, method = "arima", 
                 pimethod = "nonparametric", B = 1000)
  ts_point_fore_uni[,j] = dum$mean$y
  ts_lb_fore_uni[,j] = dum$lower$y
  ts_ub_fore_uni[,j] = dum$upper$y
}

# MAFE and MSFE

ts_point_fore_uni_mse = round(colMeans((ts_point_fore_uni - fts_return[,86:125])^2), 4)
ts_point_fore_uni_mae = round(colMeans(abs(ts_point_fore_uni - fts_return[,86:125])), 4)

# score

ts_lb_ind = ifelse(fts_return[,86:125] < ts_lb_fore_uni, 1, 0)
ts_ub_ind = ifelse(fts_return[,86:125] > ts_ub_fore_uni, 1, 0)
ts_score_fore_uni = round(colMeans((ts_ub_fore_uni - ts_lb_fore_uni) + 
                                     2/0.2 * ((ts_lb_fore_uni - fts_return[,86:125]) * ts_lb_ind + 
                                                (fts_return[,86:125] - ts_ub_fore_uni) * ts_ub_ind)), 4)


ts_point_fore_uni_var_mse = cbind(ts_point_fore_uni_mse, ts_point_fore_var_mse)
ts_point_fore_uni_var_mae = cbind(ts_point_fore_uni_mae, ts_point_fore_var_mae)
ts_interval_fore_uni_var_score = cbind(ts_score_fore_uni, ts_score_fore_var) 

colnames(ts_point_fore_uni_var_mse) = colnames(ts_point_fore_uni_var_mae) = 
  colnames(ts_interval_fore_uni_var_score) = c("ARIMA", "VAR")

savepdf("Uni_VAR_all", width = 28.5, height = 10, toplines = 0.9, pointsize = 15)
par(mfrow = c(1,3))
boxplot(ts_point_fore_uni_var_mae, ylab = "MAFE")
boxplot(ts_point_fore_uni_var_mse, ylab = "MSFE")
boxplot(ts_interval_fore_uni_var_score, ylab = "Mean interval score")
dev.off()

savepdf("Uni_VAR_all_1", width = 12, height = 10, toplines = 0.9, pointsize = 10)
boxplot(ts_point_fore_uni_var_mae, ylab = "MAFE")
dev.off()

savepdf("Uni_VAR_all_2", width = 12, height = 10, toplines = 0.9, pointsize = 10)
boxplot(ts_point_fore_uni_var_mse, ylab = "MSFE")
dev.off()

savepdf("Uni_VAR_all_3", width = 12, height = 10, toplines = 0.9, pointsize = 10)
boxplot(ts_interval_fore_uni_var_score, ylab = "Mean interval score")
dev.off()


round(colMeans(ts_point_fore_uni_var_mae), 4)
round(colMeans(ts_point_fore_uni_var_mse), 4)
round(colMeans(ts_interval_fore_uni_var_score), 4)

min_ts_arima = cbind(min(ts_point_fore_uni_var_mse[,1]), min(ts_point_fore_uni_var_mae[,1]), min(ts_interval_fore_uni_var_score[,1]))
max_ts_arima = cbind(max(ts_point_fore_uni_var_mse[,1]), max(ts_point_fore_uni_var_mae[,1]), max(ts_interval_fore_uni_var_score[,1]))
median_ts_arima = cbind(median(ts_point_fore_uni_var_mse[,1]), median(ts_point_fore_uni_var_mae[,1]), median(ts_interval_fore_uni_var_score[,1]))
mean_ts_arima = cbind(mean(ts_point_fore_uni_var_mse[,1]), mean(ts_point_fore_uni_var_mae[,1]), mean(ts_interval_fore_uni_var_score[,1]))

table_ts_arima = rbind(min_ts_arima, median_ts_arima, mean_ts_arima, max_ts_arima)
rownames(table_ts_arima) = c("Min", "Median", "Mean", "Max")
colnames(table_ts_arima) = c("MSFE", "MAFE", "Score")


min_ts_var = cbind(min(ts_point_fore_uni_var_mse[,2]), min(ts_point_fore_uni_var_mae[,2]), min(ts_interval_fore_uni_var_score[,2]))
max_ts_var = cbind(max(ts_point_fore_uni_var_mse[,2]), max(ts_point_fore_uni_var_mae[,2]), max(ts_interval_fore_uni_var_score[,2]))
median_ts_var = cbind(median(ts_point_fore_uni_var_mse[,2]), median(ts_point_fore_uni_var_mae[,2]), median(ts_interval_fore_uni_var_score[,2]))
mean_ts_var = cbind(mean(ts_point_fore_uni_var_mse[,2]), mean(ts_point_fore_uni_var_mae[,2]), mean(ts_interval_fore_uni_var_score[,2]))

table_ts_var = rbind(min_ts_var, median_ts_var, mean_ts_var, max_ts_var)
rownames(table_ts_var) = c("Min", "Median", "Mean", "Max")
colnames(table_ts_var) = c("MSFE", "MAFE", "Score")

table_ts_all = cbind(table_ts_arima[,1], table_ts_var[,1], table_ts_arima[,2], table_ts_var[,2], table_ts_arima[,3], table_ts_var[,3])
stargazer(table_ts_all, digits = 4, summary = F)

##################################
# pcdmethod = c("classical", "M")
##################################

# lambda = 2.33 and 1.81 have been examined

ts_point_fore_uni_M = ts_lb_fore_uni_M = ts_ub_fore_uni_M = matrix(,1621,40)
for(j in 1:40)
{
  fts_object = fts(1:1621, fts_return[,1:(84+j)])
  ftsm_object_order = head(which(cumsum(ftsm(fts_object, method = "M", lambda = 2.33)$varprop)>=0.9),1)
  dum = forecast(ftsm(fts_object, method = "M", order = ftsm_object_order, lambda = 2.33), h = 1, method = "arima", pimethod = "nonparametric", B = 1000)
  ts_point_fore_uni_M[,j] = dum$mean$y
  ts_lb_fore_uni_M[,j] = dum$lower$y
  ts_ub_fore_uni_M[,j] = dum$upper$y
}

# MAFE and MSFE

ts_point_fore_uni_M_mse = round(colMeans((ts_point_fore_uni_M - fts_return[,86:125])^2), 4)
ts_point_fore_uni_M_mae = round(colMeans(abs(ts_point_fore_uni_M - fts_return[,86:125])), 4)

# score

ts_lb_ind_M = ifelse(fts_return[,86:125] < ts_lb_fore_uni_M, 1, 0)
ts_ub_ind_M = ifelse(fts_return[,86:125] > ts_ub_fore_uni_M, 1, 0)
ts_score_fore_uni_M = round(colMeans((ts_ub_fore_uni_M - ts_lb_fore_uni_M) + 
                                       2/0.2 * ((ts_lb_fore_uni_M - fts_return[,86:125]) * ts_lb_ind_M + 
                                                  (fts_return[,86:125] - ts_ub_fore_uni_M) * ts_ub_ind_M)), 4)


##########
# RobRSVD 
##########

ts_point_fore_uni_RobRSVD = ts_lb_fore_uni_RobRSVD = ts_ub_fore_uni_RobRSVD = matrix(,1621,40)
for(j in 1:40)
{
  fts_object = fts(1:1621, fts_return[,1:(84+j)])
  dum = forecast(ftsm_robrsvd(fts_object), h = 1, method = "arima", B = 1000)
  ts_point_fore_uni_RobRSVD[,j] = dum$mean$y
  ts_lb_fore_uni_RobRSVD[,j] = dum$lower$y
  ts_ub_fore_uni_RobRSVD[,j] = dum$upper$y
}

# MAFE and MSFE

ts_point_fore_uni_mse_RobRSVD = round(colMeans((ts_point_fore_uni_RobRSVD  - fts_return[,86:125])^2), 4)
ts_point_fore_uni_mae_RobRSVD = round(colMeans(abs(ts_point_fore_uni_RobRSVD  - fts_return[,86:125])), 4)

# score

ts_lb_ind_RobRSVD = ifelse(fts_return[,86:125] < ts_lb_fore_uni_RobRSVD, 1, 0)
ts_ub_ind_RobRSVD = ifelse(fts_return[,86:125] > ts_ub_fore_uni_RobRSVD, 1, 0)
ts_score_fore_uni_RobRSVD = round(colMeans((ts_ub_fore_uni_RobRSVD - ts_lb_fore_uni_RobRSVD) + 
                                     2/0.2 * ((ts_lb_fore_uni_RobRSVD - fts_return[,86:125]) * ts_lb_ind_RobRSVD + 
                                                (fts_return[,86:125] - ts_ub_fore_uni_RobRSVD) * ts_ub_ind_RobRSVD)), 4)

# summary

ts_point_fore_uni_all_mse = cbind(ts_point_fore_uni_mse, ts_point_fore_uni_M_mse, ts_point_fore_uni_mse_RobRSVD)
ts_point_fore_uni_all_mae = cbind(ts_point_fore_uni_mae, ts_point_fore_uni_M_mae, ts_point_fore_uni_mae_RobRSVD)
ts_interval_fore_uni_all_score = cbind(ts_score_fore_uni, ts_score_fore_uni_M, ts_score_fore_uni_RobRSVD)
colnames(ts_interval_fore_uni_all_score) = colnames(ts_point_fore_uni_all_mse) = colnames(ts_point_fore_uni_all_mae) = c("FPCA","M-FPCA")

savepdf("TS_point_interval", width = 28.5, height = 10, toplines = 0.9, pointsize = 15)
par(mfrow = c(1,3))
boxplot(ts_point_fore_uni_all_mae, ylab = "MAFE", notch = TRUE)
boxplot(ts_point_fore_uni_all_mse, ylab = "MSFE", notch = TRUE)
boxplot(ts_interval_fore_uni_all_score, ylab = "Mean interval score", notch = TRUE)
dev.off()

savepdf("TS_point_interval_1", width = 12, height = 10, toplines = 0.9, pointsize = 10)
boxplot(ts_point_fore_uni_all_mae, ylab = "MAFE", notch = TRUE)
dev.off()

savepdf("TS_point_interval_2", width = 12, height = 10, toplines = 0.9, pointsize = 10)
boxplot(ts_point_fore_uni_all_mse, ylab = "MSFE", notch = TRUE)
dev.off()

savepdf("TS_point_interval_3", width = 12, height = 10, toplines = 0.9, pointsize = 10)
boxplot(ts_interval_fore_uni_all_score, ylab = "Mean interval score", notch = TRUE)
dev.off()

round(colMeans(ts_point_fore_uni_all_mae),4)
round(colMeans(ts_point_fore_uni_all_mse),4)
round(colMeans(ts_interval_fore_uni_all_score),4)


min_method_fpca = cbind(min(ts_point_fore_uni_all_mse[,1]), min(ts_point_fore_uni_all_mae[,1]), min(ts_interval_fore_uni_all_score[,1]))
median_method_fpca = cbind(median(ts_point_fore_uni_all_mse[,1]), median(ts_point_fore_uni_all_mae[,1]), median(ts_interval_fore_uni_all_score[,1]))
mean_method_fpca = cbind(mean(ts_point_fore_uni_all_mse[,1]), mean(ts_point_fore_uni_all_mae[,1]), mean(ts_interval_fore_uni_all_score[,1]))
max_method_fpca = cbind(max(ts_point_fore_uni_all_mse[,1]), max(ts_point_fore_uni_all_mae[,1]), max(ts_interval_fore_uni_all_score[,1]))

min_method_mfpca = cbind(min(ts_point_fore_uni_all_mse[,2]), min(ts_point_fore_uni_all_mae[,2]), min(ts_interval_fore_uni_all_score[,2]))
median_method_mfpca = cbind(median(ts_point_fore_uni_all_mse[,2]), median(ts_point_fore_uni_all_mae[,2]), median(ts_interval_fore_uni_all_score[,2]))
mean_method_mfpca = cbind(mean(ts_point_fore_uni_all_mse[,2]), mean(ts_point_fore_uni_all_mae[,2]), mean(ts_interval_fore_uni_all_score[,2]))
max_method_mfpca = cbind(max(ts_point_fore_uni_all_mse[,2]), max(ts_point_fore_uni_all_mae[,2]), max(ts_interval_fore_uni_all_score[,2]))

min_method_robrsvd = cbind(min(ts_point_fore_uni_all_mse[,3]), min(ts_point_fore_uni_all_mae[,3]), min(ts_interval_fore_uni_all_score[,3]))
median_method_robrsvd = cbind(median(ts_point_fore_uni_all_mse[,3]), median(ts_point_fore_uni_all_mae[,3]), median(ts_interval_fore_uni_all_score[,3]))
mean_method_robrsvd = cbind(mean(ts_point_fore_uni_all_mse[,3]), mean(ts_point_fore_uni_all_mae[,3]), mean(ts_interval_fore_uni_all_score[,3]))
max_method_robrsvd = cbind(max(ts_point_fore_uni_all_mse[,3]), max(ts_point_fore_uni_all_mae[,3]), max(ts_interval_fore_uni_all_score[,3]))

table_method_fpca = rbind(min_method_fpca, median_method_fpca, mean_method_fpca, max_method_fpca)
table_method_mfpca = rbind(min_method_mfpca, median_method_mfpca, mean_method_mfpca, max_method_mfpca)
table_method_robrsvd = rbind(min_method_robrsvd, median_method_robrsvd, mean_method_robrsvd, max_method_robrsvd)

table_method_all = cbind(table_method_fpca[,1], table_method_mfpca[,1], table_method_robrsvd[,1],
                         table_method_fpca[,2], table_method_mfpca[,2], table_method_robrsvd[,2],
                         table_method_fpca[,3], table_method_mfpca[,3], table_method_robrsvd[,3])

rownames(table_method_all) = c("Min", "Median", "Mean", "Max")
colnames(table_method_all) = rep(colnames(ts_point_fore_uni_all_mse),3)

stargazer(table_method_all, digits = 4, summary = F)
