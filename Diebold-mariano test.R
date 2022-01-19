###################################
# Conduct the Diebold-mariano test
###################################
library(forecast)
library(stargazer)

# Compare point forecasts using the squared loss function

## FLR vs TS
DM_sqloss_1 = dm.test(e1 = point_update_mse_all[,5], e2 = point_update_mse_all[,1], h = 1, alternative = "less")

## FLR vs BM
DM_sqloss_2 = dm.test(e1 = point_update_mse_all[,5], e2 = point_update_mse_all[,2], h = 1, alternative = "less")

## FLR vs PLS
DM_sqloss_3 = dm.test(e1 = point_update_mse_all[,5], e2 = point_update_mse_all[,3], h = 1, alternative = "less")

## FLR vs RR
DM_sqloss_4 = dm.test(e1 = point_update_mse_all[,5], e2 = point_update_mse_all[,4], h = 1, alternative = "less")


# Compare point forecasts using the absolute loss function

## FLR vs TS
DM_absloss_1 = dm.test(e1 = point_update_mae_all[,5], e2 = point_update_mae_all[,1], h = 1, alternative = "less", power = 1)

## FLR vs BM
DM_absloss_2 = dm.test(e1 = point_update_mae_all[,5], e2 = point_update_mae_all[,2], h = 1, alternative = "less", power = 1)

## FLR vs PLS
DM_absloss_3 = dm.test(e1 = point_update_mae_all[,5], e2 = point_update_mae_all[,3], h = 1, alternative = "less", power = 1)

## FLR vs RR
DM_absloss_4 = dm.test(e1 = point_update_mae_all[,5], e2 = point_update_mae_all[,4], h = 1, alternative = "less", power = 1)


# make a table summarising the results

e1_all = rep("FLR", 8)
e2_all = c("TS","BM","PLS","RR", "TS","BM","PLS","RR")
alternative_all = c(DM_sqloss_1$alternative, DM_sqloss_2$alternative, DM_sqloss_3$alternative, DM_sqloss_4$alternative,
                    DM_absloss_1$alternative, DM_absloss_2$alternative, DM_absloss_3$alternative, DM_absloss_4$alternative)
statistic_all = c(DM_sqloss_1$statistic, DM_sqloss_2$statistic, DM_sqloss_3$statistic, DM_sqloss_4$statistic,
                  DM_absloss_1$statistic, DM_absloss_2$statistic, DM_absloss_3$statistic, DM_absloss_4$statistic)
pvalue_all = round(c(DM_sqloss_1$p.value, DM_sqloss_2$p.value, DM_sqloss_3$p.value, DM_sqloss_4$p.value,
                  DM_absloss_1$p.value, DM_absloss_2$p.value, DM_absloss_3$p.value, DM_absloss_4$p.value), 4)

DM_table = data.frame("e1" = e1_all[1:4], "e2" = e2_all[1:4], "alternative" = alternative_all[1:4], "statistic" = statistic_all[1:4], "p.value" = pvalue_all[1:4],
                      "statistic" = statistic_all[5:8], "p.value" = pvalue_all[5:8])

stargazer(DM_table, digits = 4, summary = F)

