# compute pooled beta estimates and SEs for MULTIPLE IMPUTATION ANALYSIS
# and pooled odds ratios, 95% confidence intervals, and diagnostics (RIV etc.)
# pooled estimates are computed using Rubin's Rules 

library(tidyverse)

jobnum <- "1726"
M <- 20

# -----------------------------------------------------------------------------
# -------------------         MI BY CART        -------------------------------
# -----------------------------------------------------------------------------

modsOR_mi_cart <- readRDS(paste0("~cluster-scratch/SARTCORS/results",jobnum,"/modsOR_mi_df_1.RDS"))

modsOR_mi_cart %>% count(dat, imp_meth)

mi_pooled0 <- modsOR_mi_cart %>%
  group_by(outcome, predictors_all, rowname) %>%
  summarize(beta_sums = sum(beta_hat_gee1)
    , beta_pooled = beta_sums/M
    # average within imputation variance
    , W_hat_sum = sum(se_gee1^2)
    , W_hat = W_hat_sum/M)

mi_pooled <- modsOR_mi_cart %>%
  left_join(mi_pooled0, by=c("outcome", "predictors_all", "rowname")) %>%
  mutate(# between imputation variance
    diff = beta_hat_gee1 - beta_pooled
    , diff_sq = diff^2) %>%
  group_by(outcome, predictors_all, rowname, beta_pooled, W_hat) %>%
  summarize(B_hat = sum(diff_sq)/(M-1)) %>%
  mutate(# variance of pooled estimates for beta hats
    V_hat = W_hat + ((M+1)/M)*B_hat
    , se_pooled = sqrt(V_hat)) 

# compute additional measures
# https://bookdown.org/mwheymans/bookmi/rubins-rules.html
# https://bookdown.org/mwheymans/bookmi/measures-of-missing-data-information.html#eq:lambda

# full sample size
n <- 806045

mi_pooled_more <- mi_pooled %>%
  mutate(# fraction of missing information (lambda)
    # interpreted as the proportion of variation in the parameter of interest due to the missing data
    lambda = (B_hat + (B_hat/M))/V_hat
    # relative increase in variance (due to non-response)
    # interpreted as the proportional increase in the sampling variance of the parameter
    # of interest due to the missing data
    , RIV = B_hat + (B_hat/M)/W_hat
    # fraction of missing information (difference between FMI and lambda is that FMI
    # is adjusted for the fact that the number of imputed datasets generated is not
    # unlimitedly large; differs from lambda when df is small)
    # use df_adjusted (which is what is used in mice R package)
    , df_old = (M-1)/(lambda^2)
    , df_observed = (((n - M) + 1)/(n - M + 3))*(n-M)*(1-lambda)
    , df_adjusted = (df_old*df_observed)/(df_old + df_observed) 
    , FMI = (RIV + 2/(df_adjusted+3))/(1 + RIV)
    # relative efficiency
    , RE = 1/(1 + (FMI/M)) 
    # OR + CIs
    # need to use t-dist w/ adjusted df in CI
    , OR_pooled = exp(beta_pooled)
    , t_crit = qt(p=0.975, df=df_adjusted)
    , lower_pooled = exp(beta_pooled - t_crit*se_pooled)
    , upper_pooled = exp(beta_pooled + t_crit*se_pooled))


# dist'n of relative increase in variance (due to non-response)
# (proportional increase in the sampling variance of the parameter
# of interest due to the missing data)
mosaic::favstats(~RIV, data=mi_pooled_more)
mosaic::favstats(~RIV, data=filter(mi_pooled_more,str_detect(rowname,"patient_race_ethnic")))

# proportion of variation in the parameter of interest due to the missing data
mosaic::favstats(~lambda, data=mi_pooled_more)
mosaic::favstats(~lambda, data=filter(mi_pooled_more,str_detect(rowname,"patient_race_ethnic")))

# fraction of missing information
mosaic::favstats(~FMI, data=mi_pooled_more)
mosaic::favstats(~FMI, data=filter(mi_pooled_more,str_detect(rowname,"patient_race_ethnic")))

# relative efficiency
mosaic::favstats(~RE, data=mi_pooled_more)
mosaic::favstats(~RE, data=filter(mi_pooled_more,str_detect(rowname,"patient_race_ethnic")))

# ----------------------- save "final" MI pooled estimates + results ----------

saveRDS(mi_pooled_more, paste0("~cluster-scratch/SARTCORS/results",jobnum,"/modsOR_mi_pooled_cart",jobnum,".RDS"))


# -----------------------------------------------------------------------------
# -------------------         MI BY DEFAULT     -------------------------------
# -----------------------------------------------------------------------------

rm(mi_pooled0, mi_pooled, mi_pooled_more)

modsOR_mi_def <- readRDS(paste0("~cluster-scratch/SARTCORS/results",jobnum,"/modsOR_mi_df_0.RDS"))
  
modsOR_mi_def %>% count(dat, imp_meth)

mi_pooled0 <- modsOR_mi_def %>%
  group_by(outcome, predictors_all, rowname) %>%
  summarize(beta_sums = sum(beta_hat_gee1)
         , beta_pooled = beta_sums/M
         # average within imputation variance
         , W_hat_sum = sum(se_gee1^2)
         , W_hat = W_hat_sum/M)

mi_pooled <- modsOR_mi_def %>%
  left_join(mi_pooled0, by=c("outcome", "predictors_all", "rowname")) %>%
  mutate(# between imputation variance
         diff = beta_hat_gee1 - beta_pooled
         , diff_sq = diff^2) %>%
  group_by(outcome, predictors_all, rowname, beta_pooled, W_hat) %>%
  summarize(B_hat = sum(diff_sq)/(M-1)) %>%
  mutate(# variance of pooled estimates for beta hats
         V_hat = W_hat + ((M+1)/M)*B_hat
         , se_pooled = sqrt(V_hat)) 

# compute additional measures
# https://bookdown.org/mwheymans/bookmi/rubins-rules.html
# https://bookdown.org/mwheymans/bookmi/measures-of-missing-data-information.html#eq:lambda

# full sample size
n <- 806045
  
mi_pooled_more <- mi_pooled %>%
  mutate(# fraction of missing information (lambda)
    # interpreted as the proportion of variation in the parameter of interest due to the missing data
    lambda = (B_hat + (B_hat/M))/V_hat
    # relative increase in variance (due to non-response)
    # interpreted as the proportional increase in the sampling variance of the parameter
    # of interest due to the missing data
    , RIV = B_hat + (B_hat/M)/W_hat
    # fraction of missing information (difference between FMI and lambda is that FMI
    # is adjusted for the fact that the number of imputed datasets generated is not
    # unlimitedly large; differs from lambda when df is small)
    # use df_adjusted (which is what is used in mice R package)
    , df_old = (M-1)/(lambda^2)
    , df_observed = (((n - M) + 1)/(n - M + 3))*(n-M)*(1-lambda)
    , df_adjusted = (df_old*df_observed)/(df_old + df_observed) 
    , FMI = (RIV + 2/(df_adjusted+3))/(1 + RIV)
    # relative efficiency
    , RE = 1/(1 + (FMI/M)) 
    # OR + CIs
    # need to use t-dist w/ adjusted df in CI
    , OR_pooled = exp(beta_pooled)
    , t_crit = qt(p=0.975, df=df_adjusted)
    , lower_pooled = exp(beta_pooled - t_crit*se_pooled)
    , upper_pooled = exp(beta_pooled + t_crit*se_pooled))


# dist'n of relative increase in variance (due to non-response)
# (proportional increase in the sampling variance of the parameter
# of interest due to the missing data)
mosaic::favstats(~RIV, data=mi_pooled_more)
mosaic::favstats(~RIV, data=filter(mi_pooled_more,str_detect(rowname,"patient_race_ethnic")))

# proportion of variation in the parameter of interest due to the missing data
mosaic::favstats(~lambda, data=mi_pooled_more)
mosaic::favstats(~lambda, data=filter(mi_pooled_more,str_detect(rowname,"patient_race_ethnic")))

# fraction of missing information
mosaic::favstats(~FMI, data=mi_pooled_more)
mosaic::favstats(~FMI, data=filter(mi_pooled_more,str_detect(rowname,"patient_race_ethnic")))

# relative efficiency
mosaic::favstats(~RE, data=mi_pooled_more)
mosaic::favstats(~RE, data=filter(mi_pooled_more,str_detect(rowname,"patient_race_ethnic")))

# ----------------------- save "final" MI pooled estimates + results ----------

saveRDS(mi_pooled_more, paste0("~cluster-scratch/SARTCORS/results",jobnum,"/modsOR_mi_pooled_def",jobnum,".RDS"))


