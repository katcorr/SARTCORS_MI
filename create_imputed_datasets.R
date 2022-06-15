# TO BE RUN ON CLUSTER x 20
# created imputed datasets using multiple imputation

args = commandArgs(trailingOnly=TRUE)
results = args[1]
Cluster = args[2]
Process = args[3]

Process

cluster = strtoi(Cluster)
process = strtoi(Process)
cantor = (cluster + process)*(cluster + process + 1)/2 + process

# this will be the random seed to begin each mice process
print(cantor)

library(tidyverse)
library(mice)

# FULL SARTCORS dataset
#sart_full <- readRDS("Z:/SARTCORS/sart_full.RDS")
sart_full <- readRDS("sart_full.RDS")

# ------------------------------------------------------------------------------
# ---------------------- MULTIPLE IMPUTATION -----------------------------------
# ------------------------------------------------------------------------------

# -----------------------     set up data         ------------------------------

sart_forimp <- sart_full %>%
  # need to make character string vars numeric or factor to include in imputation model
  # also need to make "missing" categories actually NA so can be imputed
  # include all factors will include in analysis model + auxiliary variables
  mutate(race = factor(ifelse(patient_race_ethnic=="Missing"
                              , yes = NA_character_, no=as.character(patient_race_ethnic))
                       , levels=c("Non-Hispanic White"
                                  , "Non-Hispanic Black/African-American"
                                  , "Non-Hispanic Asian"
                                  , "Other"
                                  , "Multiple races"
                                  , "Hispanic/Latino"))
         , patient_age_at_start_sq = patient_age_at_start*patient_age_at_start
         , bmicat4f = factor(ifelse(bmicat4=="Missing"
                                    , yes = NA_character_, no=as.character(bmicat4))
                             , levels=c("Healthy weight", "Underweight", "Overweight"
                                        , "Obese"))
         , fsh_gt10f = factor(ifelse(fsh_gt10=="Missing"
                                     , yes = NA_character_, no=as.character(fsh_gt10))
                              , levels=c("<= 10", "> 10"))
         # amh is not collected for embryo banking or frozen cycles
         # only non-missing if cancelled or fresh cycle
         , amhcat3f = factor(case_when(amhcat3=="Missing" & cycle_type %in% c("Cancelled", "Fresh")
                                       ~ NA_character_
                                       , amhcat3=="Missing" & cycle_type %in% c("Embryo Banking", "Frozen")
                                       ~ "N/A"
                                       , TRUE ~ as.character(amhcat3))
                             , levels=c("1 - <4", "< 1", ">= 4", "N/A"))
         , gravidity_cat4f = as.factor(ifelse(gravidity_cat4=="Unknown"
                                              , yes = NA_character_, no=gravidity_cat4))
         , parousf = as.factor(parous)
         , any_prior_sabf = as.factor(any_prior_sab)
         , frozenf = as.factor(frozen)
         , sab_preg = case_when(preg == 0 ~ "N/A"
                                , TRUE ~ pregnancy_loss_abortion)
         , pregf = as.factor(case_when(is.na(preg) ~ "N/A"
                             , TRUE ~ as.character(preg)))
         , lbf = as.factor(case_when(is.na(lb) ~ "N/A"
                           , TRUE ~ as.character(lb)))
         , across(c(endometriosis, dx_tubal, male_infertility, dx_ovulation
                    , diminished_ovarian_reserve
                    , uterine, unexplained, any_ah, any_icsi, pgd), ~as.factor(.))) %>%
  select(state, mandate_class1, race
         , reporting_year, patient_age_at_start, patient_age_at_start_sq
         , bmicat4f, fsh_gt10f, amhcat3f, endometriosis, dx_tubal
         , male_infertility, dx_ovulation, diminished_ovarian_reserve
         , uterine, unexplained, num_transferred, any_ah, any_icsi
         , parousf, any_prior_sabf, frozenf
         , pregf, lbf
         # auxiliary vars
         , prior_fresh_cycles
         , prior_frozen_cycles, num_retrieved, pgd 
         # for checking derivations
         #, patient_race_ethnic, bmicat4, amhcat3, fsh_gt10, parous
         #, any_prior_sab, gravidity, preg, lb, cycle_type, frozen
         # vars do not want to use in imputation but need for later
         , external_patient_id, external_cycle_id, cancelled, sab_preg
  )


# split the data by mandate status to allow for interaction between race & mandate: 
# https://stefvanbuuren.name/fimd/sec-knowledge.html
# Section 6.4.2
# "Interactions involving categorical variables can be done in similar ways 
# (Van Buuren and Groothuis-Oudshoorn 2011), for example by imputing the data 
# in separate groups. One may do this in mice by splitting the dataset into two
# or more parts, run mice() on each part and then combine the imputed datasets 
# with rbind()."

# getting loggedEvents for all the non-mandated states when look at mandated data
# (even though they're not in the dataset, because it's a factor, those levels exist)
# SO create as factor AFTER splitting up data
m0_forimp <- sart_forimp %>%
  filter(mandate_class1 == 0) %>%
  mutate(statef = as.factor(state)) %>%
  select(-state)

m1_forimp <- sart_forimp %>%
  filter(mandate_class1 == 1) %>%
  mutate(statef = as.factor(state)) %>%
  select(-state)

# checks
# sart_forimp %>% count(parousf, parous)
# sart_forimp %>% count(any_prior_sabf, any_prior_sab)
# sart_forimp %>% count(frozen, cycle_type)
# sart_forimp %>% count(frozenf, frozen, cycle_type)
# sart_forimp %>% count(patient_race_ethnic, race)
# sart_forimp %>% count(bmicat4, bmicat4f)
# sart_forimp %>% count(amhcat3, amhcat3f)
# sart_forimp %>% count(fsh_gt10f, fsh_gt10)
# sart_forimp %>% count(mandate_class1, race, mandate_race)
# m1_forimp %>% count(statef)
# AH is rarely missing; the majority of "missing" are cancelled cycles where it's N/A
# mosaic::tally(any_ah ~ cycle_typef, data=sart_forimp)
# ICSI is never actually missing; it's n/a for all cycle types excepts Fresh
# mosaic::tally(any_icsi ~ cycle_typef, data=sart_forimp)
# mosaic::tally(fsh_gt10f ~ cycle_typef, data=sart_forimp) 
# mosaic::tally(fsh_gt10f ~ cycle_typef, data=sart_forimp, format="percent")
# mosaic::tally(amhcat3f ~ cycle_typef, data=sart_forimp)
# mosaic::tally(amhcat3f ~ cycle_typef, data=sart_forimp)
# mosaic::tally(amhcat3f ~ cycle_typef, data=sart_forimp, format="percent")
# mosaic::tally(pgd ~ cycle_typef, data=sart_forimp)
# mosaic::tally(pgd ~ cycle_typef, data=sart_forimp, format="percent")
#mosaic::favstats(num_transferred ~ cycle_typef, data=sart_forimp)
#sart_forimp %>% count(cycle_typef, pgd, num_transferred)

# include state in imputation model:
# https://stefvanbuuren.name/fimd/sec-missmult.html
# "If the primary interest is on the fixed effects, adding a cluster dummy is an
# easily implementable alternative, unless the missing rate is very large and/or 
# the intra-class correlation is very low and the number of records in the cluster
# is small (Drechsler 2015; L?dtke, Robitzsch, and Grund 2017)."

# ------------------------- set up imputation parameters ----------------------

# https://stefvanbuuren.name/fimd/sec-knowledge.html
# 6.4.2 Interaction terms

# update default imputation models
meth <- make.method(m0_forimp)
meth

# implement cart, which is more flexible re interactions, non-linearities etc.
# https://stefvanbuuren.name/fimd/sec-modelform.html
meth_cart <- meth
meth_cart[which(meth != "")] <- "cart"
meth_cart

# update default predictor matrix
pred <- make.predictorMatrix(m0_forimp)
pred

# don't allow patient age squared (at a current iteration) to predict 
# patient age (at the next iteration)
pred["patient_age_at_start", "patient_age_at_start_sq"] <- 0

# remove cancelled and sab as predictors since they'll be highly correlated w preg and lb
# remove patient and cycle id as predictors (only keeping in dataset so can merge later if need other vars)
pred[, c("mandate_class1", "cancelled", "sab_preg", "external_patient_id", "external_cycle_id")] <- 0
pred 

# -------------------------- perform imputations -------------------------------

# ------------------ default imputation method ---------------------------------

start_time1 <- Sys.time()
start_time1

imp_m0_def <- mice(data=m0_forimp, m=1, meth=meth, pred=pred
                   , maxit=20, seed=cantor)

save(imp_m0_def, file=paste0(results,"/",Cluster,".",Process,".imp_m0_def.RData"))

imp_m0_def$loggedEvents

imp_m1_def <- mice(data=m1_forimp, m=1, meth=meth, pred=pred
                   , maxit=20, seed=cantor)

save(imp_m1_def, file=paste0(results,"/",Cluster,".",Process,".imp_m1_def.RData"))

imp_m1_def$loggedEvents

end_time1 <- Sys.time()
end_time1

end_time1 - start_time1

# ------------------ impute all using CART ------------------------------------

start_time2 <- Sys.time()
start_time2

imp_m0_cart <- mice(data=m0_forimp, m=1, meth=meth_cart, pred=pred
                    , maxit=20, seed=cantor)

save(imp_m0_cart, file=paste0(results,"/",Cluster,".",Process,".imp_m0_cart.RData"))

imp_m0_cart$loggedEvents

imp_m1_cart <- mice(data=m1_forimp, m=1, meth=meth_cart, pred=pred
                    , maxit=20, seed=cantor)

save(imp_m1_cart, file=paste0(results,"/",Cluster,".",Process,".imp_m1_cart.RData"))

imp_m1_cart$loggedEvents

end_time2 <- Sys.time()
end_time2

end_time2 - start_time2


# ------------------------- save imputed datasets ------------------------------

imp_def_dat <- imp_m0_def %>%
  mice::complete("long") %>%
  bind_rows(imp_m1_def %>%
              mice::complete("long")) %>%
  mutate(imp_meth = "default")

imp_cart_dat <- imp_m0_cart %>%
  mice::complete("long") %>%
  bind_rows(imp_m1_cart %>%
              mice::complete("long")) %>%
  mutate(imp_meth = "cart")

imp_dat <- bind_rows(imp_def_dat, imp_cart_dat)

saveRDS(imp_dat, file=paste0(results,"/",Cluster,".",Process,".imp_dat.RDS"))
