# fit all models for clinical outcomes for MULTIPLE IMPUTATION ANALYSIS
# and save results
# estimating the ORs using the GEE-1 approach because the gee functions in
# established R packages (gee, geepack) crashes the system due to the large
# cluster sizes

args = commandArgs(trailingOnly=TRUE)
print(args)

Cluster = args[1]
Process = args[2]
# method is cart or default
method = args[3]

print(method)

library(tidyverse)
 
# GEE-1 OR function
source("~/cluster-scratch/SARTCORS/onestep_gee_Rfunction.R")  

# full dataset to get lb and preg and other fields for CLB dataset (linked external ids, etc.)
sart_full <- readRDS("sart_full.RDS") %>%
  select(reporting_year, external_patient_id, external_cycle_id
         , starts_with("linked_"), starts_with("source_")
         , lb, lb_et, preg, preg_et, sab_preg)

# --------------------------  GET IMPUTED DATASETS -----------------------------

# multiply-imputed datasets
jobnum <- "1726"
path <- paste0("results",jobnum,"/")

M <- 20
Mset <- c(0:(M-1))

for (m in Mset){
  assign(paste0("imp",m+1), readRDS(paste0(path,jobnum,".",m,".imp_dat.RDS")) %>%
                            filter(imp_meth==method) %>%
                            mutate(imp=m+1
                                   , numET4 = as.factor(case_when(num_transferred %in% c(0,1,2,3) ~ as.character(num_transferred)
                                                      , num_transferred >= 4 ~ "4+"))) %>%
                            rename(patient_race_ethnic = race
                                  , state = statef
                                  , parous = parousf
                                  , fsh_gt10 = fsh_gt10f
                                  , amhcat3 = amhcat3f
                                  , any_prior_sab = any_prior_sabf
                                  , bmicat4 = bmicat4f
                                  , frozen = frozenf) %>%
                            select(-sab_preg) %>%
                            left_join(sart_full, by=c("reporting_year", "external_patient_id", "external_cycle_id"))
         )
  
  # compute dataset for clb
  fresh_stim <- eval(as.name(paste0("imp",m+1))) %>%
    filter(reporting_year %in% c(2016:2019) & frozen != 1) %>%
    select(external_patient_id, reporting_year, external_cycle_id, 
           starts_with("linked_"), frozen, lb, everything())
  
  thaws <- eval(as.name(paste0("imp",m+1))) %>%
    filter(reporting_year %in% c(2016:2019)) %>%
    filter(linked_source_external_cycle_i_ds != "NULL") %>%
    separate(linked_source_external_cycle_i_ds, into=c(paste0("source",1:4))
             , sep=",", remove=FALSE, convert=TRUE) %>%
    select(external_patient_id, reporting_year, external_cycle_id
           , starts_with("source"), starts_with("linked_")
           , frozen, lb, everything())

  thaws_long <- thaws %>%
    select(external_patient_id, reporting_year, external_cycle_id
           , starts_with("source"), lb) %>%
    pivot_longer(cols=c(starts_with("source")), names_to="name", values_to="source_cycle_id") %>%
    filter(!is.na(source_cycle_id))
  
  clb0 <- fresh_stim %>%
    select(external_patient_id, reporting_year, external_cycle_id, lb) %>%
    mutate(source_cycle_id = external_cycle_id) %>%
    bind_rows(thaws_long) %>%
    group_by(external_patient_id, source_cycle_id) %>%
    summarize(clb0 = max(lb, na.rm=TRUE), ncyc= n()) %>%
    ungroup() %>%
    # -Inf has been assigned to cycles with no intent to transfer (NA for lb)
    # and no following cycles that used those embryos; set clb to NA 
    # (exclude from clb analysis)
    mutate(clb = ifelse(is.infinite(clb0), yes=NA_real_, no=clb0)) %>%
    # remove frozen cycles with stim cycles before 2016
    filter(source_cycle_id %in% unique(fresh_stim$external_cycle_id))
    
  # add potential covariates to dataset (all are related to stim cycle)
  clb <- clb0 %>%
    left_join(fresh_stim, by=c("external_patient_id", "source_cycle_id"="external_cycle_id")) %>%
    select(reporting_year, external_patient_id, source_cycle_id, clb, lb, everything())

  assign(paste0("imp",m+1,"_clb"), clb)
  
  rm(fresh_stim, thaws, thaws_long, clb0, clb)
}

# imp1 %>% count(lb)
# imp1 %>% count(preg)
# imp1 %>% count(sab_preg)
# imp1 %>% count(cancelled)
# imp1 %>% count(frozen)

# imp1_clb %>% count(imp_meth)
# mosaic::tally(clb ~ patient_race_ethnic, data=imp1_clb, format="percent")
# mosaic::tally(clb ~ patient_race_ethnic, data=imp10_clb, format="percent")

# # ------------------------  FIT MODELS -----------------------------------
#
# outcomes: clb (cumulative live birth), lb, lb per ET, preg, preg per ET, sab, cancelled
# # a models include mandate only;
# # b models include race/ethnicity only;
# # c models include mandate, race, and other covariates
# # d models include mandate, race, and other covariates plus number ET, ICSI, and frozen cycle 
# # e models include, mandate, race, the interaction between them, and other covariates
# # f models include mandate, race, the interaction between them, and other covariates plus number ET, ICSI, and frozen cycle

pred_adj1 <- "mandate_class1 + patient_race_ethnic +  patient_age_at_start + I(patient_age_at_start^2) + bmicat4 + parous + any_prior_sab + fsh_gt10 + amhcat3 + male_infertility + endometriosis + dx_tubal + uterine + unexplained + dx_ovulation + diminished_ovarian_reserve"
pred_adj1

pred_adj2 <- paste0(pred_adj1, " + numET4 + any_icsi + frozen")
pred_adj2

pred_int1 <- paste0(pred_adj1, "+ mandate_class1*patient_race_ethnic")
pred_int1

pred_int2 <- paste0(pred_adj2, "+ mandate_class1*patient_race_ethnic")
pred_int2

outcomes_vec <- c("clb", "lb", "lb_et", "preg", "preg_et", "sab_preg", "cancelled")


modsOR_mi_list <- list()

for (m in 1:M){

  print(m)
  
  for (i in outcomes_vec){

    print(i)
    
    if (i=="clb"){
      dat <- paste0("imp",m, "_clb")
      pred_adj2 <- str_replace_all(pred_adj2, "\\+ frozen", "")
      pred_adj2
      pred_int2 <- str_replace_all(pred_int2, "\\+ frozen", "")
      pred_int2
    }
    if (i!="clb"){
      dat <- paste0("imp",m)
    }

    ll <- length(modsOR_mi_list)

    modsOR_mi_list[[c(ll+1)]] <- onestepGEE(dataset=eval(as.name(dat))
                                             , outcome=i, cluster=state
                                             , predictors = "mandate_class1")

    modsOR_mi_list[[c(ll+1)]]$outcome <- i
    modsOR_mi_list[[c(ll+1)]]$dataset <- dat
    modsOR_mi_list[[c(ll+1)]]$imp_meth <- method

    modsOR_mi_list[[c(ll+2)]] <- onestepGEE(dataset=eval(as.name(dat))
                                             , outcome=i, cluster=state
                                             , predictors = "patient_race_ethnic")

    modsOR_mi_list[[c(ll+2)]]$outcome <- i
    modsOR_mi_list[[c(ll+2)]]$dataset <- dat
    modsOR_mi_list[[c(ll+2)]]$imp_meth <- method

    modsOR_mi_list[[c(ll+3)]] <- onestepGEE(dataset=eval(as.name(dat))
                                             , outcome=i, cluster=state
                                             , predictors = pred_adj1)

    modsOR_mi_list[[c(ll+3)]]$outcome <- i
    modsOR_mi_list[[c(ll+3)]]$dataset <- dat
    modsOR_mi_list[[c(ll+3)]]$imp_meth <- method

    if (i != "cancelled"){
      modsOR_mi_list[[c(ll+4)]] <- onestepGEE(dataset=eval(as.name(dat))
                                               , outcome=i, cluster=state
                                               , predictors = pred_adj2)
      modsOR_mi_list[[c(ll+4)]]$outcome <- i
      modsOR_mi_list[[c(ll+4)]]$dataset <- dat
      modsOR_mi_list[[c(ll+4)]]$imp_meth <- method

      modsOR_mi_list[[c(ll+5)]] <- onestepGEE(dataset=eval(as.name(dat))
                                               , outcome=i, cluster=state
                                               , predictors = pred_int1)

      modsOR_mi_list[[c(ll+5)]]$outcome <- i
      modsOR_mi_list[[c(ll+5)]]$dataset <- dat
      modsOR_mi_list[[c(ll+5)]]$imp_meth <- method


      modsOR_mi_list[[c(ll+6)]] <- onestepGEE(dataset=eval(as.name(dat))
                                               , outcome=i, cluster=state
                                               , predictors = pred_int2)

      modsOR_mi_list[[c(ll+6)]]$outcome <- i
      modsOR_mi_list[[c(ll+6)]]$dataset <- dat
      modsOR_mi_list[[c(ll+6)]]$imp_meth <- method

    }

    if (i == "cancelled"){
      modsOR_mi_list[[c(ll+4)]] <- onestepGEE(dataset=eval(as.name(dat))
                                               , outcome=i, cluster=state
                                               , predictors = pred_int1)

      modsOR_mi_list[[c(ll+4)]]$outcome <- i
      modsOR_mi_list[[c(ll+4)]]$dataset <- dat
      modsOR_mi_list[[c(ll+4)]]$imp_meth <- method

    }

    # SAVE PERIODOCIALLY/AFTER EACH OUTCOME WITHIN ONE IMPUTED DATASET
    saveRDS(modsOR_mi_list, file=paste0("~/cluster-scratch/SARTCORS/results",jobnum,"/modsOR_mi_list_",Process,".RDS"))
  }
}

# ------------------- extract most relevant and save as dataframe -------------

for (k in 1:length(modsOR_mi_list)){
  df <- data.frame(analysis = "multiple imputation"
                   , dat =  modsOR_mi_list[[k]]$dataset
                   , imp_meth = modsOR_mi_list[[k]]$imp_meth
                   , outcome = modsOR_mi_list[[k]]$outcome
                   , predictors_all =  modsOR_mi_list[[k]]$predictors
                   , beta_hat_gee1 = modsOR_mi_list[[k]]$beta_hat_gee1
                   , se_gee1 =  modsOR_mi_list[[k]]$se_gee1) %>%
    rownames_to_column()

  if (k==1){
    modsOR_mi_df <- df
  }
  if (k>1){
    modsOR_mi_df <- bind_rows(modsOR_mi_df, df)
  }
}

saveRDS(modsOR_mi_df, file=paste0("~/cluster-scratch/SARTCORS/results",jobnum,"/modsOR_mi_df_",Process,".RDS"))

