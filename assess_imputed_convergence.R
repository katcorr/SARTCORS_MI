# assess mice convergence among 20 imputed datasets (each w 20 iter)

library(tidyverse)

# compare proportions in imputed dataset to proportions observed
sart_full <- readRDS("~/cluster-scratch/SARTCORS/sart_full.RDS")

jobnum <- "1726"
path <- paste0("~/cluster-scratch/SARTCORS/results", jobnum)
M <- 20

for (m in 1:(M-1)){
  imp_dat <- readRDS(paste0(path, "/", jobnum, ".", m, ".imp_dat.RDS")) %>%
    mutate(impnum=m+1)
  
  if (m==1){
    imp_dat_all <- imp_dat
  }
  if (m>1){
    imp_dat_all <- bind_rows(imp_dat_all, imp_dat)
  }
}

names(imp_dat_all)
count(imp_dat_all, impnum, imp_meth)

# since all imputed fields are categorical, compare proportions 
meth_vec <- c("default", "cart")

race_obs <- sart_full %>%
  filter(patient_race_ethnic != "Missing") %>%
  mutate(impnum=0) %>%
  select(impnum, race = patient_race_ethnic)

race_dat <- imp_dat_all %>%
  bind_rows(race_obs %>% mutate(imp_meth="default")) %>%
  bind_rows(race_obs %>% mutate(imp_meth="cart")) %>%
  group_by(imp_meth, impnum) %>%
  count(race) %>%
  mutate(race_perc = n/sum(n)) %>%
  ungroup()

for (i in meth_vec){
  
  dat <- race_dat %>%
    filter(imp_meth == i)
  
  plot1 <- ggplot(data=dat, aes(x=impnum, y=race_perc, fill=race)) +
    geom_bar(position="fill", stat="identity") +
    labs(title=i) 
  
  print(plot1)
  
  plot2 <- ggplot(data=dat, aes(x=impnum, y=race_perc, fill=race)) +
    geom_bar(stat="identity", color="black") +
    facet_wrap(~race, scales="free_y")
  
  print(plot2)
}

# # PAROUS
parous_obs <- sart_full %>%
  filter(!is.na(parous)) %>%
  mutate(impnum=0
         , parousf = as.factor(parous)) %>%
  select(impnum, parousf)

parous_dat <- imp_dat_all %>%
  bind_rows(parous_obs %>% mutate(imp_meth="default")) %>%
  bind_rows(parous_obs %>% mutate(imp_meth="cart")) %>%
  group_by(imp_meth, impnum) %>%
  count(parousf) %>%
  mutate(parous_perc = n/sum(n)) %>%
  ungroup()

for (i in meth_vec){
  
  dat <- parous_dat %>%
    filter(imp_meth == i)
  
  plot1 <- ggplot(data=dat, aes(x=impnum, y=parous_perc, fill=parousf)) +
    geom_bar(position="fill", stat="identity") +
    labs(title=i)
  
  print(plot1)
  
  plot2 <- ggplot(data=dat, aes(x=impnum, y=parous_perc, fill=parousf)) +
    geom_bar(stat="identity", color="black") +
    facet_wrap(~parousf, scales="free_y")
  
  print(plot2)
}


# AMH
amh_obs <- sart_full %>%
  mutate(amhcat3f = factor(case_when(amhcat3=="Missing" & cycle_type %in% c("Cancelled", "Fresh")
                                     ~ NA_character_
                                     , amhcat3=="Missing" & cycle_type %in% c("Embryo Banking", "Frozen")
                                     ~ "N/A"
                                     , TRUE ~ as.character(amhcat3))
                           , levels=c("1 - <4", "< 1", ">= 4", "N/A"))) %>%
  filter(!is.na(amhcat3f)) %>%
  mutate(impnum=0) %>%
  select(impnum, amhcat3f)

amh_dat <- imp_dat_all %>%
  bind_rows(amh_obs %>% mutate(imp_meth="default")) %>%
  bind_rows(amh_obs %>% mutate(imp_meth="cart")) %>%
  group_by(imp_meth, impnum) %>%
  count(amhcat3f) %>%
  mutate(amh_perc = n/sum(n)) %>%
  ungroup()

for (i in meth_vec){
  
  dat <- amh_dat %>%
    filter(imp_meth == i)
  
  plot1 <- ggplot(data=dat, aes(x=impnum, y=amh_perc, fill=amhcat3f)) +
    geom_bar(position="fill", stat="identity") +
    labs(title=i)
  
  print(plot1)
  
  plot2 <- ggplot(data=dat, aes(x=impnum, y=amh_perc, fill=amhcat3f)) +
    geom_bar(stat="identity", color="black") +
    facet_wrap(~amhcat3f, scales="free_y")
  
  print(plot2)
}


# # ANY PRIOR SAB 
sab_obs <- sart_full %>%
  filter(!is.na(any_prior_sab)) %>%
  mutate(impnum=0
         , any_prior_sabf = as.factor(any_prior_sab)) %>%
  select(impnum, any_prior_sabf)

sab_dat <- imp_dat_all %>%
  bind_rows(sab_obs %>% mutate(imp_meth="default")) %>%
  bind_rows(sab_obs %>% mutate(imp_meth="cart")) %>%
  group_by(imp_meth, impnum) %>%
  count(any_prior_sabf) %>%
  mutate(sab_perc = n/sum(n)) %>%
  ungroup()

for (i in meth_vec){
  
  dat <- sab_dat %>%
    filter(imp_meth == i)
  
  plot1 <- ggplot(data=dat, aes(x=impnum, y=sab_perc, fill=any_prior_sabf)) +
    geom_bar(position="fill", stat="identity") +
    labs(title=i)
  
  print(plot1)
  
  plot2 <- ggplot(data=dat, aes(x=impnum, y=sab_perc, fill=any_prior_sabf)) +
    geom_bar(stat="identity", color="black") +
    facet_wrap(~any_prior_sabf, scales="free_y")
  
  print(plot2)
}


# # BMI
bmi_obs <- sart_full %>%
  filter(bmicat4 != "Missing") %>%
  mutate(impnum=0
         , bmicat4f = as.factor(bmicat4)) %>%
  select(impnum, bmicat4f)

bmi_dat <- imp_dat_all %>%
  bind_rows(bmi_obs %>% mutate(imp_meth="default")) %>%
  bind_rows(bmi_obs %>% mutate(imp_meth="cart")) %>%
  group_by(imp_meth, impnum) %>%
  count(bmicat4f) %>%
  mutate(bmi_perc = n/sum(n)) %>%
  ungroup()

for (i in meth_vec){
  
  dat <- bmi_dat %>%
    filter(imp_meth == i)
  
  plot1 <- ggplot(data=dat, aes(x=impnum, y=bmi_perc, fill=bmicat4f)) +
    geom_bar(position="fill", stat="identity") +
    labs(title=i)
  
  print(plot1)
  
  plot2 <- ggplot(data=dat, aes(x=impnum, y=bmi_perc, fill=bmicat4f)) +
    geom_bar(stat="identity", color="black") +
    facet_wrap(~bmicat4f, scales="free_y")
  
  print(plot2)
}


# # FSH
fsh_obs <- sart_full %>%
  filter(fsh_gt10 != "Missing") %>%
  mutate(impnum=0
         , fsh_gt10f = as.factor(fsh_gt10)) %>%
  select(impnum, fsh_gt10f)

fsh_dat <- imp_dat_all %>%
  bind_rows(fsh_obs %>% mutate(imp_meth="default")) %>%
  bind_rows(fsh_obs %>% mutate(imp_meth="cart")) %>%
  group_by(imp_meth, impnum) %>%
  count(fsh_gt10f) %>%
  mutate(fsh_perc = n/sum(n)) %>%
  ungroup()

for (i in meth_vec){
  
  dat <- fsh_dat %>%
    filter(imp_meth == i)
  
  plot1 <- ggplot(data=dat, aes(x=impnum, y=fsh_perc, fill=fsh_gt10f)) +
    geom_bar(position="fill", stat="identity") +
    labs(title=i)
  
  print(plot1)
  
  plot2 <- ggplot(data=dat, aes(x=impnum, y=fsh_perc, fill=fsh_gt10f)) +
    geom_bar(stat="identity", color="black") +
    facet_wrap(~fsh_gt10f, scales="free_y")
  
  print(plot2)
}

rm(list = ls(all.names = TRUE))

#quit(save="no")

# ------------------------------------------------------------------------------
# ----------------      STANDARD CONVERGENCE PLOTS  ----------------------------  
# ------------------------------------------------------------------------------

# since each field is a factor, it's taking the factor value average
# the plots above seem more relevant

# cleared environment above, so need to redefine these parameters
jobnum <- "1726"
path <- paste0("~/cluster-scratch/SARTCORS/results", jobnum)
M <- 20
meth_vec <- c("def", "cart")

for (i in meth_vec){
  
  print(i)
  all_dat <- data.frame()
  
  # Mset
  for (m in 0:(M-1)){
    
    load(paste0(path,"/", jobnum,".",m,".imp_m0_", i,".RData"))
    load(paste0(path,"/", jobnum,".",m,".imp_m1_", i,".RData"))
    
    tempdat0 <- data.frame(eval(as.name(paste0("imp_m0_", i)))$chainMean) %>%
      rownames_to_column() %>%
      filter((str_detect(rowname, "race") & !str_detect(rowname, "mandate_race")) | 
               (rowname %in% c("bmicat4f", "fsh_gt10f", "amhcat3f"
                               , "parousf", "any_prior_sabf"))) %>%
      pivot_longer(cols=-rowname, names_to = "iteration00", values_to = "value") %>%
      separate(iteration00, into=c("iteration0"), remove=FALSE) %>%
      mutate(iteration = parse_number(iteration0)
             , impnum = m + 1
             , imp_meth = i
             , mandate = 0)
    
    tempdat1 <- data.frame(eval(as.name(paste0("imp_m1_", i)))$chainMean) %>%
      rownames_to_column() %>%
      filter((str_detect(rowname, "race") & !str_detect(rowname, "mandate_race")) | 
               (rowname %in% c("bmicat4f", "fsh_gt10f", "amhcat3f"
                               , "parousf", "any_prior_sabf"))) %>%
      pivot_longer(cols=-rowname, names_to = "iteration00", values_to = "value") %>%
      separate(iteration00, into=c("iteration0"), remove=FALSE) %>%
      mutate(iteration = parse_number(iteration0)
             , impnum = m + 1
             , imp_meth = i
             , mandate = 1)
    
    all_dat <- bind_rows(all_dat, tempdat0, tempdat1)
    rm(tempdat0, tempdat1)
    rm(list=ls(pattern="^imp_"))
  }
  
  plot0 <- all_dat %>%
    filter(mandate==0) %>%
    ggplot(aes(x=iteration, y=value, group=impnum
               , color=as.factor(impnum))) +
    geom_line() +
    labs(title=paste0(i,": No mandate")) +
    facet_wrap(~rowname, scales="free")
  
  print(plot0)
  
  plot1 <- all_dat %>%
    filter(mandate==1) %>%
    ggplot(aes(x=iteration, y=value, group=impnum
               , color=as.factor(impnum))) +
    geom_line() +
    labs(title=paste0(i, ": Mandate")) +
    facet_wrap(~rowname, scales="free")
  
  print(plot1)
  
  rm(all_dat)
}