# 2022 04 09 Robustness Exhibits

rm(list = ls())
# SETUP ====

source('0_Settings.r')

# choose files to import
rdatalist = dir('intermediate/', full.names = T) %>% 
  as_tibble() %>% 
  filter(grepl('boot robust-', value) | 
grepl('boot a-priori-model simple', value)) %>% 
  pull(value)


# IMPORT STATS ====
nspec = length(rdatalist)

tab.all = tibble()
for (speci in 1:nspec){  
  
  # merge bootstrapped parameters and boot statistics in two formats
  # checking to see if it's the old (pre 2023-09) or new format
  rdata_items = load(rdatalist[speci])
  if (length(rdata_items) > 1){
    # old style
    load(rdatalist[speci])
    bootall.wide = bootpar %>% 
      left_join(bootstat, by = 'booti') %>% 
      as_tibble()
    name_cur = bootname
  }  else {    
    # new style
    load(rdatalist[speci])
    bootall.wide = boot$par %>% 
      left_join(boot$stat, by = 'booti') %>% 
      as_tibble()    
    name_cur = boot$name
  }
  
  # clean up, multiply 100 if it's nice
  bootall.long =  bootall.wide %>% 
    select(-c(mufam,pubfam)) %>% 
    pivot_longer(
      cols = c(-booti)
    ) %>% 
    rename(stat = name) %>% 
    mutate(
      value = if_else(grepl('bias',stat), value*100, value)
      , value = if_else(grepl('fdrloc',stat), value*100, value)
    )
  
  # find point estimates
  tab.point = bootall.long %>% 
    filter(booti == 1)  %>% 
    select(-booti) %>% 
    rename(point = value) 
 
  # find mean, med, standard errors
  tab.mstat = bootall.long %>% 
    group_by(stat) %>% 
    summarize(
      sd = sd(value), mean = mean(value), med = median(value)
      , p10 = quantile(value, 0.10), p90 = quantile(value, 0.90)
      , p05 = quantile(value, 0.05), p95 = quantile(value, 0.95)      
      , nboot = dplyr::n()
    )
  
  # merge
  tab.both = tab.point %>% 
    left_join(tab.mstat, by = 'stat') %>% 
    mutate(name = name_cur) %>% 
    select(name, stat, point, sd, mean, med, p05, p10, p90, p95, nboot) 
  
  
  # append
  tab.all = tab.all %>% rbind(tab.both)
  
} # for speci

# check for duplicates
namesdistinct = tab.all %>% distinct(name) %>% pull()
if (length(namesdistinct) != length(rdatalist)){
  rdatalist %>% print()
  stop('duplicate robustness tests')
}


# CLEAN UP ====

mstat1 = 'p05'
mstat2 = 'p95'

# select statistics and order
statlist = tibble(
  stat = c('h_fdr5','h_fdr1','bias_mean','fdrloc_mean','sbias_mean','nboot')
) 
statlist$statid = 1:dim(statlist)[1]
tab.small = tab.all %>% 
  inner_join(statlist, by = 'stat') 

# change to wide
tab.wide = tab.small %>% 
  select(name, statid, stat, all_of(mstat1), all_of(mstat2)) %>% 
  pivot_longer(
    cols = c(mstat1, mstat2), names_to = 'mstat'
  ) %>% 
  arrange(name,statid,mstat) %>% 
  select(-statid) %>% 
  pivot_wider(
    names_from = c(stat,mstat)
  )  %>% 
  left_join(
    tab.small  %>% distinct(name, nboot), by = 'name'
  )

# use signed bias (sbias) if mutrue can be negative   
tab.wide2 = tab.wide %>% 
  mutate(
    bias_mean_p05 = if_else(
      name %in% c('robust-tdist','robust-mix_norm'), sbias_mean_p05, bias_mean_p05
    )  
    , bias_mean_p95 = if_else(
      name %in% c('robust-tdist','robust-mix_norm'), sbias_mean_p95, bias_mean_p95
    )  
  )  %>% 
  select(name
  , starts_with('h_fdr'), starts_with('bias_mean'), starts_with('fdrloc')
  , nboot
  )

# rename 
tab.wide2 = tab.wide2  %>% 
  mutate(
    name = str_remove(name, 'robust-')
    , name = if_else(name == 'a-priori-model simple', 'baseline', name)
  )

rownames(tab.wide2) = tab.wide2$name

# reorder
tab.sorted = tab.wide2[
  c(
    'baseline'
    , 'raw-stair2'
    , 'raw-stair2-restrict'
    , 'pif_10'
    , 'pif_20'
    , 'pubpar_05'
    , 'relax_pubpar'    
    , 'logistic'        
    , 'exp'
    , 'tdist'
    , 'mix_norm'
  )
  , 
] %>% 
  add_column(blank1 = NA, .after = 3)  %>% 
  add_column(blank2 = NA, .after = 6) %>% 
  add_column(blank3 = NA, .after = 9)

tab.sorted %>% print()

# write to disk
write.csv(tab.sorted, 'output/tab-robust.csv', row.names = F)