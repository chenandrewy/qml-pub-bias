# 2022 04 09 Robustness Exhibits

rm(list = ls())
# SETUP ====

source('0_Settings.r')

# choose files to import
rdatalist = dir('intermediate/', full.names = T) %>% 
  as_tibble() %>% 
  filter(grepl('boot robust-', value)) %>% 
  pull(value)


# IMPORT STATS ====
nspec = length(rdatalist)

tab.all = tibble()
for (speci in 1:nspec){
  
  # speci = 8
  load(rdatalist[speci])
  
  # merge bootstrapped parameters and boot statistics in two formats
  bootall.wide = bootpar %>% 
    left_join(bootstat, by = 'booti') %>% 
    as_tibble()
  
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
    )
  
  # merge
  tab.both = tab.point %>% 
    left_join(tab.mstat, by = 'stat') %>% 
    mutate(name = bootname) %>% 
    select(name, stat, point, sd, mean, med) 
  
  
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

# select statistics and order
statlist = tibble(
  stat = c('h_fdr5','h_fdr1','bias_mean','fdrloc_mean')
) 
statlist$statid = 1:dim(statlist)[1]


tab.small = tab.all %>% 
  mutate(
    name = substr(name, 8, nchar(name))
  ) %>% 
  inner_join(statlist, by = 'stat') 


tab.wide = tab.small %>% 
  select(name, statid, stat, med, sd) %>% 
  pivot_longer(
    cols = c(med,sd), names_to = 'mstat'
  ) %>% 
  arrange(name,statid,mstat) %>% 
  select(-statid) %>% 
  pivot_wider(
    names_from = c(stat,mstat)
  )

rownames(tab.wide) = tab.wide$name

# reorder
tab.sorted = tab.wide[
  c(
    'raw-stair2'
    , 'raw-stair3'
    , 'raw-logistic'    
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

