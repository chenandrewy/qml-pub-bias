# 2022 04 09 Robustness Exhibits

rm(list = ls())
# SETUP ====

source('0_Settings.r')

# choose files to import
rdatalist = c(
  'intermediate/boot robust-relax_pubpar started 2022-04-09 00-23.Rdata'
  ,'intermediate/boot robust-pif_20 started 2022-04-09 01-14.Rdata'
  ,'intermediate/boot robust-exp started 2022-04-09 02-11.Rdata'
  ,'intermediate/boot robust-norm started 2022-04-09 03-04.Rdata'
)

# IMPORT STATS ====
nspec = length(rdatalist)

tab.all = tibble()
for (speci in 1:nspec){
  load(rdatalist[speci])
  
  # rename tabs to be more clear (these are |t| from filtered cz data)
  cz_filt_tabs = tabs
  
  # merge bootstrapped parameters and boot statistics in two formats
  bootall.wide = bootpar %>% 
    left_join(bootstat, by = 'booti') %>% 
    as_tibble()
  
  bootall.long = bootall.wide %>% 
    select(-c(mufam,pubfam)) %>% 
    pivot_longer(
      cols = c(-booti)
    ) %>% 
    rename(stat = name)
  
  # find point estimates
  tab.point = bootall.long %>% 
    filter(booti == 1)  %>% 
    select(-booti) %>% 
    rename(point = value) 
  
  # find standard errors
  tab.se = bootall.long %>% 
    group_by(stat) %>% 
    summarize(
      sd = sd(value)
    )
  
  # merge
  tab.both = tab.point %>% 
    left_join(tab.se, by = 'stat') %>% 
    mutate(name = bootname) %>% 
    select(name, stat, point, sd) 
  
  
  # append
  tab.all = tab.all %>% rbind(tab.both)
  
} # for speci


# CLEAN UP ====

statlist = tibble(
  stat = c('pif','mua','siga','nua','pubpar1','h_fdr5','bias_mean')
) 
statlist$id = 1:dim(statlist)[1]


tab.small = tab.all %>% 
  inner_join(statlist, by = 'stat') %>% select(-id) 


tab.wide = rbind(
tab.small %>% select(-sd) %>% 
  pivot_wider(
    names_from = stat, values_from = c(point)
  ) %>% 
  mutate(statname = 'point', statid = 1)
  
, tab.small %>% select(-point) %>% 
  pivot_wider(
    names_from = stat, values_from = c(sd)
  ) %>% 
  mutate(statname = 'sd', statid = 2)

) %>% 
  arrange(name,statid) %>% 
  select(name,statname, everything(), -statid)

tab.wide 