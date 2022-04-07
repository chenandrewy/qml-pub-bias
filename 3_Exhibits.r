# 2022 03 makes exhibits from bootall, the output of 2_CustomStats
rm(list = ls())
# ENVIRONMENT ====

source('0_Settings.r')

load('intermediate/boot single-norm-trunc started 2022-04-06 14-46.Rdata')


bootall = bootpar
bootall$booti = row.names(bootpar)

bootall = bootall %>% 
  select(-c(pubfam,mufam)) %>% 
  pivot_longer(
    -booti, names_to = 'stat', values_to = 'value'
  ) %>% 
  rbind(
    bootstat
  )


# TABLES ====

## Simple tables ====
stat_order = tibble(
  stat = c('pif','pia','mua','siga','mub','sigb','pubpar2')
)
stat_order$id = 1:nrow(stat_order)

# point estimate
point = bootall %>% 
  filter(booti == 1) %>% 
  inner_join(stat_order) %>% 
  select(-booti) %>% 
  rename(point = value) 

# bootstrap order statistics
bootsum = bootall %>% 
  inner_join(stat_order) %>% 
  group_by(id, stat) %>% 
  summarize(
    p05 = quantile(value, .05)
    , p25 = quantile(value, .25)
    , p50 = quantile(value, .50)
    , p75 = quantile(value, .75)
    , p95 = quantile(value, .95)
  )

tab.est = point %>% 
  left_join(bootsum) %>% 
  arrange(id) %>% 
  select(-id)

tab.est

write.csv(tab.est, 'output/tab-est.csv', row.names = F)

## Tables with mu_small and mu_large ====
stat_order = tibble(
  stat = c('pif','pi_small','mu_small','sig_small','mu_large','sig_large','pubpar2')
)
stat_order$id = 1:nrow(stat_order)


bootall2= bootall %>% 
  pivot_wider(
    names_from = stat, values_from = value
  ) %>% 
  mutate(
    pi_small = if_else(mua<mub, pia, (1-pia))
    , mu_small = if_else(mua<mub, mua, mub)
    , sig_small = if_else(mua<mub, siga, sigb)    
    , mu_large = if_else(mua<mub, mub, mua)
    , sig_large = if_else(mua<mub, sigb, siga)
  ) %>% 
  select(
    -c(pia,mua,siga,mub,sigb)
  ) %>% 
  pivot_longer(
    cols = -booti, names_to = 'stat', values_to = 'value'
  )

# point estimate
point = bootall2 %>% 
  filter(booti == 1) %>% 
  inner_join(stat_order) %>% 
  select(-booti) %>% 
  rename(point = value) 

bootsum = bootall2 %>% 
  inner_join(stat_order) %>% 
  group_by(id, stat) %>% 
  summarize(
    blank1 = NA
    , p05 = quantile(value, .05)
    , blank2 = NA
    , p25 = quantile(value, .25)
    , p50 = quantile(value, .50)
    , p75 = quantile(value, .75)
    , blank3 = NA    
    , p95 = quantile(value, .95)
  )

tab.est2 = point %>% 
  left_join(bootsum) %>% 
  arrange(id) %>% 
  select(-id)

tab.est2

write.csv(tab.est2, 'output/tab-est-mu-sort.csv', row.names = F)

# PLOT SKETCHES ====

## Model fit ====

# the first boot is the point estimate
est.point = bootall %>% 
  filter(booti == 1) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(
    pubfam = set.boot$model_fam$pubfam[1]
    , mufam = set.boot$model_fam$mufam[1]
  )

# call function for plotting point estimate (in 0_Environment)
hist_emp_vs_par(cz_filt$tabs,est.point)



## Prop 1: Panel A ====
bootall %>% 
  filter(stat %in% c('pif','pr_tgt_2')) %>% 
  ggplot(
    aes(x=value, fill=stat)
  ) +
  geom_histogram(position = 'dodge', alpha = 0.6)

## Prop 1: Panel B ====
bootall %>% 
  filter(stat %in% c('pif','pr_tgt_2')) %>% 
  pivot_wider(names_from = stat, values_from = value) %>% 
  mutate(
    dif = pif - pr_tgt_2
  ) %>% 
  ggplot(
    aes(x=dif)
  ) + 
  geom_histogram() +
  geom_vline(xintercept = 0) +
  annotate(geom='text', x=0.5, y=100, label = 'raise hurdle') + 
  annotate(geom='text', x=-0.5, y=100, label = 'lower hurdle')

## t-hurdle: fdr = 5% ====

pif_cut = c(1/3,2/3)
breaks = seq(0,3.6,0.2)

bootall %>% 
  filter(stat %in% c('pif','fdr5')) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(
    pif_cat = findInterval(pif, pif_cut)
    , pif_cat = as.factor(pif_cat)
  ) %>% 
  ggplot(
    aes(x=fdr5, fill=pif_cat)
  ) +
  geom_histogram(breaks = breaks) +
  geom_vline(
    xintercept = qnorm(1-0.05/2)
  )

## t-hurdle: fdr = 1% ====

pif_cut = c(1/3,2/3)
breaks = seq(0,4,0.2)

bootall %>% 
  filter(stat %in% c('pif','fdr1')) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(
    pif_cat = findInterval(pif, pif_cut)
    , pif_cat = as.factor(pif_cat)
  ) %>% 
  ggplot(
    aes(x=fdr1, fill=pif_cat)
  ) +
  geom_histogram(breaks = breaks) +
  geom_vline(
    xintercept = qnorm(1-0.01/2)
  )

