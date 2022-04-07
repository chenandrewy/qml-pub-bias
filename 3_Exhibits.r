# 2022 03 makes exhibits from bootall, the output of 2_CustomStats
rm(list = ls())
# ENVIRONMENT ====

source('0_Settings.r')

load('intermediate/boot deleteme started 2022-04-07 15-55.Rdata')


bootall.wide = bootpar %>% 
  left_join(bootstat, by = 'booti') %>% 
  as_tibble()

bootall.long = bootall.wide %>% 
  select(-c(mufam,pubfam)) %>% 
  pivot_longer(
    cols = c(-booti)
  ) %>% 
  rename(stat = name)


# TABLES ====

## Simple tables ====
stat_order = tibble(
  stat = c('pif','pia','mua','siga','mub','sigb','pubpar2')
)
stat_order$id = 1:nrow(stat_order)

# point estimate
point = bootall.long %>% 
  filter(booti == 1) %>% 
  inner_join(stat_order) %>% 
  select(-booti) %>% 
  rename(point = value) 

# bootstrap order statistics
bootsum = bootall.long %>% 
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


write.csv(tab.est, 'output/tab-est-mu-sort.csv', row.names = F)

# PLOT SKETCHES ====

## Model fit ====

# the first boot is the point estimate
est.point = bootall.wide[1,]
  mutate(
    pubfam = set.boot$model_fam$pubfam[1]
    , mufam = set.boot$model_fam$mufam[1]
  )

# call function for plotting point estimate (in 0_Environment)
hist_emp_vs_par(cz_filt$tabs,est.point)



## Prop 1: Panel A ====
bootall.long %>%
  filter(stat %in% c('pif','pr_tgt_2')) %>% 
  ggplot(
    aes(x=value, fill=stat)
  ) +
  geom_histogram(position = 'dodge', alpha = 0.6)

## Prop 1: Panel B ====
bootall.long %>% 
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

bootall.long %>% 
  filter(stat %in% c('pif','h_fdr5')) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(
    pif_cat = findInterval(pif, pif_cut)
    , pif_cat = as.factor(pif_cat)
  ) %>% 
  ggplot(
    aes(x=h_fdr5, fill=pif_cat)
  ) +
  geom_histogram(breaks = breaks) +
  geom_vline(
    xintercept = qnorm(1-0.05/2)
  )

## t-hurdle: fdr = 1% ====

pif_cut = c(1/3,2/3)
breaks = seq(0,4,0.2)

bootall.long %>% 
  filter(stat %in% c('pif','h_fdr1')) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(
    pif_cat = findInterval(pif, pif_cut)
    , pif_cat = as.factor(pif_cat)
  ) %>% 
  ggplot(
    aes(x=h_fdr1, fill=pif_cat)
  ) +
  geom_histogram(breaks = breaks) +
  geom_vline(
    xintercept = qnorm(1-0.01/2)
  )

