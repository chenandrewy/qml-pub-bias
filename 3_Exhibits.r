# ENVIRONMENT ====

source('0_Environment.r')

pubcross = fread('exhibits/pubcross.csv') 
bootall = fread('exhibits/bootall.csv')

# PLOT SKETCHES ====

## Model fit ====

# the first boot is the point estimate
est.point = bootall %>% 
  filter(booti == 1) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(
    pubfam = 'stair'# ugly hack, i'm sorry
  )

# call function for plotting point estimate (in 0_Environment)
hist_emp_vs_par(est.point)



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
  geom_histogram() +
  geom_vline(
    xintercept = qnorm(1-0.05/2)
  )

## t-hurdle: fdr = 1% ====

pif_cut = c(1/3,2/3)

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
  geom_histogram() +
  geom_vline(
    xintercept = qnorm(1-0.01/2)
  )