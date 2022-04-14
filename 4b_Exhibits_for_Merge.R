# SETUP ====

rm(list = ls())
source('0_Settings.r')

# load semi-par bootstrap for robustness
load('intermediate/boot a-priori-model semipar started 2022-04-11 08.Rdata')

# merge bootstrapped parameters and boot statistics in two formats
bootalt.wide = bootpar %>% 
  mutate(
    exp_mu = exp(mua + siga^2/2)
    , exp_sd = sqrt( (exp(siga^2)-1)*exp(2*mua+siga^2) )
    , exp_med = exp(mua)
  ) %>% 
  left_join(bootstat, by = 'booti') %>% 
  as_tibble()

bootalt.long = bootalt.wide %>% 
  select(-c(mufam,pubfam)) %>% 
  pivot_longer(
    cols = c(-booti)
  ) %>% 
  rename(stat = name)


# load main bootstrap
load('intermediate/boot a-priori-model simple started 2022-04-11 13.Rdata')

# rename tabs to be more clear (these are |t| from filtered cz data)
cz_filt_tabs = tabs

# merge bootstrapped parameters and boot statistics in two formats
bootall.wide = bootpar %>% 
  mutate(
    exp_mu = exp(mua + siga^2/2)
    , exp_sd = sqrt( (exp(siga^2)-1)*exp(2*mua+siga^2) )
    , exp_med = exp(mua)
  ) %>%   
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
  stat = c('pif','exp_mu','exp_sd','exp_med','pubpar1','mua','siga')
)
stat_order$id = 1:nrow(stat_order)

# point estimate
point = bootall.long %>% 
  filter(booti == 1) %>% 
  inner_join(stat_order) %>% 
  select(-booti) %>% 
  rename(point = value) %>% 
  mutate(blank0 = NA)

# bootstrap order statistics
bootsum = bootall.long %>% 
  inner_join(stat_order) %>% 
  group_by(id, stat) %>% 
  summarize(
    p05 = quantile(value, .05)
    , blank1 = NA
    , p25 = quantile(value, .25)
    , p50 = quantile(value, .50)
    , p75 = quantile(value, .75)
    , blank2 = NA    
    , p95 = quantile(value, .95)
  )

# alternative bootstrap
bootsumalt = bootalt.long %>% 
  inner_join(stat_order) %>% 
  group_by(id, stat) %>% 
  summarize(
    p05 = quantile(value, .05)
    , blank1 = NA
    , p25 = quantile(value, .25)
    , p50 = quantile(value, .50)
    , p75 = quantile(value, .75)
    , blank2 = NA    
    , p95 = quantile(value, .95)
  ) %>% 
  mutate(point = NA, blank0 = NA)


tab.est = point %>% 
  left_join(bootsum) %>% 
  rbind(bootsumalt) %>% 
  arrange(id,point) %>% 
  select(-id)

tab.est

write.csv(tab.est, 'output/tab-est.csv', row.names = F)

# sEMIPAR BOOT ====


# simulate residuals
set.seed = 1253
bootret = sim_cz_residuals(2000, 350)

# find correlations
tempmat = bootret %>%
  group_by(portid) %>%
  arrange(portid,date) %>%
  mutate(dateid = row_number()) %>%
  select(dateid, portid, resid) %>%
  pivot_wider(names_from = portid, values_from = resid) %>%
  select(-dateid)

bootcor = cor(tempmat, use = 'pairwise.complete.obs')
bootcor2 = bootcor[lower.tri(bootcor)]

# read empirical data
ret = fread('output/pubpanel.csv') %>% 
  filter(
    !is.na(ret), year(date) >= 1963 # focus on more observed years
  ) 


# find correlations
tempmat = ret %>%
  group_by(signalname) %>%
  arrange(signalname,date) %>%
  select(date, signalname, ret) %>%
  pivot_wider(names_from = signalname, values_from = ret) %>%
  select(-date)

empcor = cor(tempmat, use = 'pairwise.complete.obs')
empcor2 = empcor[lower.tri(empcor)]

# assemble data
cdat = data.frame(
  c = bootcor2, group = 'sim', color = niceblue
) %>% rbind(
  data.frame(
    c = empcor2 , group = 'emp', color = 'gray'
  ) 
)

# plot
edge = seq(-1,1,0.1)
plotme = cdat %>% 
  group_by(group) %>% 
  summarise(
    cmid = hist(c,edge)$mids
    , density = hist(c,edge)$density
  ) %>% 
  mutate(
    group = factor(
      group
      , levels = c('sim','emp')
      , labels = c('Cluster-Bootstrap','CZ Data')
    )
  )


ggplot(
  plotme, aes(x=cmid, y=density, group = group)
) +
  geom_line(
    aes(linetype = group, color = group), size = 2
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12)
    , axis.text = element_text(size = 10)      
    , legend.title = element_blank()
    , legend.text = element_text(size = 10)
    , legend.key.size = unit(0.1, 'cm')
    , legend.position = c(80,80)/100
    , legend.key.width = unit(1,'cm')    
    , legend.spacing.y = unit(0.000001, 'cm')
    , legend.background = element_rect(colour = 'black', fill = 'white')    
  ) +
  labs(
    x = 'Pairwise Correlation'
    , y = 'Density'
  ) +
  scale_color_manual(
    values=c(niceblue, 'gray')
  ) +
  scale_linetype_manual(values = c('solid','31'))

# save 
ggsave(
  filename = 'output/semipar-cor.pdf', width = 5, height = 4
)

