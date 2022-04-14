# 2022 03 makes exhibits from bootall, the output of 2_CustomStats
rm(list = ls())
# SETUP ====

source('0_Settings.r')

library(latex2exp) # MAybe add to settings file?

# load bootstrap
load('intermediate/boot a-priori-model simple started 2022-04-11 13.Rdata')

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

# PLOTS: T-HURDLES ====

## Model fit ====

# the first boot is the point estimate
est.point = bootall.wide[1,] %>% 
  mutate(
    pubfam = set.boot$model_fam$pubfam[1]
    , mufam = set.boot$model_fam$mufam[1]
  )

# call function for plotting point estimate (in 0_Environment)
hist_emp_vs_par(cz_filt_tabs,est.point)

## Prop 1: Panel A ====
bootall.long %>%
  filter(stat %in% c('pif','pr_tgt_2')) %>% 
  ggplot(
    aes(x = value, fill = stat)
  ) +
  geom_histogram(
    aes(y = stat(count) / sum(count)),
    color = "black",
    position = 'dodge',
    ) +
  labs(
    y = "Frequency",
    x = "Estimate"
  ) +
  theme_classic() +
  theme(
    # Font sizes
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.text = element_text(size = 22),
    
    # Tweaking legend
    legend.position = c(0.75, 0.85),
    legend.text.align = 0,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = 5, r = 20, b = 5, l = 5), 
    legend.key.size = unit(2, "cm")
    ) +
  scale_fill_manual(
    values = c("pif" = "firebrick4", "pr_tgt_2" = "dodgerblue2"),
    name = NULL,
    breaks = c("pif", "pr_tgt_2"),
    labels = c(
      TeX(r"($\hat{\pi}_F$)"),
      TeX(r"($\hat{\textit{Pr}}\left( | \ t \ | > 1.96 \right)$)")
      )
    ) +
  coord_cartesian(xlim = c(0, 0.9), ylim = c(0.01, .20))

ggsave('output/prop_illustration1.pdf', width = 4, height = 4, scale = 2)

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
  geom_vline(xintercept = 0) +
  geom_histogram(
    aes(y = stat(count) / sum(count)),
    color = "black",
    fill = "dodgerblue2",
    bins = 30
    ) +
  annotate(geom = 'text', x = 0.5, y = 0.125, label = TeX(r"(Raise hurdle \rightarrow)"), size = 10) + 
  annotate(geom = 'text', x = -0.425, y = 0.125, label = TeX(r"(\leftarrow Lower hurdle)"), size = 10) +
  labs(
    y = "Frequency",
    x = TeX(r"($\hat{\pi}_F - \hat{\textit{Pr}}\left(|\ t\ | > 1.96\right)$)")
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),

    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    ) +
  coord_cartesian(ylim = c(0.01, .15))

ggsave('output/prop_illustration2.pdf', width = 4, height = 4, scale = 2)

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
    aes(x = h_fdr5, fill = pif_cat)
  ) +
  geom_vline(
    xintercept = qnorm(1-0.05/2),
    linetype = "dashed"
  ) +
  annotate('text', x = 2.8, y = 0.425, label = 'Classical hurdle', size = 10) +
  geom_histogram(
    aes(y = stat(count) / sum(count)),
    color = "black",
    breaks = breaks
  ) +
  labs(
    y = "Frequency",
    x = "t-hurdle for FDR = 5%"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.text = element_text(size = 22),
    
    legend.position = c(0.8, 0.525),
    legend.text.align = 0,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = 5, r = 15, b = 5, l = 5),
    legend.key.size = unit(1, "cm"),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
   scale_fill_manual(
     values = c("0" = "firebrick4", "1" = "dodgerblue2", "2" = "gold1"),
     name = NULL,
     breaks = c("0", "1", "2"),
     labels = c(
       TeX(r"($\hat{\pi}_F < \ \frac{1}{3}$)"),
       TeX(r"($\hat{\pi}_F \in \left( \ \frac{1}{3} \ , \ \frac{2}{3} \right]$)"),
       TeX(r"($\hat{\pi}_F > \ \frac{1}{3}$)")
     )
   ) +
  coord_cartesian(ylim = c(0.0215, .45))

ggsave('output/t-fdr05.pdf', width = 4, height = 3, scale = 2)

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
  geom_vline(
    xintercept = qnorm(1-0.01/2),
    linetype = "dashed"
  ) + 
  geom_histogram(
    aes(y = stat(count) / sum(count)),
    color = "black",
    breaks = breaks
  ) +
  labs(
    y = "Frequency",
    x = "t-hurdle for FDR = 1%"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  annotate('text', x = 1.60, y = 0.3, label = 'Classical hurdle', size = 10) +
  scale_fill_manual(values = c("0" = "firebrick4", "1" = "dodgerblue2", "2" = "gold1")) +
  coord_cartesian(ylim = c(0.0175, .375))
  
ggsave('output/t-fdr01.pdf', width = 4, height = 3, scale = 2)

# PLOTS: STRONGLY IDENTIFIED ====

## mean shrinkage ====

pif_cut = c(1/3,2/3)
breaks = seq(0,100,1)

bootall.long %>% 
  filter(stat %in% c('pif','bias_mean')) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(
    pif_cat = findInterval(pif, pif_cut)
    , pif_cat = as.factor(pif_cat)
    , bias_mean = bias_mean*100
  )  %>% 
  ggplot(
    aes(x=bias_mean, fill=pif_cat)
  ) +
  geom_histogram(
    aes(y = stat(count) / sum(count)),
    color = "black",
    breaks = breaks
  ) +
  geom_vline(
    xintercept = 26
  ) + 
  labs(
    y = "Frequency",
    x = "Mean bias (%)"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.text = element_text(size = 22),
    
    legend.position = c(0.8, 0.525),
    legend.text.align = 0,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = 5, r = 8.5, b = 5, l = 5),
    legend.key.size = unit(1, "cm"),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  scale_fill_manual(
    values = c("0" = "firebrick4", "1" = "dodgerblue2", "2" = "gold1"),
    name = NULL,
    breaks = c("0", "1", "2"),
    labels = c(
      TeX(r"($\hat{\pi}_F < \ \frac{1}{3}$)"),
      TeX(r"($\hat{\pi}_F \in \left( \ \frac{1}{3} \ , \ \frac{2}{3} \right]$)"),
      TeX(r"($\hat{\pi}_F > \ \frac{1}{3}$)")
    )
  ) +
  annotate('text', x = 14.5, y = 0.16, label = 'McLean-Pontiff bound', size = 10) +
  coord_cartesian(xlim = c(5, 40), ylim = c(0.01, 0.18))

ggsave('output/bias-mean.pdf', width = 4, height = 3, scale = 2)

## median shrinkage ====

pif_cut = c(1/3,2/3)
breaks = seq(0,100,1)

bootall.long %>% 
  filter(stat %in% c('pif','bias_med')) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(
    pif_cat = findInterval(pif, pif_cut)
    , pif_cat = as.factor(pif_cat)
    , bias_med = bias_med*100
  )  %>% 
  ggplot(
    aes(x=bias_med, fill=pif_cat)
  ) +
  geom_histogram(
    aes(y = stat(count) / sum(count)),
    color = "black",
    breaks = breaks
  ) +
  geom_vline(
    xintercept = 26
  ) +
  labs(
    y = "Frequency",
    x = "Median bias (%)"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.text = element_text(size = 22),
    
    legend.position = c(0.8, 0.525),
    legend.text.align = 0,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = 5, r = 8.5, b = 5, l = 5),
    legend.key.size = unit(1, "cm"),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  scale_fill_manual(
    values = c("0" = "firebrick4", "1" = "dodgerblue2", "2" = "gold1"),
    name = NULL,
    breaks = c("0", "1", "2"),
    labels = c(
      TeX(r"($\hat{\pi}_F < \ \frac{1}{3}$)"),
      TeX(r"($\hat{\pi}_F \in \left( \ \frac{1}{3} \ , \ \frac{2}{3} \right]$)"),
      TeX(r"($\hat{\pi}_F > \ \frac{1}{3}$)")
    )
  ) +
  annotate('text', x = 14.5, y = 0.16, label = 'McLean-Pontiff bound', size = 10) +
  coord_cartesian(xlim = c(5, 40), ylim = c(0.01, 0.18))

ggsave('output/bias-median.pdf', width = 4, height = 3, scale = 2)

## fdr pub mean ====

pif_cut = c(1/3,2/3)
breaks = seq(0,100,2.5)

bootall.long %>% 
  filter(stat %in% c('pif','fdrloc_mean')) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(
    pif_cat = findInterval(pif, pif_cut)
    , pif_cat = as.factor(pif_cat)
    , fdrloc_mean = fdrloc_mean*100
  )  %>% 
  ggplot(
    aes(x=fdrloc_mean, fill=pif_cat)
  ) +
  geom_histogram(
    aes(y = stat(count) / sum(count)),
    color = "black",
    breaks = breaks
  ) +
  labs(
    y = "Frequency",
    x = "Mean FDR among Published (%)"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.text = element_text(size = 22),
    
    legend.position = c(0.8, 0.525),
    legend.text.align = 0,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = 5, r = 8.5, b = 5, l = 5),
    legend.key.size = unit(1, "cm"),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  scale_fill_manual(
    values = c("0" = "firebrick4", "1" = "dodgerblue2", "2" = "gold1"),
    name = NULL,
    breaks = c("0", "1", "2"),
    labels = c(
      TeX(r"($\hat{\pi}_F < \ \frac{1}{3}$)"),
      TeX(r"($\hat{\pi}_F \in \left( \ \frac{1}{3} \ , \ \frac{2}{3} \right]$)"),
      TeX(r"($\hat{\pi}_F > \ \frac{1}{3}$)")
    )
  ) +
  coord_cartesian(xlim = c(0, 40))

ggsave('output/fdr-pub-mean.pdf', width = 4, height = 3, scale = 2)

## fdr pub median ====

pif_cut = c(1/3,2/3)
breaks = seq(0,100,2.5)

bootall.long %>% 
  filter(stat %in% c('pif','fdrloc_mean')) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(
    pif_cat = findInterval(pif, pif_cut)
    , pif_cat = as.factor(pif_cat)
    , fdrloc_mean = fdrloc_mean*100
  )  %>% 
  ggplot(
    aes(x=fdrloc_mean, fill=pif_cat)
  ) +
  geom_histogram(
    aes(y = stat(count) / sum(count)),
    color = "black",
    breaks = breaks
  ) +
  labs(
    y = "Frequency",
    x = "Mean FDR among Published (%)"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.text = element_text(size = 22),
    
    legend.position = c(0.8, 0.525),
    legend.text.align = 0,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = 5, r = 8.5, b = 5, l = 5),
    legend.key.size = unit(1, "cm"),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  scale_fill_manual(
    values = c("0" = "firebrick4", "1" = "dodgerblue2", "2" = "gold1"),
    name = NULL,
    breaks = c("0", "1", "2"),
    labels = c(
      TeX(r"($\hat{\pi}_F < \ \frac{1}{3}$)"),
      TeX(r"($\hat{\pi}_F \in \left( \ \frac{1}{3} \ , \ \frac{2}{3} \right]$)"),
      TeX(r"($\hat{\pi}_F > \ \frac{1}{3}$)")
    )
  ) +
  coord_cartesian(xlim = c(0, 40))

ggsave('output/fdr-pub-median.pdf', width = 4, height = 3, scale = 2)
