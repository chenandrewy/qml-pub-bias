# SETUP ====



# rm(list = ls())
source('0_Settings.r')





## Font setup ====
# only needs to be done once
# https://fromthebottomoftheheap.net/2013/09/09/preparing-figures-for-plos-one-with-r/
# https://stackoverflow.com/questions/28863412/error-using-arial-in-eps-figure-with-extrafont-package

loadfonts(device = 'win')
loadfonts(device = 'postscript')
loadfonts()


## Load data ====

# load semi-par bootstrap for robustness
load('intermediate/boot a-priori-model semipar started 2022-04-14 22.Rdata')

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
load('intermediate/boot a-priori-model simple started 2022-04-15 03.Rdata')

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

# PLOTS: MODEL FIT ====

## estimate alternative model ====

set.alt = set.boot
set.alt$model_fam$pif = 0.5
par.guess = par.guess.all[1,] %>% mutate(pif = 0.5)

temp = estimate(set.alt, cz_filt_tabs, par.guess, print_level = 0)

est.alt = temp


## setup ====
# make f for point estimate
par.point = bootall.wide[1,] 

tabs.o = make_tabs_object(par.point)
f_point = function(tt) d(tabs.o)(tt)

# make f_obs for the alternative
tabs.o.alt = make_tabs_object(est.alt$par)
f_alt = function(tt) d(tabs.o.alt)(tt)

# integrate to find hist frequencies
delta.m = 0.1
edge = c(seq(0,15,delta.m))
mids = edge[1:length(edge)-1] + delta.m/2
freq.point = numeric(length(edge)-1)
freq.alt = freq.point
for (i in 2:length(edge)){
  a = edge[i-1]; b = edge[i]
  freq.point[i-1] = integrate(f_point, a, b)$value
  freq.alt[i-1] = integrate(f_alt, a, b)$value 
}

# assemble into tibble
dat.point = tibble(mids = mids, freq = freq.point, group = 'point')
dat.alt = tibble(mids = mids, freq = freq.alt, group = 'alt')

# make empirical data
temptabs = import_cz(dl = F) %>% pull(tabs)

delta.emp = 0.5
edge.emp = seq(0,15,delta.emp)
hemp = hist(temptabs, edge.emp, plot = F)
dat.emp = tibble(mids = hemp$mids, freq = hemp$count, group = 'emp')

# rescale
dat.all = rbind(dat.point,dat.alt) %>% rbind(dat.emp)

scale.df = dat.all %>% filter(mids >= 3) %>% 
  group_by(group) %>%  
  summarize(tot = sum(freq)) %>% 
  pivot_wider(names_from = group, values_from = tot)

dat.all = dat.all %>% 
  mutate(
    freq = if_else(group == 'point', freq*scale.df$emp/scale.df$point * delta.emp/ delta.m, freq)
    , freq = if_else(group == 'alt', freq*scale.df$emp/scale.df$alt * delta.emp/ delta.m, freq)
  )


## plot ====

dat.all %>% filter(group == 'emp') %>% 
  ggplot(
    aes(x=mids, y=freq)
  ) +
  geom_bar(aes(fill = 'emp'), stat = 'identity') +
  geom_line(
    data = dat.all %>% filter(group == 'alt')
    , aes(color = 'alt')
    , linetype = 'dashed'
    , size = 1
    ) +  
  geom_line(
    data = dat.all %>% filter(group == 'point')
    , aes(color = 'point')
    , size = 1
    ) +
  scale_fill_manual(
    values = c('emp' = 'gray'), name = NULL
    , labels = c('Chen-Zimmermann Data')
  ) +
  scale_color_manual(
    breaks = c('point','alt')
    , values = c(
      'point' = MATBLUE
      , 'alt' = MATRED)
    , labels = c(
      TeX('Point Estimate ($\\widehat{\\pi}_F = 0.01$)')
      , TeX('Estimate Fixing $\\pi_F = $0.50')
    )
  ) +
  chen_theme +
  theme(
    legend.position = c(0.70, 0.65)
    , legend.margin = margin(t = 0, r = 5, b = 5, l = 5)
    , legend.key.size = unit(1, 'cm')
    , legend.text = element_text(size = 20)
    , panel.border = element_rect(colour = "black", fill=NA, size=1)
    , axis.ticks.length=unit(0.2, "cm")
  ) + 
  xlab('|t-statistic|') + ylab('Number of Predictors') +
  coord_cartesian(xlim = c(0,15), ylim = c(0, 80))  +
  scale_x_continuous(breaks = seq(0,15,2.5))

ggsave(filename = 'output/model-fit.pdf', width = 4, height = 3, scale = 2.5, device = cairo_pdf)

# PLOTS: T-HURDLES ====


## Prop 1: Panel A ====
breaks = seq(0,1.0,0.05)
p = bootall.long %>%
  filter(stat %in% c('pif','pr_tgt_2')) %>% 
  ggplot(
    aes(x = value, fill = stat)
  ) +
  geom_histogram(
    aes(y = stat(count) / sum(count)),
    color = "black",
    position = 'identity',
    breaks = breaks
    , alpha = 0.6
  ) +
  labs(
    y = "Frequency",
    x = "Estimate"
  ) +
  chen_theme +
  scale_fill_manual(
    values = c("pif" = MATRED, "pr_tgt_2" = MATBLUE),
    name = NULL,
    breaks = c("pif", "pr_tgt_2"),
    labels = unname(c(
      TeX('$\\hat{\\pi}_F$'), 
      TeX('$\\widehat{Pr}(|t_i|>1.96)')
    ))
  ) +
  coord_cartesian(xlim = c(0, 0.9), ylim = c(0.01, .20))

p

# the device = cairo_pdf thing is super necessary
ggsave(plot = p, filename = 'output/prop_illustration1.pdf', width = 4, height = 4, scale = 1.5, device = cairo_pdf)

## Prop 1: Panel B ====


breaks = seq(-0.7,0.7,0.05)
p = bootall.long %>% 
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
    fill = 'gray',
    breaks = breaks
  ) +
  annotate(geom = 'text', x = -0.37, y = 0.132, label = TeX(r"($\leftarrow$ Lower hurdle)"), size = 8) +  
  annotate(geom = 'text', x = 0.36, y = 0.132, label = TeX(r"(Raise hurdle $\rightarrow$)"), size = 8) + 
  labs(
    y = "Frequency",
    x = TeX("$\\widehat{\\pi}_F - \\widehat{Pr}(|t_i|>1.96)$")
  ) +
  chen_theme +
  coord_cartesian(ylim = c(0.01, .15))
p 

ggsave(plot = p, filename = 'output/prop_illustration2.pdf', width = 4, height = 4, scale = 1.5, device = cairo_pdf)

## t-hurdle: fdr = 5% ====

pif_cut = c(1/3,2/3)
edge = seq(0,4.0,0.5)

ntot = sum(!is.na(bootall.long %>% filter(stat=='pif') %>% pull(value)))

plotme = bootall.long %>% 
  filter(stat %in% c('pif','h_fdr5')) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(
    pif_cat = findInterval(pif, pif_cut)
    , pif_cat = as.factor(pif_cat)
  ) %>% 
  group_by(pif_cat) %>% 
  summarize(
    counts = hist(h_fdr5, edge)$counts
    , mids = hist(h_fdr5, edge)$mids
    , freq = counts/ntot
  )   

# make plot
plotme %>% 
  ggplot(
    aes(x=mids, y=freq, group = pif_cat)
  ) +
  geom_vline(
    xintercept = qnorm(1-0.05/2),
    linetype = "dashed"
  ) +  
  geom_bar(
    aes(fill = pif_cat), stat = 'identity', width = 0.35
  )  +
  annotate('text', x = 2.4, y = 0.4, label = 'classical hurdle', size = 7) +
  labs(
    y = "Frequency",
    x = "t-hurdle"
  ) +
  chen_theme +
  scale_fill_manual(
    values = c("0" = MATBLUE, "1" = MATYELLOW, "2" = MATRED),
    name = NULL,
    breaks = c("0", "1", "2"),
    labels = unname(TeX(c(
      "$\\widehat{\\pi}_F \\leq 1/3$",
      "$\\widehat{\\pi}_F \\in  (1/3, 2/3\\]$",
      "$\\widehat{\\pi}_F  > 2/3$"
    )))
  ) +
  scale_x_continuous(breaks = seq(0,4.5, 0.5) ) +
  coord_cartesian(ylim = c(0.0215, .40), xlim = c(0,4.0)) +
  theme(
    legend.position = c(0.30, 0.65)
    , legend.margin = margin(t = 0, r = 5, b = 5, l = 5)
    , legend.key.size = unit(1, 'cm')
    , legend.text = element_text(size = 20)
    , panel.border = element_rect(colour = "black", fill=NA, size=1)
    , axis.ticks.length=unit(0.2, "cm")
  ) 


ggsave('output/t-fdr05.pdf', width = 6, height = 2.5, scale = 2.0, device = cairo_pdf)  

## t-hurdle: fdr = 1% ====

pif_cut = c(1/3,2/3)
edge = seq(0,5.0,0.5)

ntot = sum(!is.na(bootall.long %>% filter(stat=='pif') %>% pull(value)))

plotme = bootall.long %>% 
  filter(stat %in% c('pif','h_fdr1')) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(
    pif_cat = findInterval(pif, pif_cut)
    , pif_cat = as.factor(pif_cat)
  ) %>% 
  group_by(pif_cat) %>% 
  summarize(
    counts = hist(h_fdr1, edge)$counts
    , mids = hist(h_fdr1, edge)$mids
    , freq = counts/ntot
  )   

# make plot
plotme %>% 
  ggplot(
    aes(x=mids, y=freq, group = pif_cat)
  ) +
  geom_vline(
    xintercept = qnorm(1-0.01/2),
    linetype = "dashed"
  ) +  
  geom_bar(
    aes(fill = pif_cat), stat = 'identity', width = 0.35
  )  +
  annotate('text', x = 3.0, y = 0.4, label = 'classical hurdle', size = 7) +
  labs(
    y = "Frequency",
    x = "t-hurdle"
  ) +
  chen_theme +
  scale_fill_manual(
    values = c("0" = MATBLUE, "1" = MATYELLOW, "2" = MATRED),
    name = NULL,
    breaks = c("0", "1", "2"),
    labels = unname(TeX(c(
      "$\\widehat{\\pi}_F \\leq 1/3$",
      "$\\widehat{\\pi}_F \\in  (1/3, 2/3\\]$",
      "$\\widehat{\\pi}_F  > 2/3$"
    )))
  )   +
  scale_x_continuous(breaks = seq(0,4.0, 0.5) ) +
  coord_cartesian(ylim = c(0.0215, .40), xlim = c(0,4.0)) +
  theme(
    legend.position = c(0.30, 0.65)
    , legend.margin = margin(t = 0, r = 5, b = 5, l = 5)
    , legend.key.size = unit(1, 'cm')
    , legend.text = element_text(size = 20)
    , panel.border = element_rect(colour = "black", fill=NA, size=1)
    , axis.ticks.length=unit(0.2, "cm")
  ) 

ggsave('output/t-fdr01.pdf', width = 6, height = 2.5, scale = 2.0, device = cairo_pdf)  
   


# PLOTS: STRONGLY IDENTIFIED ====

## mean shrinkage ====

pif_cut = c(1/3,2/3)
edge = seq(0,50,5)

plotme0 = bootall.wide %>% 
  transmute(pif, x = bias_mean) %>% 
  mutate(x = 100*x) %>%  
  filter(!is.na(x), x>-20) # for some reason there's 1/1000 obs is negative

ntot = dim(plotme0)[1]

plotme = plotme0 %>% 
  mutate(
    pif_cat = findInterval(pif, pif_cut)
    , pif_cat = as.factor(pif_cat)
  ) %>% 
  group_by(pif_cat) %>% 
  summarize(
    counts = hist(x, edge)$counts
    , mids = hist(x, edge)$mids
    , freq = counts/ntot
  )   

# make plot
plotme %>% 
  ggplot(
    aes(x=mids, y=freq, group = pif_cat)
  ) +
  geom_vline(
    xintercept = 26,
    linetype = "dashed"
  ) +  
  annotate('text', x = 33.5, y = 0.39, label = 'McLean-Pontiff Upper Bound', size = 7) +  
  geom_bar(
    aes(fill = pif_cat), stat = 'identity', width = 4
  )  +
  labs(
    y = "Frequency"
    , x = TeX(r"(Mean $\log \,\left[ t_i / \[ideal \, t_i\] \right]$ for Published $i$ (\%) )")
  ) +
  chen_theme +
  scale_fill_manual(
    values = c("0" = MATBLUE, "1" = MATYELLOW, "2" = MATRED),
    name = NULL,
    breaks = c("0", "1", "2"),
    labels = unname(TeX(c(
      "$\\widehat{\\pi}_F \\leq 1/3$",
      "$\\widehat{\\pi}_F \\in  (1/3, 2/3\\]$",
      "$\\widehat{\\pi}_F  > 2/3$"
    )))
  )   +
  scale_x_continuous(breaks = seq(0,40, 5) ) +
  coord_cartesian(ylim = c(0.0215, .40), xlim = c(0,40)) +
  theme(
    legend.position = c(0.80, 0.65)
    , legend.margin = margin(t = 0, r = 5, b = 5, l = 5)
    , legend.key.size = unit(1, 'cm')
    , legend.text = element_text(size = 20)
    , panel.border = element_rect(colour = "black", fill=NA, size=1)
    , axis.ticks.length=unit(0.2, "cm")
  ) 

ggsave('output/bias-mean.pdf', width = 6, height = 2.5, scale = 2.0, device = cairo_pdf)  

## fdr pub ====


### HLZ's estimate ====
nsim = 1000
nfac = 2000

rho = 0.2

# use model u = c + v, Cov(c,v) = 0, so rho = sigv^2/(sigu*sigv) = sigv/sigu
# so sigv = sigu*rho = 1*rho 
# and then sigu^2 = sigc^2 + sigv^2 so sigc = sqrt(1-rho^2)

sigv = rho
sigc = sqrt(1-rho^2)

v = rnorm(nfac*nsim, 0, sigv)
c = rnorm(nfac*nsim, 0, sigc)
u = c+v

# alt model first
lambda_mu = 0.55/(15/sqrt(12)/sqrt(240)) # p 26, 29
mutrue = rexp(n = nfac*nsim, rate = 1/lambda_mu)

# hlz baseline
type = runif(nfac*nsim) > 0.444
mumix = mutrue
mumix[type == F] = 0

hlzdat = tibble(type = type, mu = mumix, t = mumix + u)  %>% 
  mutate(
    tabs = abs(t)
    , pub = (tabs > 2.57) | (tabs > 1.96 & tabs < 2.57 & runif(1) > 0.5 )
  )

hlz_estimate = hlzdat %>% 
  filter(pub) %>% 
  summarize(
    fdr = mean(!type)*100
  ) %>% 
  pull(fdr)

### Plot my stuff with HLZ====

pif_cut = c(1/3,2/3)
edge = seq(0,50,5)

plotme0 = bootall.wide %>% 
  transmute(pif, x = fdrloc_mean) %>% 
  mutate(x = 100*x) %>%  
  filter(!is.na(x))

ntot = dim(plotme0)[1]

plotme = plotme0 %>% 
  mutate(
    pif_cat = findInterval(pif, pif_cut)
    , pif_cat = as.factor(pif_cat)
  ) %>% 
  group_by(pif_cat) %>% 
  summarize(
    counts = hist(x, edge)$counts
    , mids = hist(x, edge)$mids
    , freq = counts/ntot
  )   

# make plot
plotme %>% 
  ggplot(
    aes(x=mids, y=freq, group = pif_cat)
  ) +
  geom_vline(
    xintercept = hlz_estimate
    , linetype = "dashed"
  ) +  
  annotate('text', x = 15.0, y = 0.39, label = 'HLZ Baseline Estimate', size = 7) +  
  geom_bar(
    aes(fill = pif_cat), stat = 'identity', width = 4
  )  +
  labs(
    y = "Frequency"
    , x = TeX(r"(Mean of $Pr( F_i |t_i)$ for Published $i$  (\%) )")
  ) +
  chen_theme +
  scale_fill_manual(
    values = c("0" = MATBLUE, "1" = MATYELLOW, "2" = MATRED),
    name = NULL,
    breaks = c("0", "1", "2"),
    labels = unname(TeX(c(
      "$\\widehat{\\pi}_F \\leq 1/3$",
      "$\\widehat{\\pi}_F \\in  (1/3, 2/3\\]$",
      "$\\widehat{\\pi}_F  > 2/3$"
    )))
  )   +
  scale_x_continuous(breaks = seq(0,40, 5) ) +
  coord_cartesian(ylim = c(0.0215, .40), xlim = c(0,40)) +
  theme(
    legend.position = c(0.80, 0.65)
    , legend.margin = margin(t = 0, r = 5, b = 5, l = 5)
    , legend.key.size = unit(1, 'cm')
    , legend.text = element_text(size = 20)
    , panel.border = element_rect(colour = "black", fill=NA, size=1)
    , axis.ticks.length=unit(0.2, "cm")
  ) 

ggsave('output/fdr-pub.pdf', width = 6, height = 2.5, scale = 2.0, device = cairo_pdf)  




# EXTRA EXHIBITS ====

## semi-par boot correlations ====


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

