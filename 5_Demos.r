# 2022 04 porting demo figs (intro and hlz) from Matlab back to R

# ENVIRONMENT ====

rm(list=ls())


source('0_Settings.r')

nsim = 200
nfac = 1e3

makehist = function(dat, edge = seq(0,10, 0.5)){
  tmid = edge[1:(length(edge)-1)] + 0.5*(edge[2]-edge[1])
  hdat2 = data.frame(
    bin = 1:length(tmid), tmid
  )
  
  hdat = dat %>% 
    filter(t>min(edge), t<max(edge)) %>% 
    group_by(group) %>% 
    mutate(
      bin = findInterval(t,edge)
    ) %>% 
    group_by(group,bin) %>% 
    summarize(
      n = dplyr::n()/nsim
    ) %>% 
    left_join(hdat2)  
} # end function makehist

# INTRO FIG ====

# translated from main_demo_intro.m and main-intuition.R on 2022 04

## shared setup ====

# parameters hand-selected so that pr_tgt2 = 0.5 and pr_tgt3_tgt2 = 0.5
# too lazy too code up numerical solver

simall = tibble(
  t = abs(rnorm(nfac*nsim,0,1)), group = F
) %>% 
  rbind(
    tibble(
      t  = rgamma(n= nfac*nsim, shape = 2, scale = 0.8) + 1.65, group = T
    )
  ) %>% 
  mutate(
    group = factor(
      group, levels = c(T,F)
    )
  )


## FIG: REALISTIC DEMO ====

dat = simall

fdr = dat %>% filter(t>1.96) %>% 
  summarize(fdr = round(mean(group == F)*100))

hdat= makehist(dat)

ggplot(hdat, aes(x=tmid, y=n, fill=group)) + 
  geom_bar(stat='identity', position='stack') +
  geom_vline(xintercept = 1.96, size = 1) + 
  labs(title = "", x = "|t-statistic|", y = "Number of Factors") +
  scale_fill_manual(
    labels = c("True Factor", "False Factor"), 
    values = c(MATBLUE, MATRED),
    name = ""
  ) +
  chen_theme +
  annotate(
    geom = "text",
    label = "t-hurdle = 1.96",
    x = 4, y = 375, size = 6
  ) +
  geom_segment(aes(x = 5.7, y = 95, xend = 2.3, yend = 15),
             arrow = arrow(length = unit(0.03, "npc")),
             colour = "black", size = 0.3
  ) +
  annotate(
    geom = "text",
    label = paste0("FDR = ", fdr$fdr, '%'),
    x = 7, y = 100
    , size = 6
  )  +
  scale_x_continuous(breaks = seq(0,10,2)) +
  theme(
    legend.position = c(0.75, 0.75)
    , legend.key.size = unit(0.6,'cm')
    , legend.margin = margin(t = -10, r = 5, b = 5, l = 5)
    , legend.text = element_text(size = 20)    
  ) 
  



ggsave('output/intro-realistic.pdf', width = 4, height = 4, scale = 1.5, device = cairo_pdf)



## FIG: WRONG INTUITION DEMO ====

dat = simall %>% filter(group == F)

fdr = dat %>% filter(t>1.96) %>% 
  summarize(fdr = round(mean(group == F)*100))

hdat= makehist(dat)

plot = ggplot(hdat, aes(x=tmid, y=n, fill=group)) + 
  geom_bar(stat='identity', position='stack') +
  geom_vline(xintercept = 1.96, size = 1) + 
  labs(title = "", x = "t-statistic", y = "Factors") +
  scale_fill_manual(
    labels = c("False Factor"), 
    values = c(MATRED),
    name = ""
  ) +
  chen_theme +
  annotate(
    geom = "text",
    label = "t-hurdle = 1.96",
    x = 2.9,
    y = 375
    , size = 6
  ) +
  geom_segment(aes(x = 2.7, y = 70, xend = 2.3, yend = 20),
               arrow = arrow(length = unit(0.03, "npc")),
               colour = "black", size = 0.3
  ) +
  annotate(
    geom = "text",
    label = paste0("FDR = ", fdr$fdr, '%'),
    x = 3.5, y = 80, size = 6
  ) +
  theme(
    legend.position = c(0.75, 0.75)
    , legend.key.size = unit(0.6,'cm')
    , legend.margin = margin(t = -10, r = 5, b = 5, l = 5)
    , legend.text = element_text(size = 20)    
  )   

plot

ggsave('output/intro-wrong.pdf', width = 4, height = 4, scale = 1.5, device = cairo_pdf)



# HLZ FIG ====
nsim = 100
nfac = 1000

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


dat = rbind(
  tibble(mu = mumix, t = mumix + u, group = 'hlz')  
  , tibble(mu = mutrue, t = mutrue + u, group = 'alt')
) %>% 
mutate(
  group = factor(group, levels = c('hlz','alt'))
  , t = abs(t)
)



# find tail normalization
sumdat = dat %>% 
  group_by(group) %>% 
  filter(t>2.57) %>% 
  summarize(n=dplyr::n())

npub_alt_hlz = sumdat$n[2]/sumdat$n[1]

hdat = makehist(dat, edge = seq(-4,10,0.2))

hdat_norm = hdat %>% 
  left_join(
    hdat %>% group_by(group) %>% summarise(ntot = sum(n))
  ) %>% 
  mutate(
    n = n/ntot
    , n = if_else(group=='alt', n/npub_alt_hlz, n)
  )  
  

# plot
plot = ggplot(hdat_norm, aes(x=tmid, y=n, fill=group)) + 
  geom_bar(stat='identity', position='identity', alpha = 0.5) +
  labs(title = "", x = "|t-statistic|", y = "Density") +
  scale_fill_manual(
    labels = unname(TeX(c(
      "Model 1: HLZ's Baseline ($\\widehat{\\pi}_F=0.44)$"
      ,"Model 2: HLZ's Baseline w/ $\\pi_F=0$"
      ))), 
    values = c(MATBLUE, MATRED),
    name = ""
  )  +
  theme_light() +
  theme(
    legend.position = c(0.7, 0.65)
    , legend.title = NULL
    , legend.key.size = unit(0.5,'cm')  
    , legend.text.align = 0    
    , legend.text = element_text(size=9)
    , panel.grid.minor = element_blank()
    , panel.grid.major = element_blank()
    , panel.border = element_blank()
    , axis.line = element_line(colour = 'black')
    ,     text = element_text(family = "Palatino Linotype")
  ) +
  geom_vline(xintercept = 2.57, size = 1, color = 'red') +  
  annotate(
    geom = "text", label = "Observed ->",
    x = 3.6, y = 0.118, size = 3.2, color = 'red'
  ) +
  annotate(
    geom = "text",
    label = "<- Extrapolated",
    x = 1.35, y = 0.118, size = 3.2, color = 'red'
  ) +
  xlim(0, 8)+ ylim(0,0.12) 
  

plot

ggsave('output/hlz-demo.pdf', width = 3, height = 2, scale = 1.5, device = cairo_pdf)

