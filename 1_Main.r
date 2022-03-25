rm(list = ls())

# ENVIRONMENT  ====

source('0_Environment.r')

# bootstrap timing for reference
avail_hours = 12
nboot_wanted = 1000
min_per_boot = 12*60/nboot_wanted
min_per_boot

# IMPORT PUBCROSS ====

pubcross = import_cz(dl = F)

# ESTIMATE ====

# filter 
#   remove 24 tabs < 1.96.  Modeling these may not pass cost / benefit
samp = pubcross %>% filter(tabs > 1.96)

# model parameters
par.guess = data.frame(
  pif = 0.1
  , pia = 0.25
  , mua = 2
  , siga = 0.1
  , mub = 5
  , sigb = 2
  , pubfam = 'stair' # 'stair' or 'piecelin'
  , pubpar1 = NA
  , pubpar2 = 0.75
)
par.lb = data.frame(
  pif = 0.001
  , pia = 0.001
  , mua = 0.5
  , siga = 0.001
  , mub = 1
  , sigb = 0.0001
  , pubpar1 = NA
  , pubpar2 = 0.3
)
par.ub = data.frame(
  pif = 0.99
  , pia = 0.99
  , mua = 4
  , siga = 10
  , mub = 10
  , sigb = 10
  , pubpar1 = NA  
  , pubpar2 = 1
)

## check pub_prob ====
tt = seq(1,3,0.001)
tibble(
  tt = tt, prob = pub_prob(tt,par.guess)
) %>% 
  ggplot(aes(x=tt,y=prob)) + geom_line() +theme_minimal()


## ====


# parameters to estimate
namevec = c('pif', 'pia','mua','mub','siga','sigb','pubpar2')

# test objective
negloglike(par2parvec(par.guess,namevec))

# opts
opts = list(
  algorithm = 'NLOPT_GN_CRS2_LM' # 'NLOPT_LN_BOBYQA' or 'NLOPT_GN_DIRECT_L' or 'NLOPT_GN_CRS2_LM'
  , xtol_rel = 0.001
  , print_level = 3
  , maxeval = 1000
)

opt = nloptr(
  x0 = par2parvec(par.guess, namevec)
  , lb = par2parvec(par.lb, namevec)
  , ub = par2parvec(par.ub, namevec)
  , eval_f = negloglike
  , opts = opts
)
parvec.hat = opt$solution
par.hat = parvec2par(parvec.hat,par.guess,namevec)

## plot result ====
par = par.hat


# observable
f_obs = function(tt){
  make_single_likes(tt,par)
}

# setup
edge = c(seq(0,15,1), 40)
freq = numeric(length(edge)-1)
for (i in 2:length(edge)){
  a = edge[i-1]; b = edge[i]
  freq[i-1] = integrate(f_obs, a, b)$value
}

hemp = hist(pubcross$tabs, edge, plot = F)

plotme = tibble(
  group = 'emp', tabs = hemp$mids, freq = hemp$density
) %>% rbind(
  tibble(
    group = 'mod', tabs = hemp$mids, freq = freq    
  )
) %>% 
  filter(!is.na(freq))

stat = plotme %>% 
  pivot_wider(names_from = group, values_from = freq) %>% 
  mutate(
    error = 1000*(mod - emp)
  ) 

print(stat)
print(mean(abs(stat$error)))

plotme %>% 
  ggplot(aes(x=tabs,y=freq,group=group,color=group)) +
  geom_line() +
  geom_point(size=3) +
  theme_minimal() +
  scale_color_manual(values = c('blue','red')) + 
  theme(legend.position = c(0.8,0.8)) +
  xlim(0,15)  

par

opt$objective
