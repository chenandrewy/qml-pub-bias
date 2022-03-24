# 2022 03 24 testing hand fitting of data to understand effect of parameteric assumptions

rm(list = ls())

# ENVIRONMENT  ====

source('0_Environment.r')

# bootstrap timing for reference
avail_hours = 12
nboot_wanted = 1000
min_per_boot = 12*60/nboot_wanted
min_per_boot

# DOWNLOAD DATA, MAKE PUBCROSS ====

# download signal documentation 
target_dribble = pathRelease %>% drive_ls() %>% 
  filter(name=='SignalDocumentation.xlsx')

drive_download(target_dribble, path = 'temp/deleteme.xlsx', overwrite = T)

signaldoc = left_join(
  read_excel('temp/deleteme.xlsx',sheet = 'BasicInfo')  
  , read_excel('temp/deleteme.xlsx',sheet = 'AddInfo')
) %>% 
  mutate(
    signalname = Acronym
    , pubdate = as.Date(paste0(Year, '-12-31'))
    , sampend = as.Date(paste0(SampleEndYear, '-12-31'))
    , sampstart = as.Date(paste0(SampleStartYear, '-01-01'))
  ) %>% 
  arrange(signalname) %>% 
  select(signalname, pubdate, sampend, sampstart)

# download returns and add samptype
target_dribble = pathRelease %>% drive_ls() %>%
  filter(name=='Portfolios') %>% drive_ls() %>%
  filter(name=='Full Sets OP') %>% drive_ls() %>%
  filter(name=='PredictorLSretWide.csv')

drive_download(target_dribble[1,], path = 'temp/deleteme2.csv', overwrite = T)

ret = fread('temp/deleteme2.csv') %>% 
  pivot_longer(
    -date, names_to = 'signalname', values_to = 'ret'
  ) %>% 
  filter(!is.na(ret)) %>% 
  left_join(signaldoc, by = 'signalname') %>% 
  mutate(
    samptype = case_when(
      (date >= sampstart) & (date <= sampend) ~ 'in-samp'
      , (date > sampend) & (date <= pubdate) ~ 'out-of-samp'
      , (date > pubdate) ~ 'post-pub'
      , TRUE ~ NA_character_
    )
  ) %>% 
  select(-pubdate,-sampend,-sampstart)

# summarize into pubcross
pubcross = ret %>% 
  filter(samptype == 'in-samp') %>% 
  group_by(signalname) %>% 
  summarize(
    rbar = mean(ret), vol = sd(ret), ndate = dplyr::n()
  ) %>% 
  mutate(
    tstat = rbar/vol*sqrt(ndate), tabs = abs(tstat)
  )

# filter 
#   remove 24 tabs < 1.96.  Modeling these doesn't pass cost / benefit
emp = pubcross %>% filter(tabs > 1.96)

# HAND FIT: EXPONENTIAL ====

# exponential
# pif = 0.001; scale = 2 # fits well

pif = 0.65; scale = 2.237 # 
t.o = UnivarMixingDistribution(
  Norm(0,1), Exp(1/scale) + Norm(0,1), mixCoeff = c(pif,1-pif) 
)
mod = Truncate(abs(t.o), 2)

plot_emp_vs_mod(mod, emp$tabs)

# testing: hand fit stair w/ triple mix normal ====

# model parameters
pif = 0.5; pia = 0.75; amu = 2.0; asig = 0.001; bmu = 5; bsig = 0.001; tslope = 0.5

# latent t-stats
t.o = UnivarMixingDistribution(
  Norm(0,1), Norm(amu,asig) + Norm(0,1), Norm(bmu,bsig) + Norm(0,1)
  , mixCoeff = c(pif, (1-pif)*pia, (1-pif)*(1-pia))
)
tabs.o = abs(t.o)

# selection
prob_pub = function(tt){
  (tt > 1.96 & tt <= 2.60)*tslope + (tt > 2.6)   
}

# observable
denom = tslope*(p(tabs.o)(2.6)-p(tabs.o)(1.96)) + (1-p(tabs.o)(2.6))
f_obs = function(tt) {
  d(tabs.o)(tt)*prob_pub(tt)/denom
}

# plot
edge = c(seq(0,15,1), 40)
freq = numeric(length(edge)-1)
for (i in 2:length(edge)){
  a = edge[i-1]; b = edge[i]
  freq[i-1] = integrate(f_obs, a, b)$value
}

hemp = hist(emptabs, edge, plot = F)

plotme = tibble(
  group = 'emp', tabs = hemp$mids, freq = hemp$density
) %>% rbind(
  tibble(
    group = 'mod', tabs = hemp$mids, freq = freq    
  )
)

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
  theme(legend.position = c(8,8)/10) +
  xlim(0,15)  



# testing: qml fit stair w/ simplest triple mix normal  ====

# model parameters
pif = 0.5; pia = 0.75; mua = 2; mub = 4
parvecguess = c(pif,pia,mua,mub)
parveclb =  c(0.01,0.01,02.0,0.01)
parvecub = c(0.99,0.99,10.0,10.0)

# selection (assumed)
prob_pub = function(tt){
  (tt > 1.96 & tt <= 2.60)*0.5 + (tt > 2.6)   
}

# likelihood 
negloglike = function(parvec){
  pif = parvec[1]
  pia = parvec[2]
  mua = parvec[3]
  mub = parvec[4]
  
  # latent t-stats
  t.o = UnivarMixingDistribution(
    Norm(0,1), Norm(mua,1), Norm(mub,1)
    , mixCoeff = c(pif, (1-pif)*pia, (1-pif)*(1-pia))
  )
  tabs.o = abs(t.o)
  
  # observable
  denom = tslope*(p(tabs.o)(2.6)-p(tabs.o)(1.96)) + (1-p(tabs.o)(2.6))
  
  # take on data
  tt = pubcross$tabs 
  
  onelikes = d(tabs.o)(tt)*prob_pub(tt)/denom
  return = -1*mean(log(onelikes))
  
} # end f_obs

negloglike(parvecguess) %>% print


# scope  out loglike
pigrid = seq(0.1,0.9,0.1)
for (i in 1:length(pigrid)){
  parvec = c(
    pigrid[i], pia, mua, mub
  )
  negloglike(parvec) %>% print  
}


opt = nloptr(
  x0 = parvecguess
  , lb = parveclb
  , ub = parvecub
  , eval_f = negloglike
  , opts = list(
    algorithm = 'NLOPT_LN_BOBYQA'
    , xtol_rel = 0.001
    , print_level = 3
  )
)
parvechat = opt$solution

pif = parvechat[1]
pia = parvechat[2]
mua = parvechat[3]
mub = parvechat[4]


# latent t-stats
t.o = UnivarMixingDistribution(
  Norm(0,1), Norm(mua,1), Norm(mub,1)
  , mixCoeff = c(pif, (1-pif)*pia, (1-pif)*(1-pia))
)
tabs.o = abs(t.o)

# observable
denom = tslope*(p(tabs.o)(2.6)-p(tabs.o)(1.96)) + (1-p(tabs.o)(2.6))
f_obs = function(tt) d(tabs.o)(tt)*prob_pub(tt)/denom

# plot
edge = c(seq(0,15,1), 40)
freq = numeric(length(edge)-1)
for (i in 2:length(edge)){
  a = edge[i-1]; b = edge[i]
  freq[i-1] = integrate(f_obs, a, b)$value
}

hemp = hist(emptabs, edge, plot = F)

plotme = tibble(
  group = 'emp', tabs = hemp$mids, freq = hemp$density
) %>% rbind(
  tibble(
    group = 'mod', tabs = hemp$mids, freq = freq    
  )
)

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
  theme(legend.position = c(8,8)/10) +
  xlim(0,15)  


# testing: flexible qml fit stair w/ simplest triple mix normal  ====

# model parameters
par.guess = data.frame(
  pif = 0.5, pia = 0.75, mua = 2, mub = 6
)
par.lb = data.frame(
  pif = 0.01, pia = 0.01, mua = 0.01, mub = 0.01
)
par.ub = data.frame(
  pif = 0.99, pia = 0.99, mua = 20, mub = 20
)

# function: translate par df to parvec vector form
par2parvec = function(par, namevec){
  parvec = par %>% select(all_of(namevec)) %>% as.numeric()
  return = parvec
} # end par2parvec

parvec2par = function(parvec,parbase,namevec){
  par = parbase
  par[ , namevec] = parvec
  return(par)
} # end parvec2par

# selection (assumed)
prob_pub = function(tt){
  (tt > 1.96 & tt <= 2.60)*0.5 + (tt > 2.6)   
}


# latent t-stats
make_t.o = function(par){
  t.o = UnivarMixingDistribution(
    Norm(0,1), Norm(par$mua,1), Norm(par$mub,1)
    , mixCoeff = c(par$pif, (1-par$pif)*par$pia, (1-par$pif)*(1-par$pia))
  )
  tabs.o = abs(t.o)
  return = tabs.o
}

# select names of pars to estimate
namevec = c('pif')

# likelihood 
negloglike = function(parvec){
  par = parvec2par(parvec,par.guess,namevec)
  
  # latent
  tabs.o = make_t.o(par)
  
  # observable
  denom = tslope*(p(tabs.o)(2.6)-p(tabs.o)(1.96)) + (1-p(tabs.o)(2.6))
  
  # take on data
  tt = pubcross$tabs 
  
  onelikes = d(tabs.o)(tt)*prob_pub(tt)/denom
  return = -1*mean(log(onelikes))
  
} # end f_obs

# test 1
parvec.guess = par2parvec(par.guess,namevec)
negloglike(parvec.guess) %>% print


# test 2
pigrid = seq(0.1,0.9,0.1)
for (i in 1:length(pigrid)){
  namevec = c('pif')
  parvec = pigrid[i]
  negloglike(parvec) %>% print  
}

# estimate! 
namevec = c('pif','pia')

opt = nloptr(
  x0 = par2parvec(par.guess, namevec)
  , lb = par2parvec(par.lb, namevec)
  , ub = par2parvec(par.ub, namevec)
  , eval_f = negloglike
  , opts = list(
    algorithm = 'NLOPT_LN_BOBYQA'
    , xtol_rel = 0.001
    , print_level = 3
  )
)
parvec.hat = opt$solution
par.hat = parvec2par(parvec.hat,par.guess,namevec)

## plot ====
par = par.hat

# latent 
tabs.o = make_t.o(par)

# observable
denom = tslope*(p(tabs.o)(2.6)-p(tabs.o)(1.96)) + (1-p(tabs.o)(2.6))
f_obs = function(tt) d(tabs.o)(tt)*prob_pub(tt)/denom


edge = c(seq(0,15,1), 40)
freq = numeric(length(edge)-1)
for (i in 2:length(edge)){
  a = edge[i-1]; b = edge[i]
  freq[i-1] = integrate(f_obs, a, b)$value
}

hemp = hist(emptabs, edge, plot = F)

plotme = tibble(
  group = 'emp', tabs = hemp$mids, freq = hemp$density
) %>% rbind(
  tibble(
    group = 'mod', tabs = hemp$mids, freq = freq    
  )
)

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
  theme(legend.position = c(8,8)/10) +
  xlim(0,15)  

par.hat

# testing: flexible qml fit stair w/ general triple mix normal  ====

# model parameters
par.guess = data.frame(
  pif = 0.5, pia = 0.75, mua = 2, siga = 0.001, mub = 6, sigb = 0.001
)
par.lb = data.frame(
  pif = 0.01, pia = 0.01, mua = 0.01, siga = 0.001, mub = 0.01, sigb = 0.001
)
par.ub = data.frame(
  pif = 0.99, pia = 0.99, mua = 20, siga = 20, mub = 20, sigb = 20
)

# function: translate par df to parvec vector form
par2parvec = function(par, namevec){
  parvec = par %>% select(all_of(namevec)) %>% as.numeric()
  return = parvec
} # end par2parvec

parvec2par = function(parvec,parbase,namevec){
  par = parbase
  par[ , namevec] = parvec
  return(par)
} # end parvec2par

# selection (assumed)
tslope = 0.5
prob_pub = function(tt){
  (tt > 1.96 & tt <= 2.60)*1 + (tt > 2.6)   
}

# latent t-stats
make_t.o = function(par){
  # t.o = UnivarMixingDistribution(
  #   Norm(0,1), Norm(par$mua,par$siga)+Norm(0,1), Norm(par$mub,par$sigb) + Norm(0,1)
  #   , mixCoeff = c(par$pif, (1-par$pif)*par$pia, (1-par$pif)*(1-par$pia))
  # )
  t.o = UnivarMixingDistribution(
    Norm(0,1), Norm(par$mua,1), Norm(par$mub,1)
    , mixCoeff = c(par$pif, (1-par$pif)*par$pia, (1-par$pif)*(1-par$pia))
  )  
  tabs.o = abs(t.o)
  return = tabs.o
}

# 
# # likelihood 
# negloglike = function(parvec){
#   par = parvec2par(parvec,par.guess,namevec)
#   
#   # latent
#   tabs.o = make_t.o(par)
#   
#   # observable
#   denom = tslope*(p(tabs.o)(2.6)-p(tabs.o)(1.96)) + (1-p(tabs.o)(2.6))
#   
#   # take on data
#   tt = pubcross$tabs 
#   
#   onelikes = d(tabs.o)(tt)*prob_pub(tt)/denom
#   return = -1*mean(log(onelikes))
#   
# } # end f_obs


# # estimate! 
# namevec = c('pif','pia')
# 
# opt = nloptr(
#   x0 = par2parvec(par.guess, namevec)
#   , lb = par2parvec(par.lb, namevec)
#   , ub = par2parvec(par.ub, namevec)
#   , eval_f = negloglike
#   , opts = list(
#     algorithm = 'NLOPT_LN_BOBYQA'
#     , xtol_rel = 0.001
#     , print_level = 3
#   )
# )
# parvec.hat = opt$solution
# par.hat = parvec2par(parvec.hat,par.guess,namevec)

## plot result ====
par = par.hat
par$pif = 0.01
par$pia = 0.82

par = data.frame(
  pif = 0.01, pia = 0.822, mua = 2, siga = 0.001, mub = 6, sigb = 0.001
)

# latent 
# tabs.o = make_t.o(par)
t.o = UnivarMixingDistribution(
  Norm(0,1), Norm(par$mua,par$siga)+Norm(0,1), Norm(par$mub,par$sigb) + Norm(0,1)
  , mixCoeff = c(par$pif, (1-par$pif)*par$pia, (1-par$pif)*(1-par$pia))
)
tabs.o = abs(t.o)

# observable
denom = tslope*(p(tabs.o)(2.6)-p(tabs.o)(1.96)) + (1-p(tabs.o)(2.6))
f_obs = function(tt) d(tabs.o)(tt)*prob_pub(tt)/denom

# setup
edge = c(seq(0,15,1), 40)
freq = numeric(length(edge)-1)
for (i in 2:length(edge)){
  a = edge[i-1]; b = edge[i]
  freq[i-1] = integrate(f_obs, a, b)$value
}

hemp = hist(emptabs, edge, plot = F)

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

# sanity check ====

## plot nice one ====
par = data.frame(
  pif = 0.01, pia = 0.822, mua = 2, mub = 6
)

# latent t-stats
make_t.o = function(par){
  t.o = UnivarMixingDistribution(
    Norm(0,1), Norm(par$mua,1), Norm(par$mub,1)
    , mixCoeff = c(par$pif, (1-par$pif)*par$pia, (1-par$pif)*(1-par$pia))
  )
  tabs.o = abs(t.o)
  return = tabs.o
}

# latent 
tabs.o = make_t.o(par)

# observable
denom = tslope*(p(tabs.o)(2.6)-p(tabs.o)(1.96)) + (1-p(tabs.o)(2.6))
f_obs = function(tt) d(tabs.o)(tt)*prob_pub(tt)/denom


edge = c(seq(0,15,1), 40)
freq = numeric(length(edge)-1)
for (i in 2:length(edge)){
  a = edge[i-1]; b = edge[i]
  freq[i-1] = integrate(f_obs, a, b)$value
}

hemp = hist(emptabs, edge, plot = F)

plotme = tibble(
  group = 'emp', tabs = hemp$mids, freq = hemp$density
) %>% rbind(
  tibble(
    group = 'mod', tabs = hemp$mids, freq = freq    
  )
)

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
  theme(legend.position = c(8,8)/10) +
  xlim(0,15)  

par.hat

## plot ugly one ? ====
par = data.frame(
  pif = 0.01, pia = 0.822, mua = 2, siga = 0.001, mub = 6, sigb = 0.001
)


# latent 
# tabs.o = make_t.o(par)
t.o = UnivarMixingDistribution(
  Norm(0,1), Norm(par$mua,par$siga)+Norm(0,1), Norm(par$mub,par$sigb) + Norm(0,1)
  , mixCoeff = c(par$pif, (1-par$pif)*par$pia, (1-par$pif)*(1-par$pia))
)
tabs.o = abs(t.o)

# observable
denom = tslope*(p(tabs.o)(2.6)-p(tabs.o)(1.96)) + (1-p(tabs.o)(2.6))
f_obs = function(tt) d(tabs.o)(tt)*prob_pub(tt)/denom

# setup
edge = c(seq(0,15,1), 40)
freq = numeric(length(edge)-1)
for (i in 2:length(edge)){
  a = edge[i-1]; b = edge[i]
  freq[i-1] = integrate(f_obs, a, b)$value
}

hemp = hist(emptabs, edge, plot = F)

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