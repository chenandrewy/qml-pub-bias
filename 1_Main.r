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
pubcross = pubcross %>% filter(tabs > 1.96)



# SET BASELINE BY HAND ====

source('0_Environment.r')
# lognorm: (2, 2, NA) fits well
# foldedt: (0, 0.5/0.2, 5) fits well
# exp (2, NA, NA)

# model parameters
par.hand = data.frame(
  pif = 0.4
  , truefam  = 'exp' # 'norm' or 'lognorm' or 'foldedt'
  , truepar1 = 2
  , truepar2 = NA
  , truepar3 = NA
  , pubfam  = 'logistic' # 'piecelin' or 'logistic' or 'stair' or 'trunc'
  , tmid    = 2.25
  , tslope  = 10
  , tmin    = 1
)

# make model tabs | pub object
tabspub.g.o = par2observed(par.hand)

# plot
p1 = plot_emp_vs_mod(tabspub.g.o, pubcross$tabs)
tt = seq(1,3,0.001)
p3 = tibble(
  tt = tt, prob = prob_pub(tt,par.hand)
) %>% 
  ggplot(aes(x=tt,y=prob)) + geom_line() +theme_minimal()

grid.arrange(p1,p3, nrow = 1)


# ONE ESTIMATE ====
source('0_Environment.r')
# doesn't work that great, unsurprising since because of the theorem
# pif often stays close to the initial guess

# estimation settings
est.point = list(
  parguess = par.hand %>% mutate(pif = 0.3)
  , namevec = c('pif', 'truepar1') # names of parameters to estimate  
)

# estimate
tempopts = opts.base
tempopts$xtol_rel = 0.001
tempopts$ftol_rel = 0.01
est.point = one_estimate(est.point, pubcross$tabs, tempopts)

# plot
tabspub.o = par2observed(est.point$par)
p1 = plot_2hist(
  pubcross$tabs, r(tabspub.o)(1e3), seq(0,18,1), 'data', 'model'
)

mu = create_mutrue_object(est.point$par)
p2 = tibble(mu = r(mu)(1e3)) %>% 
  ggplot(
    aes(x=mu)
  ) + 
  geom_histogram(aes(y=stat(density)))

tt = seq(0,4,0.25)
p3 = tibble(
  tt = tt, prob = prob_pub(tt,est.point$par)
) %>% 
  ggplot(aes(x=tt,y=prob)) + geom_point()

grid.arrange(p1,p2,p3, nrow = 1)

# output to console
est.point

# GRID ESTIMATE ====

source('0_Environment.r')

tic = Sys.time()

# estimation settings
est.grid = list(
  parguess = par.hand
  , namevec = c('pif', 'truepar1')
  , pifgrid = c(0.001, seq(0.05,0.95,0.05))
)

tempopts = opts.base
tempopts$print_level = 0
tempopts$xtol_rel = 0.0001
tempopts$ftol_rel = 0.001

# estimate
est.grid =  pifgrid_estimate(est.grid, pubcross$tabs, tempopts)

# plot 
tabspub.o = par2observed(est.grid$est$par)
p1 = plot_2hist(
  pubcross$tabs, r(tabspub.o)(1e3), seq(0,18,1), 'data', 'model'
)

mu = create_mutrue_object(est.grid$est$par)
p2 = tibble(mu = r(mu)(1e3)) %>%
  ggplot(
    aes(x=mu)
  ) +
  geom_histogram(aes(y=stat(density)))

p3 = est.grid$grid %>% ggplot(aes(x=pif,y=loglike)) + geom_point()

tt = seq(0,4,0.25)
p4 = tibble(
  tt = tt, prob = prob_pub(tt,est.grid$est$par)
) %>% 
  ggplot(aes(x=tt,y=prob)) + geom_point()

grid.arrange(p1,p2,p3,p4, nrow = 2)

# output to console
est.grid

toc = Sys.time()

toc - tic

# test plots ====

# plot point estimate
tempmod = par2observed(est.grid$est$par)
plot_emp_vs_mod(tempmod, pubcross$tabs)

# plot pif = 0.30 estimate
temppar = est.grid$grid %>% filter(pif == 0.3)
plot_emp_vs_mod(par2observed(temppar), pubcross$tabs)


# BOOTSTRAP ====
nboot = 100
nsignal = nrow(pubcross)
set.seed(1217)

id_all = matrix(
  sample(1:nsignal, size = nsignal*nboot, replace = T), nrow = nboot 
)

# estimation settings: use point estimate as initial guess
est.boot = list(
  parguess = par.hand
  , namevec = c('pif', 'truepar1')  
  , pifgrid = seq(0.05,0.95,0.05)
)

tempopts = opts.base 
tempopts$print_level = 0
tempopts$xtol_rel = 0.0001
tempopts$ftol_rel = 0.001

for (booti in 1:nboot){
  
  # bootstrap t-stats
  temp_tabs = pubcross$tabs[id_all[booti, ]]
  
  temp_est = pifgrid_estimate(est.boot, temp_tabs, tempopts)
  temp_bootdat = temp_est$est$par %>% 
    mutate(
      loglike = temp_est$est$loglike
      , hstar   = find_hurdle(temp_est$est$par)
    )
  
  if (booti == 1){
    bootdat = temp_bootdat
  } else {
    bootdat = rbind(bootdat, temp_bootdat)
  }
  
  print(temp_bootdat)
  
} # end for booti



# find t-hurdles ====

hlist = pubcross$tabs %>% sort()
hlist = hlist[hlist<4]




