source('0_Environment.r')

# BOOTSTRAP LIMITS ====

avail_hours = 12
nboot_wanted = 1000

min_per_boot = 12*60/nboot_wanted

min_per_boot


# DOWNLOAD DATA, MAKE PUBCROSS ====

# root of April 2021 release on Gdrive
pathRelease = 'https://drive.google.com/drive/folders/1I6nMmo8k_zGCcp9tUvmMedKTAkb9734R'

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


# FIND INITIAL GUESS BY HAND ====
# lognorm: (2, 2, NA) fits well
# foldedt: (0, 0.5/0.2, 5) fits well

# model parameters
par.hand = data.frame(
  pif = 0.4
  , truefam  = 'lognorm' # 'norm' or 'lognorm' or 'foldedt'
  , truepar1 = 2
  , truepar2 = 2
  , truepar3 = NA
  , pubfam  = 'logistic' # 'piecelin' or 'logistic
  , tmid    = 2 # prob of pub is 0.5 at the midpoint
  , tslope  = 4
)

# make model tabs | pub object
tabspub.g.o = par2observed(par.hand)

# plot
plot_2hist(
  pubcross$tabs, r(tabspub.g.o)(2e3), seq(0,15,0.5), 'emp', 'sim'
)



# ONE ESTIMATE ====
source('0_Environment.r')
# doesn't work that great, unsurprising since because of the theorem
# pif often stays close to the initial guess

# estimation settings
est.point = list(
  parguess = par.hand %>% mutate(pif = 0.6)
  , namevec = c('pif', 'truepar1','truepar2') # names of parameters to estimate  
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


grid.arrange(p1,p2, nrow = 1)

# output to console
est.point



# GRID ESTIMATE ====

source('0_Environment.r')

tic = Sys.time()

# estimation settings
est.point = list(
  parguess = par.hand
  , namevec = c('pif', 'truepar1','truepar2', 'tmid','tslope')  
  , pifgrid = seq(0.05,0.95,0.05)
)

tempopts = opts.base
tempopts$print_level = 0
tempopts$xtol_rel = 0.0001
tempopts$ftol_rel = 0.001

# estimate
est.point =  pifgrid_estimate(est.point, pubcross$tabs, tempopts)

# plot ====
tabspub.o = par2observed(est.point$est$par)
p1 = plot_2hist(
  pubcross$tabs, r(tabspub.o)(1e3), seq(0,18,1), 'data', 'model'
)

mu = create_mutrue_object(est.point$est$par)
p2 = tibble(mu = r(mu)(1e3)) %>%
  ggplot(
    aes(x=mu)
  ) +
  geom_histogram(aes(y=stat(density)))

p3 = est.point$grid %>% ggplot(aes(x=pif,y=loglike)) + geom_point()

tt = seq(0,4,0.25)
p4 = tibble(
  tt = tt, prob = prob_pub(tt,est.point$est$par)
) %>% 
  ggplot(aes(x=tt,y=prob)) + geom_point()

grid.arrange(p1,p2,p3,p4, nrow = 2)

# output to console
est.point

toc = Sys.time()

toc - tic

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
  , namevec = c('pif', 'truepar1','truepar2','tmid','tslope')  
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




