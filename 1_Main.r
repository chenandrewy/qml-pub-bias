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


# model parameters
par.g = data.frame(
  pif = 0.4
  , truefam  = 'lognorm' # see create_mutrue_object function for options
  , truepar1 = 2
  , truepar2 = 3
  , truepar3 = NA_real_
  , pubfam  = 'logistic' # 'piecelin' or 'logistic
  , tmid    = 2 # prob of pub is 0.5 at the midpoint
  , tslope  = 4
)

# make model tabs | pub object
tabspub.g.o = par2observed(par.g)

# plot
plot_2hist(
  pubcross$tabs, r(tabspub.g.o)(2e3), seq(0,15,0.5), 'emp', 'sim'
)



# ESTIMATE ====

source('0_Environment.r')

# estimation settings
namevec = c('truepar1','truepar2', 'tmid', 'tslope', 'pif') # names of parameters to estimate

est.point = list(
  parbase = par.g
  , namevec = namevec  
  , parvecguess = par2parvec(par.g, namevec)
  , parvecmin   = c(0.01, 0.01, 0.5, 1, 0.01)
  , parvecmax   = NULL
)

# estimate
est.point = one_estimate(est.point, pubcross$tabs)

# plot
tabspub.o = par2observed(est$par)
p1 = plot_2hist(
  pubcross$tabs, r(tabspub.o)(1e3), seq(0,18,1), 'data', 'model'
)

mu = create_mutrue_object(est$par)
p2 = tibble(mu = r(mu)(1e3)) %>% 
  ggplot(
    aes(x=mu)
  ) + 
  geom_histogram(aes(y=stat(density)))


grid.arrange(p1,p2, nrow = 1)

# output to console
est.point

# PLOT OBJECTIVE ====

tempest = est.point
tempest$namevec = 'pif'
tabspub = pubcross$tabs

# -1* the log likelihood
negloglike = function(parvec){
  par = parvec2par(parvec,tempest$parbase, tempest$namevec)
  
  # create observed t-stat object
  tabspub.o = par2observed(par)
  
  onelike = d(tabspub.o)(tabspub)
  onelike[onelike <= 0 ] = 1e-6 # might be a better way to do this
  return = -1*mean(log( onelike ))
  
} # end negloglike

negloglike(0.5) %>% print()

pif = seq(0.05,0.55,0.01)
ll_list = pif*NA
for (i in 1:length(pif)){
  ll_list[i] = -1*negloglike(pif[i])
}

plot(pif, ll_list)

# BOOTSTRAP ====
nboot = 10
nsignal = nrow(pubcross)
set.seed(1217)

id_all = matrix(
  sample(1:nsignal, size = nsignal*nboot, replace = T), nrow = nboot 
)
namevec = c('truepar1','truepar2', 'tmid', 'tslope', 'pif') # names of parameters to estimate

# estimation settings: use point estimate as initial guess
est.boot = list(
  parbase = est.point$par
  , namevec = namevec  
  , parvecguess = par2parvec(est.point$par, namevec)
  , parvecmin   = c(0.01, 0.01, 0.5, 1, 0.01)
  , parvecmax   = NULL
)


for (booti in 1:nboot){
  
  # bootstrap t-stats
  temp_tabs = pubcross$tabs[id_all[booti, ]]
  
  temp_est = one_estimate(est.boot, temp_tabs)
  temp_bootdat = temp_est$par %>% mutate(loglike = temp_est$loglike)
  
  if (booti == 1){
    bootdat = temp_bootdat
  } else {
    bootdat = rbind(bootdat, temp_bootdat)
  }
  
  print(bootdat)
  
} # end for booti

