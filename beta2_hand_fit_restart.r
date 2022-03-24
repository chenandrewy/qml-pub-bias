# 2022 03 24 testing hand fitting of data to understand effect of parameteric assumptions
# was running into mysterious bugs, so starting over

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


# testing: flexible qml fit stair w/ general triple mix normal  ====

# bootstrap or not
nport = nrow(pubcross)
id = 1:nport # no boot
id = sample(1:nport, nport, replace = T) # boot
samp = pubcross[id,]

hist(samp$tabs)



## estimate! ====
namevec = c('pif', 'pia','mua','mub','siga','sigb')

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

# latent 
tabs.o = make_t.o(par)

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