# ENVIRONMENT ====
library(tidyverse)

# for data download
library(data.table)
library(googledrive)
library(readxl)

# for estimation
# 'algorithm'='NLOPT_LN_COBYLA' or 'NLOPT_LN_BOBYQA'

library(nloptr)
opts.base = list(
  'algorithm'='NLOPT_LN_BOBYQA'
  , 'xtol_rel' = 1e-2
  , 'ftol_rel' = 1e-1
  , print_level = 3 # 0, 1, 2, 3
  , maxtime = 20
)
library(distr) 
library(ggplot2)
library(gridExtra)

dir.create('intermediate/', showWarnings = F)
dir.create('exhibits/', showWarnings = F)

# root of April 2021 release on Gdrive
pathRelease = 'https://drive.google.com/drive/folders/1I6nMmo8k_zGCcp9tUvmMedKTAkb9734R'

# trigger login
# target_dribble = pathRelease %>% drive_ls() 


# GLOBAL VARIABLS ====


# FUNCTIONS ====

import_cz = function(dl = F){
  
  if (dl) {
    # download signal documentation 
    target_dribble = pathRelease %>% drive_ls() %>% 
      filter(name=='SignalDocumentation.xlsx')
    
    drive_download(target_dribble, path = 'intermediate/deleteme.xlsx', overwrite = T)
    
    signaldoc = left_join(
      read_excel('intermediate/deleteme.xlsx',sheet = 'BasicInfo')  
      , read_excel('intermediate/deleteme.xlsx',sheet = 'AddInfo')
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
    
    drive_download(target_dribble[1,], path = 'intermediate/deleteme2.csv', overwrite = T)
    
    ret = fread('intermediate/deleteme2.csv') %>% 
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
    
    write.csv(x = pubcross, file = 'exhibits/pubcross.csv', row.names = F)
    
  } # if dl 
    
  pubcross = fread('exhibits/pubcross.csv') # loads pubcross
  
  return = pubcross
  
} # end download_cz

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

# prob of publication
pub_prob = function(tt,par){
  
  if (par$pubfam == 'piecelin'){
    
    prob = 0.5 + par$pubpar2*(tt-par$pubpar1)
    prob[prob<0] = 0; prob[prob>1] = 1
    return = prob    
    
  } else if (par$pubfam == 'logistic'){
    
    # this family may be required if mutrue is student's t
    prob = 1/(1+exp(-par$pubpar2*(tt-par$pubpar1)))
    
  } else if (par$pubfam == 'stair'){
    # does not work
    prob = (tt > 1.96 & tt <= 2.58)*par$pubpar2 +(tt>2.58)*1
    # prob = (tt > 1.96)*1
    
  } else if (par$pubfam == 'trunc'){
    prob = 1*(tt > par$pubpar1)
    
  } # end if par$pubfam
  
  
} # end pub_prob




make_single_likes = function(tt,par){
  
  # latent density
  t.o = UnivarMixingDistribution(
    Norm(0,1), Norm(par$mua,sqrt(par$siga^2+1)), Norm(par$mub,sqrt(par$sigb^2+1))
    , mixCoeff = c(par$pif, (1-par$pif)*par$pia, (1-par$pif)*(1-par$pia))
  )
  tabs.o = abs(t.o)
  
  # observed density
  #   find denominator for conditional density
  # denom = par$pubpar2*(p(tabs.o)(2.6)-p(tabs.o)(1.96)) + (1-p(tabs.o)(2.6))
  denom = integrate(
    function (tt) d(tabs.o)(tt)*pub_prob(tt,par), 0, Inf
  )$value
  return = d(tabs.o)(tt)*pub_prob(tt,par)/denom
  
} # individual likes



# likelihood
#   important: par.base is a global
negloglike = function(parvec){
  # setup
  par = parvec2par(parvec,par.base,namevec)
  
  # individual likes
  singlelikes =  make_single_likes(samp$tabs, par)
  
  singlelikes[singlelikes<=0] = 1e-6
  return = -1*mean(log(singlelikes))
  
} # end f_obs



hist_emp_vs_par = function(par){
  
  # observable
  # this is slow
  # f_obs = function(tt){make_single_likes(tt,par)} 
  
  # this is faster
  # latent density
  t.o = UnivarMixingDistribution(
    Norm(0,1), Norm(par$mua,sqrt(par$siga^2+1)), Norm(par$mub,sqrt(par$sigb^2+1))
    , mixCoeff = c(par$pif, (1-par$pif)*par$pia, (1-par$pif)*(1-par$pia))
  )
  tabs.o = abs(t.o)
  
  # observed density
  #   find denominator for conditional density
  # denom = par$pubpar2*(p(tabs.o)(2.6)-p(tabs.o)(1.96)) + (1-p(tabs.o)(2.6))
  denom = integrate(
    function (tt) d(tabs.o)(tt)*pub_prob(tt,par), 0, Inf
  )$value
  
  f_obs = function(tt) d(tabs.o)(tt)*pub_prob(tt,par)/denom
  
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
  
  # stat = plotme %>% 
  #   pivot_wider(names_from = group, values_from = freq) %>% 
  #   mutate(
  #     error = 1000*(mod - emp)
  #   ) 
  
  # print(stat)
  # print(mean(abs(stat$error)))
  
  plotme %>% 
    ggplot(aes(x=tabs,y=freq,group=group,color=group)) +
    geom_line() +
    geom_point(size=3) +
    theme_minimal() +
    scale_color_manual(values = c('blue','red')) + 
    theme(legend.position = c(0.8,0.8)) +
    xlim(0,15)  
  
} # end hist_emp_vs_par