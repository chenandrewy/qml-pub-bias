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
dir.create('temp/')

# root of April 2021 release on Gdrive
pathRelease = 'https://drive.google.com/drive/folders/1I6nMmo8k_zGCcp9tUvmMedKTAkb9734R'

# trigger login
target_dribble = pathRelease %>% drive_ls() 


# GLOBAL VARIABLS ====


# FUNCTIONS ====

import_cz = function(dl = F){
  
  if (dl) {
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
    
    save(pubcross, file = 'temp/pubcross.Rdata')
    
  } # if dl 
    
  load(file = 'temp/pubcross.Rdata') # loads pubcross
  
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

# selection (assumed)
prob_pub = function(tt, par){
  (tt > 1.96 & tt <= 2.60)*par$tslope + (tt > 2.6)   
}

# latent t-stats
make_t.o = function(par){
  t.o = UnivarMixingDistribution(
    Norm(0,1), Norm(par$mua,sqrt(par$siga^2+1)), Norm(par$mub,sqrt(par$sigb^2+1))
    , mixCoeff = c(par$pif, (1-par$pif)*par$pia, (1-par$pif)*(1-par$pia))
  )
  tabs.o = abs(t.o)
  return = tabs.o
}


# likelihood
negloglike = function(parvec){
  par = parvec2par(parvec,par.guess,namevec)
  
  # latent
  tabs.o = make_t.o(par)
  
  # observable
  #   take on data
  tt = samp$tabs
  
  #   find single likelihoods
  denom = par$tslope*(p(tabs.o)(2.6)-p(tabs.o)(1.96)) + (1-p(tabs.o)(2.6))
  onelikes = d(tabs.o)(tt)*prob_pub(tt,par)/denom
  
  onelikes[onelikes<=0] = 1e-6
  return = -1*mean(log(onelikes))
  
} # end f_obs
