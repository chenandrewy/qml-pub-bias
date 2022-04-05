# LIBRARIES AND PATHS ====
library(tidyverse)

# for data download
library(data.table)
library(googledrive)
library(readxl)
library(writexl)
library(nloptr)
library(distr) 
library(ggplot2)
library(gridExtra)

dir.create('intermediate/', showWarnings = F)
dir.create('output/', showWarnings = F)


# GLOBALS ====

# root of March 2022 
pathRelease = 'https://drive.google.com/drive/folders/1O18scg9iBTiBaDiQFhoGxdn4FdsbMqGo'

## model families ====

fam.mixnorm = data.frame(
  desc = c('base','lb','ub')
  , mufam = 'mix-norm'
  , pubfam = 'stair' # 'stair' or 'piecelin'    
  , pif     = c(0.5, 0.001, 0.999)
  , pia  = c(0.7, 0.001, 0.999)
  , mua  = c(3, 0.5    , 10)
  , siga = c(3, 0.001, 10)
  , mub  = c(6, 1.0    , 10)
  , sigb = c(6, 0.001, 10)
  , pubpar1 = c(NA, NA, NA)
  , pubpar2 = c(1/2, 1/3, 2/3)
) 


fam.t = data.frame(
  desc = c('base','lb','ub')
  , mufam = 't'
  , pubfam = 'stair' # 'stair' or 'piecelin'    
  , pif     = c(0.5, 0.001, 0.999)
  , mua  = c(3, 0.001, 10)
  , siga = c(3, 0.01, 10) # siga = 0.001 leads to integration errors
  , nua  = c(4, 2.1, 1000)
  , pubpar1 = c(NA, NA, NA)
  , pubpar2 = c(1/2, 1/3, 2/3)
) 

## NLOPT settings  ====

opts.crs = function(
  print_level = 0, maxeval = 200, xtol_rel = 0.001, ranseed = 1142
){
  list(
    algorithm = 'NLOPT_GN_CRS2_LM' 
    , xtol_rel = xtol_rel
    , maxeval = maxeval
    , ranseed = ranseed
    , print_level = print_level    
  )
} # end opts.crs  


opts.qa = function(
  print_level = 0, maxeval = 1000, xtol_rel = 1e-3, ftol_rel = 1e-8
){
  list(
    algorithm = 'NLOPT_LN_BOBYQA' 
    , xtol_rel = xtol_rel
    , ftol_rel = ftol_rel
    , maxeval = maxeval
    , print_level = print_level    
  )
} # end opts.qa


# FUNCTIONS ====



## Estimation ====

# function: translate par df to parvec vector form
par2parvec = function(par, namevec){
  parvec = par %>% select(all_of(namevec)) %>% as.numeric()
  return = parvec
} # end par2parvec

parvec2par = function(parvec,parbase,namevec){
  par = parbase
  par[ , namevec] = parvec %>% c() %>% t() %>% as.numeric
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


# make latent density object with distr
make_tabs_object = function(par){
  
  # first make latent t object
  #   mixture normal case
  if (par$mufam == 'mix-norm'){
    
    t.o = UnivarMixingDistribution(
      Norm(0,1), Norm(par$mua,sqrt(par$siga^2+1)), Norm(par$mub,sqrt(par$sigb^2+1))
      , mixCoeff = c(par$pif, (1-par$pif)*par$pia, (1-par$pif)*(1-par$pia))
    )
    
    #   generalized student's t
  } else if (par$mufam == 't'){
    
    tempscale = par$siga*sqrt((par$nua-2)/par$nua)
    mutrue.o = tempscale*Td(df = par$nua) + par$mua
    
    t.o = UnivarMixingDistribution(
      Norm(0,1), mutrue.o
      , mixCoeff = c(par$pif, (1-par$pif))
    )
    
  } # end if par$mufam
  
  # take absolute value
  tabs.o = abs(t.o)
  
} # end make_tabs_object

# make mutrue object with distr
#   should only be used in sim
#   separate from tabs object to improve accuracy (hopefully)
make_mutrue_object = function(par){
  
  
  # first make latent t object
  #   mixture normal case
  if (par$mufam == 'mix-norm'){
    
    mutrue.o = UnivarMixingDistribution(
      Norm(par$mua,par$siga), Norm(par$mub,par$sigb)
      , mixCoeff = c(par$pia, (1-par$pia))
    )
    
    #   generalized student's t
  } else if (par$mufam == 't'){
    
    tempscale = par$siga*sqrt((par$nua-2)/par$nua)
    mutrue.o = tempscale*Td(df = par$nua) + par$mua
    
  } # end if par$mufam
  
  
} # end make_tabs_object

negloglike = function(tabs,par){
  
  # individual likes
  #   latent density
  tabs.o = make_tabs_object(par)
  
  #   observed density
  denom = integrate(
    function (tt) d(tabs.o)(tt)*pub_prob(tt,par), 0, Inf
  )$value
  
  singlelikes = d(tabs.o)(tabs)*pub_prob(tabs,par)/denom
  
  # clean up and take mean    
  singlelikes[singlelikes<=0] = 1e-6
  return = -1*mean(log(singlelikes))
  
} # end f_obs

random_guess = function(est.set, nguess, seed = NULL){
  
  set.seed(seed)
  
  # first replicate base n times
  par.guess.all = est.set$model_fam %>% filter(desc == 'base') %>% select(-desc) 
  par.guess.all=  do.call("rbind", replicate(nguess, par.guess.all, simplify = FALSE))
  
  # then randomize est_names columns uniform between lb and ub
  par.lb = est.set$model_fam %>% filter(desc == 'lb')  
  par.ub = est.set$model_fam %>% filter(desc == 'ub')  
  
  for (name in est.set$est_names){
    if (par.lb %>% pull(name) %>% is.numeric()){
      par.guess.all[[name]] = runif(nguess, par.lb %>% pull(name), par.ub %>% pull(name))
    } else {
      par.guess.all[[name]] = par.lb %>% pull(name)
    }
  } # for name
  
  return = par.guess.all
  
} # end random_tuess


estimate = function(est.set, tabs, par.guess, print_level = 0){
  # est.set$opt_method= 'two-stage' or 'crs-only' or 'pif-grid'
  
  # optimization methods
  if (est.set$opt_method == 'two-stage'){
    
    # likelihood that gets optimized
    minme = function(parvec){
      
      par = parvec2par(parvec,par.guess,est.set$est_names)
      return = negloglike(tabs, par)
      
    } # end minme
    
    # test likelihood at guess
    parvecguess = par2parvec(par.guess,est.set$est_names)
    minme(parvecguess)    
    
    # optimize
    opt = nloptr(
      x0 = par2parvec(par.guess, est.set$est_names)
      , lb = par2parvec(est.set$model_fam %>% filter(desc == 'lb'), est.set$est_names)
      , ub = par2parvec(est.set$model_fam %>% filter(desc == 'ub'), est.set$est_names)
      , eval_f = minme
      , opts = opts.crs(print_level = print_level)
    )
    
    # optimize again 
    opt = nloptr(
      x0 = opt$solution
      , lb = par2parvec(est.set$model_fam %>% filter(desc == 'lb'), est.set$est_names)
      , ub = par2parvec(est.set$model_fam %>% filter(desc == 'ub'), est.set$est_names)
      , eval_f = minme
      , opts = opts.qa(print_level = print_level)
    )
    
    
    # store
    parvec.hat = opt$solution
    parhat = parvec2par(parvec.hat,par.guess,est.set$est_names)
    
    
  } else if (est.set$opt_method == 'crs-only'){
    
    # likelihood that gets optimized
    minme = function(parvec){
      
      par = parvec2par(parvec,par.guess,est.set$est_names)
      return = negloglike(tabs, par)
      
    } # end minme
    
    # test likelihood at guess
    parvecguess = par2parvec(par.guess,est.set$est_names)
    minme(parvecguess)    
    
    opt = nloptr(
      x0 = par2parvec(par.guess, est.set$est_names)
      , lb = par2parvec(est.set$model_fam %>% filter(desc == 'lb'), est.set$est_names)
      , ub = par2parvec(est.set$model_fam %>% filter(desc == 'ub'), est.set$est_names)
      , eval_f = minme
      , opts = opts.crs(print_level = print_level, maxeval = 1000)
    )        
    
    # store
    parvec.hat = opt$solution
    parhat = parvec2par(parvec.hat,par.guess,est.set$est_names)
    
    
  } else if (est.set$opt_method == 'pif-grid'){
    
    # put 0.5 first
    pif_grid = c(0.5, seq(0.05,0.95,0.05) ) %>% unique
    
    # loop
    parhat_grid = tibble()
    opt_grid = list()
    temp.par.guess = par.guess 
    for (i in 1:length(pif_grid)){

      # fix pif
      temp.par.guess = temp.par.guess %>% mutate(pif = pif_grid[i])
      
      # remove pif from stuff that is optimized directly
      temp_est_names = est.set$est_names[!est.set$est_names == 'pif']
      
      # likelihood that gets optimized, with fixed pif
      minme = function(parvec){
        
        par = parvec2par(parvec,temp.par.guess,temp_est_names)
        return = negloglike(tabs, par)
        
      } # end minme      

      # test likelihood at guess
      parvecguess = par2parvec(temp.par.guess,temp_est_names)
      minme(parvecguess)    
      
      # optimize
      temp_opt = nloptr(
        x0 = par2parvec(temp.par.guess, temp_est_names)
        , lb = par2parvec(est.set$model_fam %>% filter(desc == 'lb'), temp_est_names)
        , ub = par2parvec(est.set$model_fam %>% filter(desc == 'ub'), temp_est_names)
        , eval_f = minme
        , opts = opts.qa(print_level = print_level)
      )         
      
      # enforce bounds
      parvechat = pmax(temp_opt$solution, par2parvec(est.set$model_fam %>% filter(desc == 'lb'), temp_est_names)) 
      parvechat = pmin(parvechat, par2parvec(est.set$model_fam %>% filter(desc == 'ub'), temp_est_names))
        
      temp_parhat = parvec2par(parvechat,temp.par.guess,temp_est_names) %>% 
        mutate(negloglike = temp_opt$objective)
      
      # store
      parhat_grid = rbind(
        parhat_grid
        , temp_parhat
      )
      opt_grid[[i]] = temp_opt
      
      # update guess
      temp.par.guess = temp_parhat
      
    } # for i in pif_grid
    
    # find best
    ibest = which.min(parhat_grid$negloglike)
    parhat = parhat_grid[ibest,] %>% select(-negloglike)
    opt = opt_grid[[ibest]]
    
  } # end if est.set$opt_method
  
  
  return = list(
    par = parhat
    , negloglike = opt$objective
    , opt = opt
  )
  
} # end estimate

toc = Sys.time()


## Data generation ====

import_cz = function(dl = F){
  
  if (dl) {
    
    # download signal documentation 
    target_dribble = pathRelease %>% drive_ls() %>% 
      filter(name=='SignalDoc.csv')
    
    drive_download(target_dribble, path = 'intermediate/deleteme.csv', overwrite = T)
    
    signaldoc = fread('intermediate/deleteme.csv') %>% 
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
    
    write.csv(x = pubcross, file = 'output/pubcross.csv', row.names = F)
    
  } # if dl 
  
  pubcross = fread('output/pubcross.csv') # loads pubcross
  
  return = pubcross
  
} # end download_cz




# function for simulating truth, but only cross-sectional
sim_pubcross_ar1 = function(par, rho, nport){
  
  # cross-section of mus  
  type = runif(nport) > par$pif
  
  #   start with everything false
  mu = numeric(nport) 
  
  #   then add true stuff
  mutrue.o = make_mutrue_object(par)
  mu[type] = r(mutrue.o)(sum(type))
  
  # cross-section of observables
  #   correlated noise
  ep = numeric(nport)
  for (i in 1:nport){
    if (i==1){
      ep[i] = rnorm(1,0,1)
    } else{
      ep[i] = rho*ep[i-1] + sqrt(1-rho^2)*rnorm(1,0,1)
    }
  } # end for i
  
  # cross-sectional stats
  cross = tibble(
    type = type
    , mu = mu
    , tstat = mu + ep
    , tabs = abs(tstat)
    , pub = runif(nport) <= pub_prob(tabs, par)
  )
  
  # published subset
  pubcross = cross %>% filter(pub)
  
  return = pubcross
  
} # end sim_pubcross_ar1

## Stats and Plotting ====

make_stats = function(par, fdrlist = c(10, 5, 1)){
  
  # latent density 
  tabs.o = make_tabs_object(par)
  
  # find fdr
  tlist = c(1.96, seq(0,5,0.1)) %>% sort()
  pr_discovery = NA*numeric(length(tlist))
  fdr = NA*numeric(length(tlist))
  for (ti in 1:length(tlist)){
    tt = tlist[ti]
    pr_discovery[ti] = integrate(function(ttt) d(tabs.o)(ttt), tt, Inf)$value
    fdr[ti] = 2*pnorm(-tt)/pr_discovery[ti]*par$pif*100
    if (fdr[ti]<=min(fdrlist) & tlist[ti] >= 1.96){
      break
    }
  } # end for ti
  
  # find hurdles
  hurdlelist = NA*numeric(length(fdrlist))
  for (fdri in 1:length(fdrlist)){
    hurdlelist[fdri] = min(tlist[which(fdr<fdrlist[fdri])])
  }
  
  # save into stat
  stat = tibble(
    booti = booti, stat = paste0('fdr', fdrlist), value = hurdlelist
  ) %>% 
    rbind(
      tibble(
        booti = booti, stat = 'pr_tgt_2', value = pr_discovery[tlist == 1.96]
      )
    )
  
  return = stat
  
} # end make_stats


hist_emp_vs_par = function(tabs,par){
  
  # latent density
  tabs.o = make_tabs_object(par)
  
  # observed density
  #   find denominator for conditional density
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
  
  hemp = hist(tabs, edge, plot = F)
  
  plotme = tibble(
    group = 'emp', tabs = hemp$mids, freq = hemp$density
  ) %>% rbind(
    tibble(
      group = 'mod', tabs = hemp$mids, freq = freq    
    )
  ) %>% 
    filter(!is.na(freq))
  
  xmin = 0; xmax = 15
  plotme %>% 
    filter(tabs >= xmin, tabs <= xmax) %>% 
    ggplot(aes(x=tabs,y=freq,group=group,color=group)) +
    geom_line() +
    geom_point(size=3) +
    theme_minimal() +
    scale_color_manual(values = c('blue','red')) + 
    theme(legend.position = c(0.8,0.8)) +
    xlim(xmin,xmax)  
  
} # end hist_emp_vs_par
