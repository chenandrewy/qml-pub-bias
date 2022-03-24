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


# FUNCTIONS ====
# translate par0$true* into a distr object
#   must be a less ugly way to do this but I dunno what it is
#   to find a list of distributions, use ??distr or see page 7 of long manual (distribution classes)
create_mutrue_object = function(par){
  
  if (par$truefam == 'norm'){
    
    str.o = 'Norm(par$truepar1, par$truepar2)'
    
  } else if (par$truefam == 'lognorm'){
    
    # truepar1 = E(mu) and truepar2 = sd(mu), use wikipedia formulas to translate
    mu = log( par$truepar1^2/sqrt(par$truepar2^2+par$truepar1^2) )
    sig = sqrt(log(1+par$truepar2^2/par$truepar1^2))
    str.o = paste0('exp(Norm(mu,sig))')
    
  } else if (par$truefam == 'foldedt'){
    
    str.o = 'abs(par$truepar1 + par$truepar2*Td(par$truepar3))'
    
  } else if (par$truefam == 'exp'){
    str.o = 'Exp(1/par$truepar1)'
  }
  
  return = eval(parse(text = str.o))
}

# prob of publication
prob_pub = function(tt,par){
  
  if (par$pubfam == 'piecelin'){
    
    prob = 0.5 + par$tslope*(tt-par$tmid)
    prob[prob<0] = 0; prob[prob>1] = 1
    return = prob    
    
  } else if (par$pubfam == 'logistic'){
    
    # this family may be required if mutrue is student's t
    prob = 1/(1+exp(-par$tslope*(tt-par$tmid)))
    
  } else if (par$pubfam == 'stair'){
    # does not work
    prob = (tt > 1.96 & tt <= 2.58)*par$tslope +(tt>2.58)*1
    # prob = (tt > 1.96)*1
    
  } else if (par$pubfam == 'trunc'){
    prob = 1*(tt > par$tmid)
    
  } # end if par$pubfam
  
  
} # end prob_pub


# function for creating object of |t| cond on publication
par2observed = function(par){
  # create tabspub.o = |t| cond on publication in steps
  #   1) make all t
  t.o =UnivarMixingDistribution(
    Norm(0,1), create_mutrue_object(par) + Norm(0,1)
    , mixCoeff = c(par$pif, 1-par$pif)
  ) 
  #   2) take abs
  tabs.o = abs(t.o)
  
  #   3) condition on publication
  # doc here: https://www.rdocumentation.org/packages/distr/versions/2.8.0/topics/AbscontDistribution
  
  if (par$pubfam != 'trunc'){
    tabspub.o = AbscontDistribution(
      d=function(tt) d(tabs.o)(tt)*prob_pub(tt,par)
      , withStand = 1
      , low = par$tmin
      , low1 = par$tmin
    )
    
  } else {
    tabspub.o = Truncate(tabs.o, lower = par$tmid)
  } # if par$pubfam 
  
  return = tabspub.o
  
} # par2observed

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

# lower and upper boudns of search space
parveclim = function(namevec){
  par.lb = data.frame(
    pif = 0.001, truepar1 = 0, truepar2 = 0, truepar3 = 0, tmid = 0, tslope = 0
  )
  par.ub = data.frame(
    pif = 1, truepar1 = 100, truepar2 = 100, truepar3 = 100, tmid = 3, tslope = 100
  )
  
  out = list(
    lb = par.lb %>% select(all_of(namevec)) %>% as.numeric()
    , ub = par.ub %>% select(all_of(namevec)) %>% as.numeric()
  )
  
  return = out
} # end parveclim

# function for turning two series into a histogram plot.  
# This should be default in r but why is this so hard
plot_2hist = function(x1,x2,edge = seq(0,15,0.5), x1name = 'x1', x2name = 'x2'){
  datall = tibble(
    group = x1name, x = x1
  ) %>% rbind(
    tibble(
      group = x2name, x = x2
    )
  )
  
  hall = datall %>% 
    filter(x>min(edge), x<max(edge)) %>% 
    group_by(group) %>% 
    summarise(
      tmid = hist(x,edge)$mid
      , density = hist(x,edge)$density
    ) 
  
  hall %>% 
    ggplot(aes(x=tmid, y=density)) +
    geom_bar(
      aes(fill = group), stat = 'identity', position = 'dodge'
    )
  
} # end plot 2 hist




one_estimate = function(est,tabspub,opts=opts.base){
  
  # -1* the log likelihood
  negloglike = function(parvec){
    par = parvec2par(parvec,est$parguess,est$namevec)
    
    # create observed t-stat object
    tabspub.o = par2observed(par)
    
    onelike = d(tabspub.o)(tabspub)
    onelike[onelike <= 0 ] = 1e-6 # might be a better way to do this
    out = -1*mean(log( onelike ))
    return = out
    
  } # end negloglike
  
  parvecguess = par2parvec(est$parguess, est$namevec)
  
  # optimize
  opt = nloptr(
    x0 = parvecguess
    , lb = parveclim(est$namevec)$lb
    , ub = parveclim(est$namevec)$ub
    , eval_f = negloglike
    , opts = opts
  )
  
  # pack output nicely along with settings
  est$loglike = -opt$objective
  est$par    = parvec2par(opt$solution, est$parguess, est$namevec) 
  
  return = est
  
} # end one_estimate

# estimate pi_f on a grid
pifgrid_estimate = function(est,tabspub,opts = opts.base){
  
  # estimate at all pifs on a grid, coarsely
  for (i in 1:length(est$pifgrid)){
    
    # modify one_estimate settings
    tempest = est
    tempest$namevec = tempest$namevec[tempest$namevec != 'pif'] # remove pif from namevec
    tempest$parguess$pif = est$pifgrid[i] # replace pif guess with current grid
    
    tempest = one_estimate(tempest,tabspub,opts)
    temppar = tempest$par %>% 
      mutate(loglike = tempest$loglike)
    
    if ( i == 1 ){
      estgrid = temppar
    } else {
      estgrid = rbind(estgrid,temppar)
    } # if i == 1
    
    # feedback
    # temppar %>% print()
    
  } # for i
  
  # find the pif that has the best loglike and estimate again, but finely
  ibest = which.max(estgrid$loglike)
  tempest = est
  tempest$parguess = estgrid[ibest, ] %>% select(-loglike)
  estfinal = one_estimate(tempest,tabspub,opts)
  
  # pack up and output
  out = list(
    grid = estgrid
    , est = estfinal
  )
  return = out
  
} # end grid_estimate


# calculate t-hurdle from model
find_hurdle = function(par){
  hlist = seq(0,5,0.01)
  
  t.o =UnivarMixingDistribution(
    Norm(0,1), create_mutrue_object(par) + Norm(0,1)
    , mixCoeff = c(par$pif, 1-par$pif)
  ) 
  
  fdrdat = tibble(
    h = hlist
    , prob_h_f = 2*pnorm(-hlist)
    , prob_h = 1-p(t.o)(hlist)
    , fdr_h = prob_h_f/prob_h*par$pif
  )    
  
  hstar = fdrdat %>% filter(fdr_h <= 0.05) %>% summarize(hstar = min(h)) %>% pull(hstar)
  
  return = hstar
  
} # end find_hurdle

plot_emp_vs_mod = function(mod,emptabs){
  
  edge = c(seq(0,15,1), 40)
  freq = numeric(length(edge)-1)
  for (i in 2:length(edge)){
    a = edge[i-1]; b = edge[i]
    freq[i-1] = p(mod)(b) - p(mod)(a)
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
  
} # end plot_emp_vs_mod
