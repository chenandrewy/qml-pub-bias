# ENVIRONMENT ====
library(tidyverse)

# for data download
library(data.table)
library(googledrive)
library(readxl)

# for estimation
library(nloptr)
opts = list(
  'algorithm'='NLOPT_LN_COBYLA'
  , 'xtol_rel' = 1e-2
  , 'ftol_rel' = 1e-1
  , print_level = 3
  , maxtime = 20
)
library(distr) 
library(ggplot2)
library(gridExtra)
dir.create('temp/')

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
  tabspub.o = AbscontDistribution(
    d=function(tt) d(tabs.o)(tt)*prob_pub(tt,par)
    , withStand = 1
    , low = 0
  )
  
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




one_estimate = function(est,tabspub){
  
  # -1* the log likelihood
  negloglike = function(parvec){
    par = parvec2par(parvec,est$parbase,est$namevec)
    
    # create observed t-stat object
    tabspub.o = par2observed(par)
    
    onelike = d(tabspub.o)(tabspub)
    onelike[onelike <= 0 ] = 1e-6 # might be a better way to do this
    return = -1*mean(log( onelike ))
    
  } # end negloglike
  
  # optimize
  opt = nloptr(
    x0 = est$parvecguess
    , lb = est$parvecmin
    , ub = est$parvecmax
    , eval_f = negloglike
    , opts = opts
  )
  
  # pack output nicely along with settings
  est$parvec = opt$solution
  est$loglike = -opt$objective
  est$par    = parvec2par(est$parvec, est$parbase, est$namevec) 
  
  return = est
  
} # end one_estimate