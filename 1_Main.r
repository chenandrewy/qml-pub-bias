rm(list = ls())

# ENVIRONMENT  ====

source('0_Environment.r')

# bootstrap timing for reference
avail_hours = 12
nboot_wanted = 1000
min_per_boot = 12*60/nboot_wanted
min_per_boot

# IMPORT PUBCROSS ====

pubcross = import_cz(dl = F)

# ESTIMATION SETUP ====

# filter 
#   remove 24 tabs < 1.96.  Modeling these may not pass cost / benefit
samp = pubcross %>% filter(tabs > 1.96)

# bootstrap
nboot = 100
nsignal = nrow(samp)

# predraw indexes
set.seed(833)
id_all = matrix(
  sample(1:nsignal, size = nsignal*nboot, replace = T), nrow = nsignal 
)
id_all[,1] = 1:nsignal # fix first column is always the empirical data


# model parameters
par.guess = data.frame(
  pif = 0.1
  , pia = 0.25
  , mua = 2
  , siga = 0.1
  , mub = 5
  , sigb = 2
  , pubfam = 'stair' # 'stair' or 'piecelin'
  , pubpar1 = NA
  , pubpar2 = 0.75
)
par.lb = data.frame(
  pif = 0.001
  , pia = 0.001
  , mua = 0.5
  , siga = 0.001
  , mub = 1
  , sigb = 0.0001
  , pubpar1 = NA
  , pubpar2 = 0.3
)
par.ub = data.frame(
  pif = 0.99
  , pia = 0.99
  , mua = 4
  , siga = 10
  , mub = 10
  , sigb = 10
  , pubpar1 = NA  
  , pubpar2 = 1
)


# parameters to estimate
namevec = c('pif', 'pia','mua','mub','siga','sigb','pubpar2')


# opts
opts = list(
  algorithm = 'NLOPT_GN_CRS2_LM' # 'NLOPT_LN_BOBYQA' or 'NLOPT_GN_DIRECT_L' or 'NLOPT_GN_CRS2_LM'
  , xtol_rel = 0.001
  , print_level = 0
  , maxeval = 20
)

## checks ====

# test objective
negloglike(par2parvec(par.guess,namevec))

# check pub prob
tt = seq(1,3,0.001)
tibble(
  tt = tt, prob = pub_prob(tt,par.guess)
) %>% 
  ggplot(aes(x=tt,y=prob)) + geom_line() +theme_minimal()

# add more later
namevec

opts


# ESTIMATE ====

for (booti in 1:nboot){
  
  # optimize
  opt = nloptr(
    x0 = par2parvec(par.guess, namevec)
    , lb = par2parvec(par.lb, namevec)
    , ub = par2parvec(par.ub, namevec)
    , eval_f = negloglike
    , opts = opts
  )
  
  # store
  parvec.hat = opt$solution
  est = parvec2par(parvec.hat,par.guess,namevec)
  est$loglike = -1*opt$objective
  
  if (booti == 1){
    bootdat = est
  } else {
    bootdat = rbind(bootdat,est)
  } 
  
  # feedback
  print(booti)
  print(est)
  
  p1 = bootdat %>% ggplot(aes(x=pif)) + geom_histogram() + theme_minimal()
  p2 = hist_emp_vs_par(est)
  
  grid.arrange(p1,p2,nrow = 1)

} # for booti

