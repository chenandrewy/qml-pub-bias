rm(list = ls())

# ENVIRONMENT  ====

source('0_Environment.r')

# bootstrap timing for reference (about 1 min is desired)
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
nboot = 1000 # use 1 to just estimate the point
nsignal = nrow(samp)

# opts
opts = list(
  algorithm = 'NLOPT_GN_CRS2_LM' # 'NLOPT_LN_BOBYQA' or 'NLOPT_GN_DIRECT_L' or 'NLOPT_GN_CRS2_LM'
  , xtol_rel = 0.001
  , print_level = 0
  , maxeval = 1000 # 1000 for full
)

# model family
#   par.base is used in loglike for anything that isn't listed in namevec
par.base = data.frame(
  pif = 0.5
  , pia = 0.5
  , mua = 1
  , siga = 1
  , mub = 2
  , sigb = 2
  , pubfam = 'stair' # 'stair' or 'piecelin'  
  , pubpar1 = NA
  , pubpar2 = 2/3
)
par.lb = data.frame(
  pif = 0.001
  , pia = 0.001
  , mua = 0.5
  , siga = 0.001
  , mub = 1
  , sigb = 0.0001
  , pubfam = 'stair' # 'stair' or 'piecelin'  
  , pubpar1 = NA
  , pubpar2 = 1/3
)
par.ub = data.frame(
  pif = 0.99
  , pia = 0.99
  , mua = 4
  , siga = 10
  , mub = 10
  , sigb = 10
  , pubfam = 'stair' # 'stair' or 'piecelin'  
  , pubpar1 = NA  
  , pubpar2 = 1/2
)


# parameters to estimate
namevec = c('pif', 'pia','mua','mub','siga','sigb','pubpar2')


## checks (might want to add more later) ====

# test objective
test_loglike = negloglike(par2parvec(par.base,namevec)) %>% print()

# check pub prob
tt = seq(1,3,0.001)
tibble(
  tt = tt, prob = pub_prob(tt,par.base)
) %>% 
  ggplot(aes(x=tt,y=prob)) + geom_line() +theme_minimal()

namevec
opts


# ESTIMATE ====


# predraw indexes
set.seed(833)
id_all = matrix(
  sample(1:nsignal, size = nsignal*nboot, replace = T), nrow = nsignal 
)
id_all[,1] = 1:nsignal # fix first column is always the empirical data

# predraw initial guesses (uniform between par.lb and par.ub)
namelist = colnames(par.ub)
par.guess = list()
for (name in namelist){
  if (par.lb %>% pull(name) %>% is.numeric()){
    par.guess[[name]] = runif(nboot, par.lb %>% pull(name), par.ub %>% pull(name))
  } else {
    par.guess[[name]] = par.lb %>% pull(name)
  }
} # for name
par.guess = as_tibble(par.guess)

start_time = Sys.time()

for (booti in 1:nboot){
  
  tic = Sys.time()
  
  # optimize
  opt = nloptr(
    x0 = par2parvec(par.guess[booti,], namevec)
    , lb = par2parvec(par.lb, namevec)
    , ub = par2parvec(par.ub, namevec)
    , eval_f = negloglike
    , opts = opts
  )
  
  # store
  parvec.hat = opt$solution
  est = parvec2par(parvec.hat,par.base,namevec)
  est$loglike = -1*opt$objective
  
  if (booti == 1){
    bootdat = est
  } else {
    bootdat = rbind(bootdat,est)
  } 
  
  toc = Sys.time()
  
  # feedback
  print(paste0('booti = ', booti))
  print(toc - tic)
  print(est)
  
  p1 = bootdat %>% ggplot(aes(x=pif)) + geom_histogram() + theme_minimal()
  p2 = hist_emp_vs_par(est)
  
  grid.arrange(p1,p2,nrow = 1)
  
  # store
  if (booti%%100 == 0){
    filename = paste0('bootdat started ', start_time)
    filename = gsub('[:]', '-', filename)
    save.image(paste0('output/', filename, '.Rdata'))
  } # if booti
  

} # for booti

# BETA: FIND HURDLES ====

booti = 1
for (booti in 1:nrow(bootdat)){

  par = bootdat[booti,]
  
  # latent density
  t.o = UnivarMixingDistribution(
    Norm(0,1), Norm(par$mua,sqrt(par$siga^2+1)), Norm(par$mub,sqrt(par$sigb^2+1))
    , mixCoeff = c(par$pif, (1-par$pif)*par$pia, (1-par$pif)*(1-par$pia))
  )
  tabs.o = abs(t.o)

  f_obs = function(tt) d(tabs.o)(tt)*pub_prob(tt,par)/denom
  
  tlist = seq(0,5,0.1)
  fdr = 100 + numeric(length(tlist))
  for (ti in 1:length(tlist)){
    tt = tlist[ti]
    pr_discovery_false = 2*pnorm(-tt)
    pr_discovery = integrate(function(ttt) d(tabs.o)(ttt), tt, Inf)$value
    fdr[ti] =  pr_discovery_false/pr_discovery*par$pif*100
  } # end for ti
  
  hurdle = tlist[which(fdr < 5)[1]]
  stat = tibble(booti = booti, hurdle = hurdle )
  
  if (booti == 1){
    bootstat = stat
  } else{
    bootstat = rbind(bootstat,stat)
  }
  
  print(booti)
  
} # for booti


hist(bootstat$hurdle)