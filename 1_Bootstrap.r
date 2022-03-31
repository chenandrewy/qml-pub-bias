rm(list = ls())

# ENVIRONMENT  ====

source('0_Environment.r')

# bootstrap timing for reference (about 1 min is desired)
avail_hours = 12
nboot_wanted = 1000
min_per_boot = 12*60/nboot_wanted
min_per_boot

# IMPORT PUBCROSS ====

pubcross = import_cz(dl = T)

# ESTIMATION SETUP ====

# name for documentation
bootname = 'benchmark-rerun'

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
  , pubpar2 = 1/2
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
  , pubpar2 = 2/3
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
  if (booti%%1 == 0){
    filename = paste0('bootdat ', bootname, ' started ', start_time)
    filename = gsub('[:]', '-', filename)
    save.image(paste0('intermediate/', filename, '.Rdata'))
  } # if booti
  

} # for booti

