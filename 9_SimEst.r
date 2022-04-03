# 2022 04 Andrew: estimations on simulated data with correlations

# Based off experiments/qml5*.ret

# ENVIRONMENT ====
source('0_Environment.r')

# function for simulating truth, but only cross-sectional
sim_pubcross_pure = function(par, nport){
  
  # cross-section of mus  
  type = runif(nport) > par$pif
  
  #   start with everything false
  mu = numeric(nport) 
  
  #   then add true stuff
  mutrue.o = UnivarMixingDistribution(
    Norm(par$mua,par$siga), Norm(par$mub,par$sigb)
    , mixCoeff = c(par$pia, (1-par$pia))
  )
  mu[type] = r(mutrue.o)(sum(type))
  
  # cross-section of observables
  #   correlated noise
  ep = numeric(nport)
  for (i in 1:nport){
    if (i==1){
      ep[i] = rnorm(1,0,1)
    } else{
      ep[i] = par$rho*ep[i-1] + sqrt(1-par$rho^2)*rnorm(1,0,1)
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
  
} # end sim_pubcross_pure


## Simulation settings ====

nport = 1000

par.truth.base = data.frame(
  pif = 0.1
  , pia = 0.6
  , mua = 0.8
  , siga = 2
  , mub = 0.6
  , sigb = 6
  , pubfam = 'trunc'
  , pubpar1 = 0
  , pubpar2 = NA
  , rho = 0.5
)

piflist = seq(0.1,0.9,length.out = 3)
rholist = seq(0.1,0.9,length.out = 3)

## Estimation settings ====
# simulations should use pubfam = 'trunc' and pubpar1 = 0 to avoid
# LN_BOBYQA is a lot better than the linear version
# for BOBYQA, xtol_rel needs to be 0.001 or smaller

# opts
opts = list(
  algorithm = 'NLOPT_LN_BOBYQA' # 'NLOPT_LN_BOBYQA' or 'NLOPT_GN_CRS2_LM'
  , xtol_rel = 0.001
  , print_level = 0
  , maxeval = 1000 # 1000 for full
)

# model family
par.base = par.truth
par.lb = data.frame(
  pif   = 0.001
  , pia = 0.001
  , mua  = 0.5
  , siga = 0.001
  , mub  = 0.5
  , sigb = 0.001
)
par.ub = data.frame(
  pif   = 0.999
  , pia = 0.999
  , mua  = 10
  , siga = 10
  , mub  = 10
  , sigb = 10
)


# parameters to estimate
namevec = c('pif', 'pia','mua','siga','mub','sigb')
#namevec = c('pif','pia','mua')


# SIMULATE ESTIMATIONS====
# simulations should use pubfam = 'trunc' and pubpar1 = 0 to avoid
# weak identification problems


parveclist = expand.grid(piflist, rholist) %>% 
  as_tibble() %>% 
  rename(pif = 1, rho = 2)

set.seed(1103)

# predraw initial guesses (uniform between par.lb and par.ub)
namelist = colnames(par.ub)
par.guess = list()
for (name in namelist){
  if (par.lb %>% pull(name) %>% is.numeric()){
    par.guess[[name]] = runif(nrow(parveclist), par.lb %>% pull(name), par.ub %>% pull(name))
  } else {
    par.guess[[name]] = par.lb %>% pull(name)
  }
} # for name
par.guess = as_tibble(par.guess)



estlist = tibble()
for (i in 1:nrow(parveclist)){
  # # debug ==== 
  # i = 4
  # # opts
  # opts = list(
  #   algorithm = 'NLOPT_LN_BOBYQA' # 'NLOPT_LN_BOBYQA' or 'NLOPT_GN_DIRECT_L' or 'NLOPT_GN_CRS2_LM'
  #   , xtol_rel = 1
  #   , print_level = 3
  #   , maxeval = 1000 # 1000 for full
  # )  
  # namevec = c('pif','pia','mua','siga','mub','sigb')
  
  
  # load parvec into truth
  par.truth = par.truth.base
  par.truth[ , names(parveclist)] = parveclist[i,]
  
  # simulate
  set.seed(106)
  samp = sim_pubcross_pure(par.truth, nport)

  # likelihood
  #   important: par.base is a global, and so is samp
  par.base = par.truth  
  negloglike = function(parvec){
    # setup
    par = parvec2par(parvec,par.base,namevec)
    
    # individual likes
    singlelikes =  make_single_likes(samp$tabs, par)
    
    singlelikes[singlelikes<=0] = 1e-6
    return = -1*mean(log(singlelikes))
    
  } # end f_obs

  # optimize
  opt = nloptr(
    x0 = par2parvec(par.guess[i, ], namevec)
    , lb = par2parvec(par.lb, namevec)
    , ub = par2parvec(par.ub, namevec)
    , eval_f = negloglike
    , opts = opts
  )
  
  # store
  parvec.hat = opt$solution
  est = parvec2par(parvec.hat,par.base,namevec)
  est$loglike = -1*opt$objective
  estlist = rbind(estlist,est)
  
  # feedback  
  p1 = hist_emp_vs_par(samp,est) + 
    ggtitle('est')
  p2 = hist_emp_vs_par(samp,par.truth) +
    ggtitle('truth')
  grid.arrange(p1,p2, nrow = 1)
  
  opt$objective
  negloglike(par2parvec(par.truth,namevec)) %>% print
  
  rbind(
    parvec.hat, parvec.truth = par2parvec(par.truth,namevec)
  ) %>% print()
  
  # ====
  
} # for i



# AGGREGATE RESULTS ====

parveclist
estlist[ , namevec]

estlist2 = estlist %>% 
  mutate(
    pi_small = if_else(mua<mub, pia, (1-pia))
    , mu_small = if_else(mua<mub, mua, mub)
    , sig_small = if_else(mua<mub, siga, sigb)    
    , mu_large = if_else(mua<mub, mub, mua)
    , sig_large = if_else(mua<mub, sigb, siga)
  ) %>% 
  mutate(
    pia = pi_small, mua = mu_small, siga = sig_small
    , mub = mu_large, sigb = sig_large
  ) %>% 
  select(
    -ends_with('small'), -ends_with('large')
  )
  

stat_long = parveclist %>% 
  rename_at(vars(everything()), function(x) paste0(x, '_truth')) %>% 
  cbind(
    estlist2 %>% select(all_of(namevec))
  ) 

stat_pif_hat = stat_long %>% 
  pivot_wider(
    id_cols = c(pif_truth,rho_truth,pif)
    , names_from = rho_truth, values_from = pif, names_prefix = 'rho_truth_'   
  )

stat_long

par.truth.base

stat_pif_hat


# ====

stop()

est
par.truth.base

par.truth.base %>% 
  select(-pubfam) %>% 
  pivot_longer(
    cols = everything(), names_to = 'name', values_to = 'truth'
  ) %>%
  left_join(
    est %>% 
      select(-pubfam) %>% 
      pivot_longer(
      cols = everything(), names_to = 'name', values_to = 'estimate'
    )
    , by = 'name'
  )

# CHECK ONE EST ====

# opts
opts = list(
  algorithm = 'NLOPT_LN_BOBYQA' # 'NLOPT_LN_BOBYQA' or 'NLOPT_GN_DIRECT_L' or 'NLOPT_GN_CRS2_LM'
  , xtol_rel = 0.001
  , print_level = 3
  , maxeval = 200 # 1000 for full
)
namevec = c('pif','pia','mua')

# choose parameters
par.truth = par.truth.base 
par.truth = par.truth %>% 
  mutate(
    pubfam = 'stair', pubpar2 = 0.5, pif = 0.4
  )

# simulate
samp = sim_pubcross_pure(par.truth, nport = 5000)
samp = samp[1:200, ]

# likelihood
par.base = par.truth
#   important: par.base is a global, and so is samp
negloglike = function(parvec){
  # setup
  par = parvec2par(parvec,par.base,namevec)
  
  # individual likes
  singlelikes =  make_single_likes(samp$tabs, par)
  
  singlelikes[singlelikes<=0] = 1e-6
  return = -1*mean(log(singlelikes))
  
} # end f_obs

# optimize
opt = nloptr(
  x0 = par2parvec(par.base, namevec)
  , lb = par2parvec(par.lb, namevec)
  , ub = par2parvec(par.ub, namevec)
  , eval_f = negloglike
  , opts = opts
)

# store
parvec.hat = opt$solution
est = parvec2par(parvec.hat,par.base,namevec)
est$loglike = -1*opt$objective
estlist = rbind(estlist,est)


# feedback  
p1 = hist_emp_vs_par(samp,est)
p2 = hist_emp_vs_par(samp,par.truth)
grid.arrange(p1,p2,nrow=1)

opt$objective
negloglike(par2parvec(par.truth,namevec)) %>% print

rbind(
  parvec.hat, parvec.truth = par2parvec(par.truth,namevec)
) %>% print()