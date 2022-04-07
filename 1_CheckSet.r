# 2022 04 Andrew: estimations on simulated data with correlations

# Based off experiments/qml5*.ret


# ENVIRONMENT ====
rm(list = ls())
source('0_Settings.r')
tic_all = Sys.time()

# SETTINGS ====

## estimation settings  ====
set.check = list(
  model_fam = fam.lognormraw
  , est_names = c('pif','mua','siga')
  , opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
)

# force no selection for the sim (due to weak id)
set.sim = set.check
set.sim$model_fam = set.sim$model_fam %>% 
  mutate(
    pubfam = 'trunc', pubpar1 = 0, pubpar2 = NA
  )
set.sim$est_names = set.sim$est_names[set.sim$est_names != 'pubpar2']


## simulation settings ====

piflist = seq(0.1,0.9,length.out = 3)
rholist = seq(0.1,0.9,length.out = 3)
nportlist = c(10000, 200)

# TEST POINT ESTIMATE ====
tic = Sys.time()

cz_filt_tabs = import_cz(dl=F) %>% filter(tabs >= 1.96) %>% pull(tabs)

est.point = estimate(
  est.set = set.check
  , tabs = cz_filt_tabs
  , par.guess = random_guess(set.check,1,1243)
  , print_level = 3
)

est.point$par %>% print()
est.point$negloglike %>% print()

p_point = hist_emp_vs_par(cz_filt_tabs,est.point$par)
p_point

toc = Sys.time()
toc - tic


# TEST ON MANY SIMULATIONS ====
# simulations should use pubfam = 'trunc' and pubpar1 = 0 to avoid
# weak identification problems

parveclist = expand.grid(nportlist, piflist, rholist) %>% 
  as_tibble() %>% 
  rename(nport = 1, pif = 2, rho = 3)

set.seed(1103)

# predraw initial guesses
par.guess.all = random_guess(set.sim, nrow(parveclist), seed = 135)


estlist = tibble()
for (i in 1:nrow(parveclist)){

  ## simulate ===
  # load parvec 
  par.truth = set.sim$model_fam %>% filter(desc == 'base')
  par.truth[ , names(parveclist[ , 2:3])] = parveclist[i, 2:3]
  nport = parveclist[i, 1] %>% as.numeric()
  
  # simulate
  set.seed(106)
  samp = sim_pubcross_ar1(par.truth, par.truth$rho, nport)

  # estimate ===
  par.guess = par.guess.all[i, ]
  est = estimate(set.sim, samp$tabs, par.guess)
  parhat = est$par
  
  # store
  estlist = rbind(estlist,parhat) 
  
  # feedback
  p1 = hist_emp_vs_par(samp$tabs,parhat) + 
    ggtitle('est')
  p2 = hist_emp_vs_par(samp$tabs,par.truth) +
    ggtitle('truth')
  grid.arrange(p1,p2, nrow = 1)
  
  rbind(
    parvec.hat = par2parvec(parhat, set.sim$est_names)
    , parvec.truth = par2parvec(par.truth, set.sim$est_names)
  ) %>% print()
  
  # restart ===
  
} # for i


# AGGREGATE RESULTS ====


parveclist
estlist[ , set.sim$est_names]

if (set.check$model_fam$mufam[1] == 'mix-norm'){
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
  
} else {
  estlist2 = estlist
} # if mufam


stat_long = parveclist %>% 
  rename_at(vars(everything()), function(x) paste0(x, '_truth')) %>% 
  cbind(
    estlist2 %>% select(all_of(set.sim$est_names))
  ) %>% 
  arrange(nport_truth, pif_truth)

stat_pif_hat = stat_long %>% 
  pivot_wider(
    id_cols = c(nport_truth, pif_truth,rho_truth,pif)
    , names_from = rho_truth, values_from = pif, names_prefix = 'rho_truth_'   
  ) %>% arrange(
    nport_truth, pif_truth
  )

stat_long

stat_pif_hat


toc_all = Sys.time()

toc_all - tic_all