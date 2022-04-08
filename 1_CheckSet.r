# 2022 04 Andrew: estimations on simulated data with correlations

# Based off experiments/qml5*.ret


# ENVIRONMENT ====
rm(list = ls())
source('0_Settings.r')
tic_all = Sys.time()

# SETTINGS ====

## estimation settings  ====
set.check = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200)
  , opt_list2 = opts.qa(xtol_rel = 1e-3)  
  , model_fam = data.frame(
    mufam   = 'lognormraw' # 'mix-norm' or 't' or 'lognorm'
    , pif     = c(NA, 0.01, 0.99)
    , mua     = c(NA,    0, 2) # mua = 0 => median = 1.0
    , siga    = c(NA, 0.05, 1) # siga > 1 leads to crazy variances
    , pubfam  = 'stair' # 'stair' or 'piecelin' or 'trunc'  
    , pubpar1 = c(NA,  1/3, 2/3)
    , row.names  = c('base','lb','ub')  
  )  
)

## simulation settings ====

piflist = seq(0.1,0.9,length.out = 3)
rholist = seq(0.1,0.9,length.out = 3)
nportlist = c(200)
n_miniboot = 10


# TEST POINT ESTIMATE ====
tic = Sys.time()

cz_filt_tabs = import_cz(dl=F) %>% filter(tabs >= 1.96) %>% pull(tabs)

est.point = estimate(
  est.set = set.check
  , tabs = cz_filt_tabs
  , par.guess = random_guess(set.check,1,1243)
  , print_level = 3
)
est.point$bias = shrink_samp(cz_filt_tabs,est.point$par) # for now assume signs don't matter

est.point$par %>% print()
est.point$negloglike %>% print()

p_point = hist_emp_vs_par(cz_filt_tabs,est.point$par)
p_bias = est.point$bias %>% as_tibble() %>% ggplot(aes(x=value)) + geom_histogram() +
  xlab('bias')

grid.arrange(p_point,p_bias, nrow = 1)

toc = Sys.time()
toc - tic

# TEST ONE SIMULATION ====
nshrink = 1000

temp.truth = est.point$par %>% mutate(pif = 0.9)
  
sim.one = sim_pubcross_ar1(temp.truth, rho = 0.2, nport = 10000)

nshrink = min(dim(sim.one)[1], nshrink)

est.sim.one = estimate(
  est.set = set.check 
  , tabs = sim.one$tabs
  , par.guess = random_guess(set.check,1,1243)
  , print_level = 3
)

est.sim.one$bias = shrink_samp(sim.one$tstat[1:nshrink],est.point$par) 

rbind(
  temp.truth %>% mutate(group = 'truth')
  ,est.sim.one$par  %>% mutate(group = 'est')
) 

# plot
p_point = hist_emp_vs_par(sim.one$tabs,est.sim.one$par)
plotme = cbind(
  sim.one[1:nshrink,], est.sim.one$bias %>% as_tibble() %>% rename(bias = value)
) %>% 
  mutate(
    muhat = tstat*(1-bias)
  ) %>% 
  as_tibble()

temptext = paste0('Ebias = ', round(mean(est.sim.one$bias),2))
p_bias = plotme %>% 
  ggplot(aes(x=muhat,y=mu)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  annotate("text", x=10,y=1,label = temptext, size = 5)

grid.arrange(p_point,p_bias, nrow = 1)



# TEST ON MANY SIMULATIONS ====

## many simulations setup ====
# simulate truth as point estimate without selection bias to 
# deal with weak identification (no need to worry about that here)
par.truth = est.point$par %>% 
  mutate(
    pubfam = 'trunc', pubpar1 = 0
  )

set.sim = set.check 
set.sim$model_fam = set.sim$model_fam %>% 
  mutate(
    pubfam = 'trunc', pubpar1 = 0
  )


parveclist = expand.grid(nportlist, piflist, rholist) %>% 
  as_tibble() %>% 
  rename(nport = 1, pif = 2, rho = 3) %>% 
  arrange(nport, pif, rho)

set.seed(1103)

# predraw initial guesses
par.guess.all = random_guess(set.sim, nrow(parveclist), seed = 135)

estlist = tibble()
truthlist = tibble()

## loop over many simulations ====
for (truthi in 1:nrow(parveclist)){

  ## simulate one ====

  # load parvec 
  par.truth[ , names(parveclist[ , 2:3])] = parveclist[truthi, 2:3]
  nport = parveclist[truthi, 1] %>% as.numeric()
  
  # simulate
  set.seed(958)
  
  for (booti in 1:n_miniboot){
    print(paste0('truthi = ', truthi, ' of ', nrow(parveclist)))
    print(paste0('booti = ', booti))
    
    samp = sim_pubcross_ar1(par.truth, par.truth$rho, nport)
    
    # estimate ===
    par.guess = par.guess.all[truthi, ]
    est = estimate(set.sim, samp$tabs, par.guess)

    # store
    estlist = rbind(
      estlist
      , est$par %>% mutate(
        negloglike = est$negloglike, booti = booti, truthi = truthi
      )
    )     
    truthlist = rbind(
      truthlist
      , par.truth %>% mutate(nport = nport)
    )
    
    # feedback
    p1 = hist_emp_vs_par(samp$tabs,est$par) + 
      ggtitle(paste0('truthi = ', truthi, ', est'))
    p2 = hist_emp_vs_par(samp$tabs,par.truth) +
      ggtitle('truth')
    grid.arrange(p1,p2, nrow = 1)
    
    print('parhat')
    print(est$par)
    print('par.truth')
    print(par.truth)    
    
  } # for booti

  
  
  # restart ===
  
} # for truthi


# AGGREGATE RESULTS ====

est_names = set.sim$model_fam['base', is.na(set.sim$model_fam['base',])] %>% colnames

# put all data together
dat_long = truthlist %>% 
  select(nport, pif, rho) %>% 
  rename_at(vars(everything()), function(x) paste0(x, '_truth')) %>% 
  cbind(
    estlist %>% select(all_of(est_names))
  ) %>% 
  arrange(nport_truth, pif_truth, rho_truth) 

# find bootstrapped stats
temp1 = dat_long %>% 
  group_by(nport_truth, pif_truth, rho_truth) %>% 
  summarize_at(
    vars(c(pif,mua,siga)), mean
  ) %>%
  mutate(
    bootstat = 'mean'
  ) 
temp2 = dat_long %>% 
  group_by(nport_truth, pif_truth, rho_truth) %>% 
  summarize_at(
    vars(c(pif,mua,siga)), sd
  ) %>%
  mutate(
    bootstat = 'sd'
  )   

stat_long = rbind(temp1,temp2) %>% 
  arrange(bootstat, nport_truth, pif_truth, rho_truth)


bootmean_pif = stat_long %>% 
  filter(bootstat == 'mean') %>% 
  pivot_wider(
    id_cols = c(nport_truth, pif_truth,rho_truth,pif)
    , names_from = rho_truth, values_from = pif, names_prefix = 'rho_truth_'   
  ) %>% arrange(
    nport_truth, pif_truth
  )

bootsd_pif = stat_long %>% 
  filter(bootstat == 'sd') %>% 
  pivot_wider(
    id_cols = c(nport_truth, pif_truth,rho_truth,pif)
    , names_from = rho_truth, values_from = pif, names_prefix = 'rho_truth_'   
  ) %>% arrange(
    nport_truth, pif_truth
  )

# output to console
stat_long
bootmean_pif
bootsd_pif

toc_all = Sys.time()

toc_all - tic_all