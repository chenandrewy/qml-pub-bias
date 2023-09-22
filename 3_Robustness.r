# one script to run all robustness tests 2022 04 08
# outputs Rdata to intermediate/ 
# does not overwrite unless run twice in the same day

# 20 sec per bootstrap, 3 hours per model, 30 hours total

# SETUP  ====
rm(list = ls())
source('0_Settings.r')

# ESTIMATION SETUP ====

# filter or alternative filter
# keeping tabs < 0 (like tabs = 0.05) leads to hard to interpret shrinkage (dividing by 0.05)
cz_filt = import_cz(dl=F) %>% filter(tabs > 1.96)
cz_filt_less = import_cz(dl=F)  %>% filter(tabs > 0.5)

# bootstrap
nbootrobust = 500 

# raw t-stat, staircase pub prob (not really raw anymore)
set.raw_stair2 = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200)
  , opt_list2 = opts.qa(xtol_rel = 1e-3)  
  , model_fam = data.frame(
    mufam   = 'lognormraw' # 'mix-norm' or 't' or 'lognorm'
    , pif     = c(NA, 0.01, 0.99)
    , mua     = c(NA,    0, 2) # mua = 0 => median = 1.0
    , siga    = c(NA, 0.05, 1) # siga > 1 leads to crazy variances
    , pubfam  = 'stair2' # two step stair
    , pubpar1 = c(NA,  0 , 1) # prob for t in 2 to 2.6
    , pubpar2 = c(NA,  0 , 1) # prob for 1.5 to 2
    , pubpar3 = c(NA,  0 , 1) # prob for 0 to 1.5
    , row.names  = c('base','lb','ub')  
  )  
  , name = 'robust-raw-stair2'  
  , boot_method = 'simple' # simple or semipar  
)


# raw t-stat, staircase pub prob, restricted (not really raw anymore)
set.raw_stair2_restrict = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200)
  , opt_list2 = opts.qa(xtol_rel = 1e-3)  
  , model_fam = data.frame(
    mufam   = 'lognormraw' # 'mix-norm' or 't' or 'lognorm'
    , pif     = c(NA, 0.01, 0.99)
    , mua     = c(NA,    0, 2) # mua = 0 => median = 1.0
    , siga    = c(NA, 0.05, 1) # siga > 1 leads to crazy variances
    , pubfam  = 'stair2' # two step stair
    , pubpar1 = c(NA,  1/3 , 2/3) # prob for t in 2 to 2.6
    , pubpar2 = c(NA,  0   , 1/3) # prob for 1.5 to 2
    , pubpar3 = c(NA,  0   , 1/3) # prob for 0 to 1.5
    , row.names  = c('base','lb','ub')  
  )  
  , name = 'robust-raw-stair2-restrict'  
  , boot_method = 'simple' # simple or semipar  
)



# relax the pubpar constraint
set.relax_pubpar = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200)
  , opt_list2 = opts.qa(xtol_rel = 1e-3)  
  , model_fam = data.frame(
    mufam   = 'lognormraw' # 'mix-norm' or 't' or 'lognorm'
    , pif     = c(NA, 0.01, 0.99)
    , mua     = c(NA,    0, 2) # mua = 0 => median = 1.0
    , siga    = c(NA, 0.05, 1) # siga > 1 leads to crazy variances
    , pubfam  = 'stair' # 'stair' or 'piecelin' or 'trunc'  
    , pubpar1 = c(NA,  1/3, 1)
    , row.names  = c('base','lb','ub')  
  )  
  , name = 'robust-relax_pubpar'  
  , boot_method = 'simple' # simple or semipar  
)


# pubpar = 0.5
set.pubpar_05 = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200)
  , opt_list2 = opts.qa(xtol_rel = 1e-3)  
  , model_fam = data.frame(
    mufam   = 'lognormraw' # 'mix-norm' or 't' or 'lognorm'
    , pif     = c(NA, 0.01, 0.99)
    , mua     = c(NA,    0, 2) # mua = 0 => median = 1.0
    , siga    = c(NA, 0.05, 1) # siga > 1 leads to crazy variances
    , pubfam  = 'stair' # 'stair' or 'piecelin' or 'trunc'  
    , pubpar1 = c(0.5,  0.5, 0.5)
    , row.names  = c('base','lb','ub')  
  )  
  , name = 'robust-pubpar_05'  
  , boot_method = 'simple' # simple or semipar    
)


# force pif > 10%
set.pif_10 = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200)
  , opt_list2 = opts.qa(xtol_rel = 1e-3)  
  , model_fam = data.frame(
    mufam   = 'lognormraw' # 'mix-norm' or 't' or 'lognorm'
    , pif     = c(NA, 0.1, 0.99)
    , mua     = c(NA,    0, 2) # mua = 0 => median = 1.0
    , siga    = c(NA, 0.05, 1) # siga > 1 leads to crazy variances
    , pubfam  = 'stair' # 'stair' or 'piecelin' or 'trunc'  
    , pubpar1 = c(NA,  1/3, 2/3)
    , row.names  = c('base','lb','ub')  
  )  
  , name = 'robust-pif_10'    
  , boot_method = 'simple' # simple or semipar    
)


# force pif > 20%
set.pif_20 = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200)
  , opt_list2 = opts.qa(xtol_rel = 1e-3)  
  , model_fam = data.frame(
    mufam   = 'lognormraw' # 'mix-norm' or 't' or 'lognorm'
    , pif     = c(NA, 0.2, 0.99)
    , mua     = c(NA,    0, 2) # mua = 0 => median = 1.0
    , siga    = c(NA, 0.05, 1) # siga > 1 leads to crazy variances
    , pubfam  = 'stair' # 'stair' or 'piecelin' or 'trunc'  
    , pubpar1 = c(NA,  1/3, 2/3)
    , row.names  = c('base','lb','ub')  
  )  
  , name = 'robust-pif_20'    
  , boot_method = 'simple' # simple or semipar    
)


# exponential mutrue
set.exp = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200)
  , opt_list2 = opts.qa(xtol_rel = 1e-3)  
  , model_fam = data.frame(
    mufam   = 'exp' 
    , pif     = c(NA, 0.01, 0.99)
    , mua     = c(NA, 0.1, 6) 
    , pubfam  = 'stair' # 'stair' or 'piecelin' or 'trunc'  
    , pubpar1 = c(NA,  1/3, 2/3)
    , row.names  = c('base','lb','ub')  
  )  
  , name = 'robust-exp'
  , boot_method = 'simple' # simple or semipar    
)


# t-distributed mutrue (a la CZ)
set.tdist = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200, print_level = 3)
  , opt_list2 = opts.qa(xtol_rel = 1e-3, print_level = 3)  
  , model_fam = data.frame(
    mufam   = 't' # 'mix-norm' or 't' or 'lognorm'
    , pif     = c(NA, 0.01, 0.99)
    , mua     = c(NA, 1, 10) 
    , siga    = c(NA, 0.1, 6) # fix integration roundoff error
    , nua     = c(NA, 2.5, 50)    
    , pubfam  = 'stair' # 'stair' or 'piecelin' or 'trunc'  
    , pubpar1 = c(NA,  1/3, 2/3)
    , row.names  = c('base','lb','ub')  
  )  
  , name = 'robust-tdist'    
  , boot_method = 'simple' # simple or semipar    
)

# mixture normal
set.mix_norm = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200, print_level = 3)
  , opt_list2 = opts.qa(xtol_rel = 1e-3, print_level = 3)  
  , model_fam = data.frame(
    mufam   = 'mix-norm'
    , pif     = c(NA, 0.01, 0.99)
    , pia     = c(NA, 0.001, 0.999)
    , mua     = c(NA, 1, 10) 
    , siga    = c(NA, 0.1, 8) 
    , mub     = c(NA, 1, 10) 
    , sigb    = c(NA, 0.1, 8)     
    , pubfam  = 'stair' # 'stair' or 'piecelin' or 'trunc'  
    , pubpar1 = c(NA,  1/3, 2/3)
    , row.names  = c('base','lb','ub')  
  )  
  , name = 'robust-mix_norm'    
  , boot_method = 'simple' # simple or semipar    
)


# logistic pub prob
set.logistic = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200)
  , opt_list2 = opts.qa(xtol_rel = 1e-3)  
  , model_fam = data.frame(
    mufam   = 'lognormraw' # 'mix-norm' or 't' or 'lognorm'
    , pif     = c(NA, 0.01, 0.99)
    , mua     = c(NA,    0, 2) # mua = 0 => median = 1.0
    , siga    = c(NA, 0.05, 1) # siga > 1 leads to crazy variances
    , pubfam  = 'logistic' 
    , pubpar1 = c(NA,  0 , 3) # midpoint
    , pubpar2 = c(NA,  1 , 20)    # slope
    , row.names  = c('base','lb','ub')  
  )  
  , name = 'robust-logistic'  
  , boot_method = 'simple' # simple or semipar  
)




# compile
# setlist = list(
#   set.tdist
#   , set.mix_norm  
#   , set.raw_stair2_restrict
#   , set.raw_stair2
#   , set.relax_pubpar  
#   , set.pubpar_05
#   , set.logistic  
#   , set.pif_10  
#   , set.pif_20  
#   , set.exp
# )


# debug
setlist = list(
  set.raw_stair2_restrict
  , set.raw_stair2
)
nbootrobust = 100
# end debug

namelist = sapply(setlist,"[[",'name') 
print(namelist)

unique(namelist) 

# LOOP OVER ALL SETTINGS! ====

bootlistall = list()
for (seti in 1:length(setlist)){
  
  # run one ====
  print(
    paste0('robustness seti = ', seti, ' | ', setlist[[seti]]$name)
  )
  
  # load the right tabs
  if (grepl('raw', setlist[[seti]]$name)){
    temp_cz = cz_filt_less 
  } else {
    temp_cz = cz_filt
  }
  
  # shuffle data and estimate
  #   this call saves Rdata to intermediate/boot [name] [date].RData
  bootlistall[[seti]] = bootstrap(
    temp_cz$tabs
    , setlist[[seti]]
    , nbootrobust
    , setlist[[seti]]$name
  )
  # end run one ====
  
} # for seti

