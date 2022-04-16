# one script to run all robustness tests 2022 04 08

# SETUP  ====
rm(list = ls())
source('0_Settings.r')

# ESTIMATION SETUP ====

# filter or not
cz_filt = import_cz(dl=F) %>% filter(tabs > 1.96)
cz_raw = import_cz(dl=F)  

# bootstrap
nbootrobust = 100 

# mixture normal

# raw t-stat, staircase pub prob w/ cuts at 1.0 and 2.0
set.raw_stair3 = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200)
  , opt_list2 = opts.qa(xtol_rel = 1e-3)  
  , model_fam = data.frame(
    mufam   = 'lognormraw' # 'mix-norm' or 't' or 'lognorm'
    , pif     = c(NA, 0.01, 0.99)
    , mua     = c(NA,    0, 2) # mua = 0 => median = 1.0
    , siga    = c(NA, 0.05, 1) # siga > 1 leads to crazy variances
    , pubfam  = 'stair3' # two step stair alt
    , pubpar1 = c(NA,  0 , 1) # prob for t in 2 to 2.6
    , pubpar2 = c(NA,  0 , 1) # prob for 1.0 to 2
    , pubpar3 = c(NA,  0 , 1) # prob for 0 to 1.5
    , row.names  = c('base','lb','ub')  
  )  
  , name = 'robust-raw-stair3'  
  , boot_method = 'simple' # simple or semipar  
)


# raw t-stat, staircase pub prob
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



# raw t-stat, logistic pub prob
set.raw_logistic = list(
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
  , name = 'robust-raw-logistic'  
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
    , siga    = c(NA, 0.1, 8) 
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
setlist = list(
  set.raw_logistic
  , set.raw_stair2
  , set.raw_stair3
  , set.relax_pubpar  
  , set.pubpar_05
  , set.pif_20  
  , set.exp
  , set.tdist
  , set.mix_norm
  , set.logistic
)

namelist = sapply(setlist,"[[",'name') 
print(namelist)

unique(namelist) 

# BOOTSTRAP! ====

bootlistall = list()
for (seti in 1:length(setlist)){
  
  # run one ====
  print(
    paste0('robustness seti = ', seti, ' | ', setlist[[seti]]$name)
  )
  
  # load the right tabs
  if (grepl('raw', setlist[[seti]]$name)){
    temp_cz = cz_raw 
  } else {
    temp_cz = cz_filt
  }
  
  temp_cz
  
  bootlistall[[seti]] = bootstrap(
    temp_cz$tabs
    , setlist[[seti]]
    , nbootrobust
    , setlist[[seti]]$name
  )
  # end run one ====
  
} # for seti

