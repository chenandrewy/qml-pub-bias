# one script to run all robustness tests 2022 04 08

# SETUP  ====
rm(list = ls())
source('0_Settings.r')

# ESTIMATION SETUP ====

# name for documentation
bootname = 'a-priori-model'

# filter 
#   remove 24 tabs < 1.96.  Modeling these does not pass cost / benefit
cz_filt = import_cz(dl=F) %>% filter(tabs > 1.96)

# bootstrap
nbootrobust = 200 

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


# compile
setlist = list(
  set.relax_pubpar  
  , set.pubpar_05
  , set.pif_20  
  , set.exp
  , set.tdist
)

# BOOTSTRAP! ====

bootlistall = list()
for (seti in 1:length(setlist)){
  
  # run one ====
  print(
    paste0('robustness seti = ', seti, ' | ', setlist[[seti]]$name)
  )
  seti = 5
  bootlistall[[seti]] = bootstrap(
    cz_filt$tabs
    , setlist[[seti]]
    , nbootrobust
    , setlist[[seti]]$name
  )
  # end run one ====
  
} # for seti

