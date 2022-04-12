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
nboot = 1000 # use 1 to just estimate the point

# estimation settings 
#   copied over from 1_CheckSet.r after we find that it works
set.boot = list(
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
  , boot_method = 'semipar' # simple or semipar
)

# alt estimation settings
set.alt = set.boot
set.alt$boot_method = 'simple'
nbootalt = 1000

# BOOTSTRAP FULL ====

bootdat = bootstrap(cz_filt$tabs, set.boot, nboot, paste0(bootname, ' semipar'))


# BOOTSTRAP alt  ====

altdat = bootstrap(cz_filt$tabs, set.alt, nbootalt, paste0(bootname, ' simple'))
