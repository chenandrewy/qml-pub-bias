rm(list = ls())

# ENVIRONMENT  ====

source('0_Settings.r')

# bootstrap timing for reference (about 1 min is desired)
avail_hours = 12
nboot_wanted = 1000
min_per_boot = 12*60/nboot_wanted
min_per_boot

# ESTIMATION SETUP ====

# name for documentation
bootname = 'lognorm-raw'

# filter 
#   remove 24 tabs < 1.96.  Modeling these does not pass cost / benefit
cz_filt = import_cz(dl=F) %>% filter(tabs > 1.96)

# bootstrap
nboot = 1000 # use 1 to just estimate the point
nsignal = nrow(cz_filt)

# estimation settings 
#   opt_method= 'two-stage' or 'crs-only' or 'pif-grid'
set.boot = list(
  model_fam = fam.lognormraw
  , est_names = c('pif','mua','siga','pubpar2')
  , opt_method = 'two-stage'
)

# initialize bootstrap
id_all = matrix(
  sample(1:nsignal, size = nsignal*nboot, replace = T), nrow = nsignal 
)
id_all[,1] = 1:nsignal # fix first column is always the empirical data

# initialize random guesses
par.guess.all = random_guess(set.boot, nboot, seed = 155)

# BOOTSTRAP ====

start_time = Sys.time()
bootpar = tibble()
bootstat = tibble()
for (booti in 1:nboot){
  
  tic = Sys.time()
  
  ## estimate one ====
  # booti = 3, 5, 6, 13
  
  bootsamp = cz_filt[id_all[,booti], ]
  par.guess = par.guess.all[booti,]
  
  est = estimate(set.boot,bootsamp$tabs,par.guess)
  stat = make_stats(est$par)
  
  
  # store ==== 
  bootpar = rbind(bootpar, est$par %>% mutate(booti = booti))
  bootstat = rbind(bootstat, stat)
  
  
  if (booti%%50 == 0){
    filename = paste0('boot ', bootname, ' started ', start_time)
    filename = gsub('[:]', '-', filename)
    filename = substr(filename, 1, nchar(filename)-3)
    save.image(paste0('intermediate/', filename, '.Rdata'))
  } # if booti
  
  # feedback
  print(paste0('booti = ', booti))
  toc = Sys.time()
  print(toc - tic)
  print(est$par)
  
  p.fit = hist_emp_vs_par(bootsamp$tabs,est$par) +
    ggtitle(
      paste0(
        'pif = ', round(est$par$pif,2), ', ', substr(est$opt$message, 7, 20)
        )
    )
  p.hurdle = bootstat %>% filter(stat == 'fdr5') %>% 
    ggplot(aes(x=value)) + geom_histogram(aes(y=stat(density))) + theme_minimal() +
    xlab('t-hurdle 5%')
  p.pif = bootpar %>% ggplot(aes(x=pif)) + geom_histogram(aes(y=stat(density))) + theme_minimal()
  
  grid.arrange(p.fit,p.pif,p.hurdle,nrow = 1)  
  

} # for booti

