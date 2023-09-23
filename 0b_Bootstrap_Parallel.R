# bootstrap parallel

# seems to require its own file so that I can source 0_Settings inside the foreach loop
bootstrap_parallel = function(tabs,set.boot,nboot,ncores = 4){
  
  nsignal = length(tabs)
  
  # initialize bootstrap
  set.seed(1053)
  id_all = matrix(
    sample(1:nsignal, size = nsignal*nboot, replace = T), nrow = nsignal 
  )
  id_all[,1] = 1:nsignal # fix first column is always the empirical data
  
  # initialize random guesses
  par.guess.all = random_guess(set.boot, nboot, seed = 155)  
  
  start_time = Sys.time()
  cl = makePSOCKcluster(ncores)
  registerDoParallel(cl)
  bootdat = foreach(booti = 1:nboot
                    , .combine = rbind
                    , .packages = c('tidyverse','nloptr','distr')) %dopar% {
                      
                      source('0_Settings.r')
                      
                      print(paste0('booti = ', booti))    
                      
                      tic = Sys.time()
                      
                      # estimate one ===
                      
                      # if first round, use empirical data
                      if (booti == 1){
                        boottabs = tabs
                        
                        # if second round, do either simple boot  
                      } else if (set.boot$boot_method == 'simple'){
                        
                        boottabs = tabs[id_all[,booti]]
                        
                        
                        # or semi-parameteric boot
                      } else if (set.boot$boot_method == 'semipar') {
                        
                        simcross = sim_pubcross(par = bootpar[1, ], nport = 5000, ndate = 350, seed = id_all[1, booti])
                        tempn = min(dim(simcross)[1], length(tabs)) # use nport that is smaller of simcross or emp nport
                        boottabs = simcross$tabs[1:tempn]
                        
                      } # end if booti == 1
                      
                      par.guess = par.guess.all[booti,]    
                      est = estimate(set.boot,boottabs,par.guess)
                      
                      # make stats
                      stat = make_stats(est$par)
                      temp = make_stats_pub(boottabs, est$par)
                      stat = stat %>% 
                        mutate(
                          bias_mean = mean(temp$bias)
                          , bias_med = median(temp$bias)
                          , fdrloc_mean = mean(temp$fdr_tabs)
                          , fdrloc_med = median(temp$fdr_tabs)        
                          , sbias_mean = mean(temp$bias_signed)
                          , sbias_med = median(temp$bias_signed)        
                        )
                      
                      # feedback
                      toc = Sys.time()
                      print(toc - tic)
                      print(est$par)
                      
                      # store === 
                      bootdat_one = cbind(est$par, stat) %>% 
                        mutate(booti = booti)
                      
                      return(bootdat_one)
                      
                    } # for booti
  
  # separate parameters from stats (results)
  name_par = names(par.guess.all)
  name_stat = setdiff(names(bootdat), c(name_par, 'booti'))
  
  boot = list(
    par = bootdat %>% select(booti, all_of(name_par)) %>% setDT()
    , stat = bootdat %>% select(booti, all_of(name_stat)) %>% setDT()
  )
  
  return(boot)
  
} # end function bootstrap
