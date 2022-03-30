# ENVIRONMENT ====

source('0_Environment.r')

# load('intermediate/bootdat started 2022-03-27 00-27-22.Rdata') # pif in 0.33,1.0, inc
# load('intermediate/bootdat started 2022-03-26 18-23-40.Rdata') # pif in 0.33,0.5, inc
load('intermediate/bootdat started 2022-03-25 16-09-49 pif in 0.33 0.67.Rdata') # BENCHMARK


# CREATE BOOTSTAT (find custom stats) ====

fdrlist = c(10, 5, 1)

booti = 1
for (booti in 1:nrow(bootdat)){
  
  par = bootdat[booti,]
  
  # latent density
  t.o = UnivarMixingDistribution(
    Norm(0,1), Norm(par$mua,sqrt(par$siga^2+1)), Norm(par$mub,sqrt(par$sigb^2+1))
    , mixCoeff = c(par$pif, (1-par$pif)*par$pia, (1-par$pif)*(1-par$pia))
  )
  tabs.o = abs(t.o)
  
  f_obs = function(tt) d(tabs.o)(tt)*pub_prob(tt,par)/denom
  
  # find fdr
  tlist = c(1.96, seq(0,5,0.1)) %>% sort()
  pr_discovery = NA*numeric(length(tlist))
  fdr = NA*numeric(length(tlist))
  for (ti in 1:length(tlist)){
    tt = tlist[ti]
    pr_discovery[ti] = integrate(function(ttt) d(tabs.o)(ttt), tt, Inf)$value
    fdr[ti] = 2*pnorm(-tt)/pr_discovery[ti]*par$pif*100
    if (fdr[ti]<=min(fdrlist) & tlist[ti] >= 1.96){
      break
    }
  } # end for ti
  
  # find hurdles
  hurdlelist = NA*numeric(length(fdrlist))
  for (fdri in 1:length(fdrlist)){
    hurdlelist[fdri] = min(tlist[which(fdr<fdrlist[fdri])])
  }
  
  # save into stat
  stat = tibble(
    booti = booti, stat = paste0('fdr', fdrlist), value = hurdlelist
  ) %>% 
    rbind(
      tibble(
          booti = booti, stat = 'pr_tgt_2', value = pr_discovery[tlist == 1.96]
      )
    )

  # save into bootstat
  if (booti == 1){
    bootstat = stat
  } else{
    bootstat = rbind(bootstat,stat)
  }
  
  print(booti)
  
} # for booti

# COMBINE WITH BOOTDAT ====
bootall = bootdat
bootall$booti = row.names(bootdat)

bootall = bootall %>% 
  select(-c(pubfam)) %>% 
  pivot_longer(
    -booti, names_to = 'stat', values_to = 'value'
  ) %>% 
  rbind(
    bootstat
  )

# SAVE TO EXHIBITS/====
dir.create('exhibits/', showWarnings = F)
write.csv(x = bootall, file = 'exhibits/bootall.csv', row.names = F)