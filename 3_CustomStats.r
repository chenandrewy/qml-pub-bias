# 2022 03 converts parameters saved in bootpar to some custom stats and then
# saves to disk for use in 3_Exhibits.r
rm(list = ls())
# ENVIRONMENT ====

source('0_Settings.r')


# load('intermediate/bootdat started 2022-03-25 16-09-49 pif in 0.33 0.67.Rdata') # first benchmark, restricted
# bootpar = bootdat
# bootpar$mufam = 'mix-norm'

# load('intermediate/bootpar benchmark-rerun started 2022-04-01 15-28-20.Rdata') # benchmark with no restrictions
# load('intermediate/bootpar mu-lb-0.5 started 2022-04-02 14-41-32.Rdata')
# load('intermediate/bootpar two-stage started 2022-04-03 16-39.Rdata')

# set.boot$model_fam
# set.boot$est_names

# single-norm mua >= 2.0
load('intermediate/boot single-norm-trunc started 2022-04-06 14-46.Rdata')

# CREATE BOOTSTAT (find custom stats) ====


for (booti in 1:nrow(bootpar)){
  
  stat = make_stats(bootpar[booti,])

  # save into bootstat
  if (booti == 1){
    bootstat = stat
  } else{
    bootstat = rbind(bootstat,stat)
  }
  
  print(booti)
  
} # for booti

# COMBINE WITH bootpar ====
bootall = bootpar
bootall$booti = row.names(bootpar)

bootall = bootall %>% 
  select(-c(pubfam,mufam)) %>% 
  pivot_longer(
    -booti, names_to = 'stat', values_to = 'value'
  ) %>% 
  rbind(
    bootstat
  )

# SAVE TO OUTPUT/====
write.csv(x = bootall, file = 'output/bootall.csv', row.names = F)

