# check normality

# ====

# download returns and add samptype
target_dribble = pathRelease %>% drive_ls() %>%
  filter(name=='Portfolios') %>% drive_ls() %>%
  filter(name=='Full Sets OP') %>% drive_ls() %>%
  filter(name=='PredictorLSretWide.csv')

drive_download(target_dribble[1,], path = 'intermediate/deleteme2.csv', overwrite = T)

ret = fread('intermediate/deleteme2.csv') %>% 
  pivot_longer(
    -date, names_to = 'signalname', values_to = 'ret'
  ) %>% 
  filter(!is.na(ret))


# ====


ndate = 300
nboot = 5000

r = sample(ret$ret, size = 300*nboot, replace = T) %>% 
  matrix(nrow = ndate, ncol = nboot)

rbar = apply(r, 2, mean)
vol  = apply(r, 2, sd)
tstat = rbar/vol*sqrt(ndate) 

tcenter = tstat - mean(tstat)

sd(tcenter)

qqnorm(tcenter)