# 2024 01 Andrew simulating supplementary results

# Environment ---------------------------------------------------------------
rm(list = ls())
source('0_Settings.r')

# function for simulating factors
simulate_factors = function(lambda=4, tmin = 2, zmin=1, rho=0, N=20000){
    set.seed(1)

    # mu and noise
    mu = rnorm(N, 0, sqrt(lambda))
    ep_c = rnorm(N, 0, sqrt(rho))
    ep_t = ep_c + rnorm(N, 0, sqrt(1-rho))
    ep_z = ep_c + rnorm(N, 0, sqrt(1-rho))

    # evidence    
    t = mu + ep_t # main
    z = mu + ep_z # everything else
    
    # selection bias
    pub = (t > tmin) & (z > zmin)
    
    # assemble into data frame
    alldat = data.table(
        mu = mu, t = t, z = z, ep_t = ep_t, ep_z = ep_z, pub = pub
    )

    return(alldat)
} # end function

# function for estimating lambda
estimate_lambda = function(tbar_pub=3.5, Nest=20000, tminest=2, zminest=-Inf){
    zerome = function(ltemp){
        Et_pub = simulate_factors(lambda = ltemp, tmin=tminest, 
            zmin = zminest, N=Nest)[pub == TRUE, mean(t)]
        return(Et_pub - tbar_pub)
    }
    # find lambda
    lhat = uniroot(zerome, c(0.1, 20))$root
    return(lhat)
}

# function for estimating mu if you have both t and z
# this isn't really necessary, but it's here for completeness
estimate_mu_t_z = function(l, t, z){
    S12 = matrix(l, nrow = 1, ncol = 2)
    S22 = c(l+1, l, l, l+1) %>% matrix(nrow = 2)
    slopemat = S12 %*% solve(S22)
    muhat = slopemat %*% t(cbind(t, z))
    muhat = muhat[1, ]
    return(muhat)
}

# function for estimating mu if you only have t
estimate_mu_t = function(l, t){    
    slopemat = l/(l+1)
    muhat = slopemat * t
    return(muhat)
}

# show how zmin affects estimation of lambda --------------------------------
# E(t|pub) increases in zmin, so the lambdahat ends up too large

lambda_list = seq(1, 10, 0.2)
zmin_list = c(-Inf, 1)
par_list = expand_grid(lambda = lambda_list, zmin = zmin_list) %>% 
    mutate(Et_pub = NA_real_)
Nsim = 100*1000
Et_pub_data = 3.5

# loop over settings
for (i in 1:nrow(par_list)){    
    print(paste0(i, ' of ', nrow(par_list)))
    par0 = par_list[i,]
    temp = simulate_factors(lambda=par0$lambda, zmin=par0$zmin, N=Nsim)
    par_list$Et_pub[i] = mean(temp[pub == TRUE, t])
} # end for i in 1:nrow(par_list)

# solve for lambda for each zmin
lambda_est = par_list %>% 
    mutate(err = (Et_pub - Et_pub_data)^2) %>% 
    group_by(zmin) %>% 
    arrange(zmin,err) %>% 
    slice(1) %>% 
    select(zmin, everything())

# plot
p = par_list %>% mutate(zmin = factor(zmin, levels = zmin_list
        , labels=c('Assuming no z requirement', 'Actual: z > 1 required'))) %>% 
    ggplot(aes(x = lambda, y = Et_pub, color = zmin, linetype=zmin)) +
    geom_line(size=1.5) +
    labs(x = TeX("$\\lambda$"), y = 'Mean published t-stat') +
    chen_theme +
    scale_color_manual(values = c('gray', MATRED), name = 'Min z for publication') +
    scale_linetype_manual(values = c(1, 2), name = 'Min z for publication') +
    theme(legend.position = c(8,3)/10, legend.text = element_text(size = 16),
        legend.key.size = unit(1, 'cm')) +
    # horizontal line for data
    geom_hline(yintercept = Et_pub_data, linetype = 2, color = 'black') +
    annotate('text', x = 2, y = 3.55, label = 'Data', size = 6, color = 'black') +
    # vertical line for lambda estimates
    geom_vline(data=lambda_est[1,], aes(xintercept = lambda), color='gray') +
    annotate('text', x = 6.3, y = 2.65, label=TeX("$\\lambda_{est}$"), 
        size = 11, color = 'gray50')+
    # vertical line for actual lambda
    geom_vline(data=lambda_est[2,], aes(xintercept = lambda), color=MATRED) +
    annotate('text', x=4.5, y=2.65, label=TeX("$\\lambda_{act}$"),
        size = 11, color = MATRED)     

ggsave('output/demo-other-evidence-lambda-bias.pdf',p, width = 4, height = 2.2, 
    scale = 2.5, device = cairo_pdf)

# examine one simulation --------------------------------------------------
# check var([ep_t ep_z]) is close to [1 rho; rho 1]
# check: slope of mu on muhat should be 1 if you use the right lambda
# check: bias should be zero if zmin0 = -Inf
# check: lhat > l0 if zmin0 is not -Inf but is assumed to be

# user
l0 = 3
zmin0 = 1
tmin0 = 2
rho0 = 0.5

# simulate
alldat = simulate_factors(lambda = l0, zmin = zmin0, tmin = tmin0, rho = rho0)
pubmean = alldat %>% filter(pub == TRUE) %>% summarise(t = mean(t), z = mean(z))

# check variance
alldat %>% select(ep_t, ep_z) %>% cov()

# check shrinkage
alldat$muhat0_t = estimate_mu_t(l0, alldat$t)
lm(mu ~ muhat0_t, data = alldat)$coef

# estimate lambda
lhat = estimate_lambda(tbar_pub = pubmean$t)
lhat

# estimate shrinkage
alldat$muhat_t = estimate_mu_t(lhat, alldat$t)
lm(mu ~ muhat_t, data = alldat)$coef

# compare in terms of mean shrinkage
pubmean = alldat %>% filter(pub == TRUE) %>% 
    summarise(
        mu = mean(mu),
        t = mean(t),
        muhat_t = mean(muhat_t)
    ) %>%
    mutate(
        s0 = 1-mu/t,
        shat = 1-muhat_t/t,        
        bias = shat - s0
    )
pubmean


# examine many simulations ------------------------------------------------
lambda_list = c(seq(0.5, 3, 0.5), seq(6, 12, 4))
zmin_list = c(-Inf, 0, 1, 2)
rho_list = c(0, 0.6)
tmin0 = 2.5
par_list = expand_grid(lambda = lambda_list, zmin = zmin_list, rho = rho_list)
Nsim = 100*1000
Nest = Nsim

# loop over settings
pubmean_list = tibble()
for (i in 1:nrow(par_list)){    

    print(paste0(i, ' of ', nrow(par_list)))
    par0 = par_list[i,]
    
    # simulate data
    alldat = simulate_factors(lambda=par0$lambda, zmin=par0$zmin, rho=par0$rho, 
        N=Nsim, tmin = tmin0)
    pubdat = alldat[pub == TRUE]
    
    # estimate lambda
    lhat = estimate_lambda(tbar_pub = mean(pubdat$t), tminest=tmin0, zminest=-Inf, Nest=Nest)

    # estimate shrinkage
    pubdat$muhat = estimate_mu_t(lhat, pubdat$t)

    # summarize and save
    pubmean = pubdat %>% summarize(
        mu = mean(mu),
        t = mean(t),
        muhat = mean(muhat)
    )
    pubmean_list = bind_rows(pubmean_list, pubmean)
} # end for i in 1:nrow(par_list)

# find shrinkage
pubmean_list = pubmean_list %>% mutate(
    shrink = 1 - mu/t, shrinkhat = 1 - muhat/t, bias = shrinkhat - shrink
)


# plot shrinkage
xlimmax = 0.8
p = par_list %>% bind_cols(pubmean_list) %>% 
    mutate(zmin = as.factor(zmin), rho = as.factor(rho)) %>%
    ggplot(aes(x = shrinkhat, y = shrink
        , color = rho, shape = zmin))  +
    geom_abline(intercept = 0, slope = 1, linetype = 1, color = 'black') +
    geom_point(size = 4.5) +
    scale_color_manual(values = c(MATBLUE,MATYELLOW),
        name = TeX("$\\rho = Cor(t_i,z_i|\\mu_i)$")) +
    scale_shape_manual(values = c(16, 17, 15, 3), name = TeX("Min $z_i$ for pub")) +    
    labs(x = 'estimated shrinkage', y = 'actual shrinkage') +
    chen_theme +
    coord_cartesian(xlim = c(0, xlimmax), ylim = c(0, xlimmax)) +
    theme(
        legend.position = c(8.6,3.8)/10, 
        legend.text = element_text(size = 16),
        legend.margin = margin(0.2, 0.2, 0.2, 0.2, 'cm'),
        legend.key.size = unit(0.7, 'cm'),
        legend.title = element_text(size = 16),
        # axis lines
        panel.grid.major = element_line(color = 'gray90'),
    )     
ggsave('output/demo-other-evidence-shrink.pdf',p, width = 4, height = 2.2, 
    scale = 2.5, device = cairo_pdf)

