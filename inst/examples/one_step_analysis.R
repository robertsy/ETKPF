## simple one step analysis with a Gaussian random field

set.seed(1)
q <- 100
d <- q
k <- 20



# generate data -----------------------------------------------------------
gp_data <- gp_simulate(q=q)

xb <- gp_ens0(k, gp_data)
gp_plot(gp_data, xb, 'Background')


y <- gp_data$y; H <- gp_data$H; R <- gp_data$R


myseed <- 2

## EnKF analysis:
set.seed(myseed)
eps <- matrix(rnorm(k*k), k,k)
enkf_fit <- etkpf_analysis(xb, y, H, R, 
                           weps_type = 0, adaptive = 0, gam=1, eps=eps)
gp_plot(gp_data, enkf_fit$xa, 'EnKF analysis')

## ETKF analysis:
set.seed(myseed)
etkf_fit <- etkpf_analysis(xb, y, H, R, 
                           weps_type = 1, adaptive = 0, gam=1, eps=eps, etkpf_foo = etkpf_es_f90)
gp_plot(gp_data, etkf_fit$xa, 'ETKF analysis')



## ETKPF-ess50
set.seed(myseed)
etkpf_ess50_fit <- etkpf_analysis(xb, y, H, R, 
                           weps_type = 1, adaptive = 1, ess0=0.5, ess1=0.5, etkpf_foo = etkpf_es_f90)
gp_plot(gp_data, etkpf_ess50_fit$xa, paste('ETKPF-ess50 ( gam =',etkpf_ess50_fit$mod$gam,')'))


## ETKPF-minMSE
set.seed(myseed)
etkpf_minMSE_fit <- etkpf_analysis(xb, y, H, R, 
                                  weps_type = 0, adaptive = 2, etkpf_foo = etkpf_es_f90)
gp_plot(gp_data, etkpf_minMSE_fit$xa, paste('ETKPF-minMSE ( gam =',etkpf_minMSE_fit$mod$gam,')'))


## ETKPF-minES
set.seed(myseed)
etkpf_minES_fit <- etkpf_analysis(xb, y, H, R, 
                                   weps_type = 0, adaptive = 3, etkpf_foo = etkpf_es)
gp_plot(gp_data, etkpf_minES_fit$xa, paste('ETKPF-minES ( gam =',etkpf_minES_fit$mod$gam,')'))



## ETKPF-minCV(MSE)
set.seed(myseed)
etkpf_minMSECV_fit <- etkpf_analysis(xb, y, H, R, 
                                   weps_type = 1, adaptive = 4,
                                   eps_cv=eps,
                                   etkpf_foo = etkpf_es_f90)
gp_plot(gp_data, etkpf_minMSECV_fit$xa, paste('ETKPF-minMSECV ( gam =',etkpf_minMSECV_fit$mod$gam,')'))


## ETKPF-minCV(ES)
set.seed(myseed)
etkpf_minESCV_fit <- etkpf_analysis(xb, y, H, R, 
                                    weps_type = 1, adaptive = 5, 
                                    eps_cv=eps,
                                    etkpf_foo = etkpf_es_f90)
gp_plot(gp_data, etkpf_minESCV_fit$xa, paste('ETKPF-minESCV ( gam =',etkpf_minESCV_fit$mod$gam,')'))







