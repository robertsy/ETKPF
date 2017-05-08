context("Test equivalence of f90 and R implementations")


set.seed(1)
q <- 100
d <- q
k <- 10


l <- 10
taper <- ring_GC(q, l)

# generate data -----------------------------------------------------------
gp_data <- gp_simulate(q=q)

xb <- gp_ens0(k, gp_data)
gp_plot(gp_data, xb, 'Background')

H <- gp_data$H
R <- gp_data$R
y <- gp_data$y


# Ensemble space elements -------------------------------------------------

xbbar <- rowMeans(xb)
Xb <- sweep(xb, 1, xbbar)
Yb <- H%*%xb
Yb <- sweep(Yb, 1, apply(Yb, 1, mean))

ytRy <- t(Yb)%*%solve(R, Yb)
ee <- eigen(ytRy)
R_evc <- ee$vectors
R_evl <- ee$values
## change order as in fortran code:
R_evc <- R_evc[, order(R_evl)]
R_evl <- R_evl[order(R_evl)]
Cvec <- t(Yb)%*%solve(R, y - H%*%xbbar)

#fix cases with negative eigenvalues
#which should theoretically not happen
rvl <- rep(0, length(R_evl))
pos_mask <- R_evl > 1e-5#.Machine$double.eps #1e-5
rvl[ pos_mask ] <- R_evl [ pos_mask ]



# Wmu matrices computation --------------------------------------------------
gam <- 0.7
## to test routines internally:
n <- k; R_evl <- rvl; C <- Cvec

Wmu_R   <- get_Wmu    (R_evc, rvl, k, gam)
Wmu_f90 <- get_Wmu_f90(R_evc, rvl, k, gam)
test_that("Wmu is the same",expect_equal( Wmu_R, Wmu_f90 ))

wmubar_R   <- get_wmubar    (R_evc, rvl,Cvec, k, gam)
wmubar_f90 <- get_wmubar_f90(R_evc, rvl,Cvec, k, gam)
test_that("wmubar is the same",expect_equal( wmubar_R, wmubar_f90 ))




# Walpha computation ------------------------------------------------------

alphai_R   <- get_alpha    (R_evc, rvl,Cvec, k, gam)
alphai_f90 <- get_alpha_f90(R_evc, rvl,Cvec, k, gam)
test_that("get_alpha is the same",expect_equal( alphai_R, alphai_f90 ))

unif <- runif(1)
ind_resample_R   <- bal_sample    (alphai_R$w, k, unif)
ind_resample_f90 <- bal_sample_f90(alphai_R$w, k, unif)
test_that("bal_sample is the same",expect_equal( ind_resample_R, ind_resample_f90 ))

ind_in <- sample(1:k, k, replace = T)
ind_resample_R   <- reorder_ind    (ind_in,k)
ind_resample_f90 <- reorder_ind_f90(ind_in,k)
test_that("reorder is the same",expect_equal( ind_resample_R, ind_resample_f90 ))



# Adaptive gamma ----------------------------------------------------------

gam_an_R   <- adaptive_gamma_ess    (R_evc, rvl, Cvec, k, 0.1, 0.3, imin=0, imax=40, delta_gam=1/40)
gam_an_f90 <- adaptive_gamma_ess_f90(R_evc, rvl, Cvec, k, 0.1, 0.3, imin=0, imax=40, delta_gam=1/40)
test_that("adaptive_gamma_ess is the same",expect_equal( gam_an_R, gam_an_f90 ))


gam_minerror_R   <- adaptive_gamma_minSCORE    (R_evc, rvl, Cvec, k,  imin=0, imax=40, delta_gam=1/40, unif, 
                                                score_type='MSE')
gam_minerror_f90 <- adaptive_gamma_minMSE_f90(R_evc, rvl, Cvec, k,  imin=0, imax=40, delta_gam=1/40, unif)
test_that("adaptive_gamma_minMSE is the same",expect_equal( gam_minerror_R[1:3], gam_minerror_f90 ))


gam_minerror_R   <- adaptive_gamma_minSCORE    (R_evc, rvl, Cvec, k,  imin=0, imax=40, delta_gam=1/40, unif, 
                                                score_type='ES')
gam_minerror_f90 <- adaptive_gamma_minES_f90(R_evc, rvl, Cvec, k,  imin=0, imax=40, delta_gam=1/40, unif)
test_that("adaptive_gamma_minES is the same",expect_equal( gam_minerror_R[1:3], gam_minerror_f90 ))


EE <- matrix( rnorm(k*k), k, k)
gam_minerror_R   <- adaptive_gamma_minSCORE    (R_evc, rvl, Cvec, k,  imin=0, imax=40, delta_gam=1/40, 
                                                unif=unif, score_type='MSECV', eps=EE)
# print(gam_minerror_R$gam)
# plot(gam_minerror_R$gammas, gam_minerror_R$errors)
gam_minerror_f90 <- adaptive_gamma_minMSECV_f90(R_evc, rvl, Cvec, k,  imin=0, imax=40, delta_gam=1/40, 
                                                unif=unif, eps=EE)
# print(gam_minerror_f90$gam)

test_that("adaptive_gamma_minMSECV is the same",expect_equal( gam_minerror_R[1:3], gam_minerror_f90 ))



EE <- matrix( rnorm(k*k), k, k)
gam_minerror_R   <- adaptive_gamma_minSCORE    (R_evc, rvl, Cvec, k,  imin=0, imax=40, delta_gam=1/40, unif, 
                                                score_type='ESCV', eps=EE)
# print(gam_minerror_R$gam)
gam_minerror_f90 <- adaptive_gamma_minESCV_f90(R_evc, rvl, Cvec, k,  imin=0, imax=40, delta_gam=1/40, unif, eps=EE)
# print(gam_minerror_f90$gam)

test_that("adaptive_gamma_minESCV is the same",expect_equal( gam_minerror_R[1:3], gam_minerror_f90 ))


# # gam_mincverror_R   <- adaptive_gamma_mincverror    (R_evc, rvl, Cvec, k,  imin=0, imax=40, delta_gam=1/40, unif)
# # gam_mincverror_f90 <- adaptive_gamma_mincverror_f90(R_evc, rvl, Cvec, k,  imin=0, imax=40, delta_gam=1/40, unif)
# # test_that("adaptive_gamma_mincverror is the same",expect_equal( gam_mincverror_R, gam_mincverror_f90 ))
# 
# 
# 
# eps <- matrix(rnorm(k*k),k,k)
# unif <- runif(1)
# gam_mincvfake_R   <- adaptive_gamma_mincvfake    (R_evc, rvl, Cvec, k,  imin=0, imax=40, delta_gam=1/40, unif, eps=eps)
# gam_mincverror_f90 <- adaptive_gamma_mincvfake_f90(R_evc, rvl, Cvec, k,  imin=0, imax=40, delta_gam=1/40, unif, eps=eps)
# test_that("adaptive_gamma_mincvfake is the same",expect_equal( gam_mincvfake_R, gam_mincverror_f90 ))
# 
# rvl0 <- rep(0,k)
# gam_mincvfake_R   <- adaptive_gamma_mincvfake    (R_evc, rvl0, Cvec, k,  imin=0, imax=40, delta_gam=1/40, unif, eps=eps)
# gam_mincverror_f90 <- adaptive_gamma_mincvfake_f90(R_evc, rvl0, Cvec, k,  imin=0, imax=40, delta_gam=1/40, unif, eps=eps)
# test_that("adaptive_gamma_mincvfake is the same with rvl=0",expect_equal( gam_mincvfake_R, gam_mincverror_f90 ))




# CRPS --------------------------------------------------------------------
nn <- 100
mu <- rnorm(nn)
obs <- 0
sig2 <- 1#rep(1, nn)
w <- abs(rnorm(nn))
w <- w/sum(w)
crps_R <- f_crps_wgm(mu, sig2, obs, w)
crps_f90 <- f_crps_wgm_f90(mu, sig2, obs, w)
all.equal(crps_R, crps_f90, tolerance = .0001)
test_that("f_crps_wgm is the same",expect_equal( crps_R, crps_f90, tolerance = .0001))






# Weps computation --------------------------------------------------------
eps <- matrix( rnorm(k*k), k, k)

Weps_R   <- get_Weps_stochastic    (R_evc, rvl, k, gam, eps)
Weps_f90 <- get_Weps_stochastic_f90(R_evc, rvl, k, gam, eps)
all.equal(Weps_R, Weps_f90)

test_that("get_Weps_stochastic is the same",
          expect_equal( Weps_R, Weps_f90 ))


## Riccati:
unif <- runif(1)
wi <- alphai_R$w
ind_resample_R   <- bal_sample    (wi, k, unif)
ind_resample   <- ind_resample_R$index
ind_resample   <- reorder_ind    (ind_resample,k)

Wmu   <- get_Wmu    (R_evc, rvl, k, gam)

A <- Wmu[,ind_resample] - rowMeans(Wmu[,ind_resample])
A <- t(A)
## initial condition:
beta <- 1.1 * sqrt(sum(A^2))
C = A + beta*diag(n)
E = 2*diag(n)
X1_R   <-  lyapunov    (C,t(C), E, n)
X1_f90 <-  lyapunov_f90(C,E, n)
all.equal(X1_R, X1_f90)
test_that("Riccati initial conditions are the same",
          expect_equal( X1_R, X1_f90 ))


Weps_ricc_R   <- get_Weps_riccati    (R_evc, rvl, k, gam, Wmu, ind_resample)
Weps_ricc_f90 <- get_Weps_riccati_f90(R_evc, rvl, k, gam, Wmu, ind_resample)
test_that("get_Weps_riccati is the same",
          expect_equal( Weps_ricc_R, Weps_ricc_f90 ))


## special case where one particle resampled:
ind_one <- rep(1,k)
Weps_ricc_R   <- get_Weps_riccati    (R_evc, rvl, k, gam, Wmu, ind_one)
Weps_ricc_f90 <- get_Weps_riccati_f90(R_evc, rvl, k, gam, Wmu, ind_one)
all.equal(Weps_ricc_R, Weps_ricc_f90)
test_that("get_Weps_riccati with ind=1 is the same",
          expect_equal( Weps_ricc_R, Weps_ricc_f90 ))

# 
# ## Riccati with exact covariance (not stable on COSMO)
# Weps_ricce_R   <- get_Weps_riccati_exact    (R_evc, rvl, k, gam, Wmu, ind_resample, wi)
# Weps_ricce_f90 <- get_Weps_riccati_exact_f90(R_evc, rvl, k, gam, Wmu, ind_resample, wi)
# test_that("get_Weps_riccati_exact with ind=1 is the same",
#           expect_equal( Weps_ricce_R, Weps_ricce_f90 ))
# 
# 


# ETKPF analysis ----------------------------------------------------------
unif <- runif(1)
EE <- matrix( rnorm(k*k), k, k)
# epsbar <- apply(EE, 1, mean)
# EE <- EE - epsbar %*% t(rep(1,n))

# EE <- matrix(0, k,k)

ad_vec <- 0:5
weps_vec <- 0:1
cc <- 1
for (ad_par in ad_vec){
  for (weps_par in weps_vec){

    
    EE <- matrix( rnorm(k*k), k, k)
    EEcv  <- matrix( rnorm(k*k), k, k)
    
    fit_R <- etkpf_es(R_evc, rvl, Cvec, k,
                      adaptive = ad_par,
                      weps_type = weps_par, 
                      eps_cv=EEcv,eps=EE,
                      unif = unif, center=cc)

    fit_f90 <- etkpf_es_f90(R_evc, rvl, Cvec, k,
                            adaptive = ad_par,
                            weps_type = weps_par, 
                            eps=EE, eps_cv = EEcv,
                            unif = unif, center=cc)

# rvl[rvl <= 1e-4] <- 0
    Weps_R   <- get_Weps_stochastic    (R_evc, rvl, k, gam=0.5, EE)
    Weps_f90 <- get_Weps_stochastic_f90(R_evc, rvl, k, gam=0.5, EE)
    all.equal(Weps_R, Weps_f90)

    # all.equal(fit_R$Weps, Weps_R)
    # all.equal(fit_f90$Weps, Weps_f90)

    all.equal(fit_R, fit_f90)

    test_that(paste("etkpf_es is the same with adaptive", ad_par, "weps_type", weps_par) ,expect_equal( fit_R, fit_f90 ))

  }
}



