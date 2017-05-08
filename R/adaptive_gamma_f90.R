## fortran implementation of adaptive gamma functions
## see adaptive_gamma.R for the documentation


adaptive_gamma_ess_f90 <- function(R_evc, R_evl, C, n, ess0, ess1, imin=1, imax=15, delta_gam=1/imax){
  alphai <- numeric(n)
  gam <- 0
  output <- .Fortran('adaptive_gamma_ess',
                     R_evc = as.double(R_evc),
                     R_evl = as.double(R_evl),
                     C     = as.double(C),
                     n     = as.integer(n),
                     alphai= alphai,
                     gam   = gam,
                     ess0  = as.double(ess0),
                     ess1  = as.double(ess1),
                     imin_p= as.double(imin),
                     imax_p= as.double(imax))
  
  
  return( list(gam=output$gam, w=output$alphai, ess=1/sum(output$alphai^2)/n) )
}



adaptive_gamma_minMSE_f90 <- function(R_evc, R_evl, C, n,
                                        imin=0, imax=40, delta_gam=1/imax,
                                        unif=runif(1)){
  alphai <- numeric(n)
  gam <- 0
  output <- .Fortran('adaptive_gamma_minMSE',
                     R_evc = as.double(R_evc),
                     rvl   = as.double(R_evl),
                     C     = as.double(C),
                     n     = as.integer(n),
                     alphai= alphai,
                     gam   = gam,
                     imin  = as.double(imin),
                     imax  = as.double(imax),
                     unif  = as.double(unif))
  
  
  return( list(gam=output$gam, w=output$alphai, ess=1/sum(output$alphai^2)/n) )
}


adaptive_gamma_minES_f90 <- function(R_evc, R_evl, C, n,
                                      imin=0, imax=40, delta_gam=1/imax,
                                      unif=runif(1)){
  alphai <- numeric(n)
  gam <- 0
  output <- .Fortran('adaptive_gamma_minES',
                     R_evc = as.double(R_evc),
                     rvl   = as.double(R_evl),
                     C     = as.double(C),
                     n     = as.integer(n),
                     alphai= alphai,
                     gam   = gam,
                     imin  = as.double(imin),
                     imax  = as.double(imax),
                     unif  = as.double(unif))
  
  
  return( list(gam=output$gam, w=output$alphai, ess=1/sum(output$alphai^2)/n) )
}



adaptive_gamma_minMSECV_f90 <- function(R_evc, R_evl, C, n,
                                         imin=0, imax=40, delta_gam=1/imax,
                                         unif=runif(1),
                                         eps=matrix(rnorm(n*n),n,n),
                                         ...){
  alphai <- numeric(n)
  gam <- 0
  output <- .Fortran('adaptive_gamma_minMSECV',
                     R_evc = as.double(R_evc),
                     rvl   = as.double(R_evl),
                     C     = as.double(C),
                     n     = as.integer(n),
                     alphai= alphai,
                     gam   = gam,
                     imin  = as.double(imin),
                     imax  = as.double(imax),
                     unif  = as.double(unif),
                     eps   = as.double(eps))

  return( list(gam=output$gam, w=output$alphai, ess=1/sum(output$alphai^2)/n) )
}




adaptive_gamma_minESCV_f90 <- function(R_evc, R_evl, C, n,
                                        imin=0, imax=40, delta_gam=1/imax,
                                        unif=runif(1),
                                        eps=matrix(rnorm(n*n),n,n),
                                        ...){
  alphai <- numeric(n)
  gam <- 0
  output <- .Fortran('adaptive_gamma_minESCV',
                     R_evc = as.double(R_evc),
                     rvl   = as.double(R_evl),
                     C     = as.double(C),
                     n     = as.integer(n),
                     alphai= alphai,
                     gam   = gam,
                     imin  = as.double(imin),
                     imax  = as.double(imax),
                     unif  = as.double(unif),
                     eps   = as.double(eps))
  
  
  return( list(gam=output$gam, w=output$alphai, ess=1/sum(output$alphai^2)/n) )
}






## needed for ES score
f_crps_wgm_f90 <- function(mu, sig2, obs, w){
  
  crps = 0
  n <- length(mu)
  if (length(sig2) > 1) {
    warning('f90 function does not support different sig2')
    sig2 <- sig2[1]
  }
  
  output <- .Fortran('f_crps',
                     mu    = as.double(mu),
                     sig2  = as.double(sig2),
                     w     = as.double(w),
                     obs   = as.double(obs),
                     n     = as.integer(n),
                     crps  = crps)
  
  return( output$crps)
}


