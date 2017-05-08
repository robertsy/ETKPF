## fortran implementation of basic etkpf functions
## see etkpf_util_R.R for the documentation

# Wmu part ---------------------------------------------------------------


get_Wmu_f90 <- function(R_evc, R_evl, n, gam){

  wmu <- matrix(0,n,n)
  output <- .Fortran('get_Wmu',
                     R_evc = as.double(R_evc),
                     R_evl = as.double(R_evl),
                     n     = as.integer(n),
                     Wmu   = wmu,                 ## important to not pass as.double, no idea why...
                     gam   = as.double(gam))

  return(output$Wmu)
}




get_wmubar_f90 <- function(R_evc, R_evl, C, n, gam){
  wmubar <- numeric(n)
  output <- .Fortran('get_wmubar',
                     R_evc = as.double(R_evc),
                     R_evl = as.double(R_evl),
                     C     = as.double(C),
                     n     = as.integer(n),
                     wmubar= wmubar,
                     gam   = as.double(gam))
  return(output$wmubar)
}



# Walpha part --------------------------------------------------------------


get_alpha_f90 <- function(R_evc, R_evl, C, n, gam){
  alphai <- numeric(n)
  ess <- 0
  output <- .Fortran('get_alpha',
                     R_evc = as.double(R_evc),
                     R_evl = as.double(R_evl),
                     C     = as.double(C),
                     n     = as.integer(n),
                     alphai= alphai,
                     gam   = as.double(gam),
                     ess   = as.double(ess))
  return(list(w=output$alphai,  ess=output$ess))
}



bal_sample_f90 <- function(w, R=length(w), unif=runif(1)){
  n <- length(w)
  ind_resample <- as.integer(numeric(n))
  output <- .Fortran('bal_sample',
                     pai  = as.double(w),
                     n    = as.integer(n),
                     ind  = ind_resample,
                     unif = as.double(unif))
  ind_resample <- output$ind
  Ni <- numeric(n)
  Ni[sort(unique(ind_resample))] <- table(ind_resample)
  return(list(N=Ni,  index=ind_resample) )
}



# Rearrange indices such that they match 1:n
reorder_ind_f90 <- function(ind_resample, n){
  ind_reordered <- as.integer(numeric(n))
  output <- .Fortran('reorder_ind',
                     ind_in = as.integer(ind_resample),
                     out    = ind_reordered,
                     n      = as.integer(n))
  return(output$out)
}




# Weps part: --------------------------------------------------------------


get_Weps_stochastic_f90 <- function( R_evc, R_evl, n, gam, eps ){
  ## where eps is a nxn matrix of iid N(0,1)
  ## center eps (done in fortran code)
  # epsbar <- apply(eps, 1, mean)
  # eps <- eps - epsbar %*% t(rep(1,n))

  Weps <- matrix(0,n,n)
  output <- .Fortran('get_Weps_stoch',
                     R_evc = as.double(R_evc),
                     R_evl = as.double(R_evl),
                     n     = as.integer(n),
                     Weps  = Weps,                 ## important to not pass as.double, no idea why...
                     gam   = as.double(gam),
                     eps   = as.double(eps))

  return(output$Weps)
}




get_Weps_riccati_f90 <- function( R_evc, R_evl, n, gam, Wmu, ind_resample,
                                  tol=10^-9, maxit=20){
  Weps <- matrix(0,n,n)
  output <- .Fortran('get_Weps_riccati',
                     R_evc = as.double(R_evc),
                     R_evl = as.double(R_evl),
                     n     = as.integer(n),
                     Weps  = Weps,                 ## important to not pass as.double, no idea why...
                     gam   = as.double(gam),
                     Wmu   = as.double(Wmu),
                     ind_resample = as.integer(ind_resample))

  return(output$Weps)
}






lyapunov_f90 <- function(A,C,n){
  ## solves for AX + XA' = C
  X <- matrix(0,n,n)
  INFO <- 0
  output <- .Fortran('lyap',
                     X     = X,
                     A     = as.double(A),
                     C     = as.double(C),
                     n     = as.integer(n),
                     INFO  = as.integer(INFO))

  return(output$X)
}
