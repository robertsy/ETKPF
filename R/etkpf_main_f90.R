# fortran implementation of main routine of etkpf


#' ETKPF analysis in ensemble space (in fortran90)
#' compute the W matrices of the analysis of ETKPF
#' @inheritParams etkpf_es
etkpf_es_f90 <- function(R_evc, R_evl, C, n,
                         weps_type,  # 0: stochastic (in ensemble space)
                         # 1: Riccati solver
                         adaptive,   # adaptive choice of gamma:
                         # 0: not adaptive
                         # 1: in ess bounds
                         # 2: min MSE
                         # 3: min ES
                         # 4: min CV(MSE)
                         # 5: min CV(ES)
                         center=1,   # parametrize such that (necessary for cov inflation):
                         # xabar = xbbar + Xb wm
                         # xa = xabar + Xb W
                         unif=runif(1),
                         gam=0.5,    # default gam value if not adaptive
                         eps=matrix(rnorm(n*n),n,n),   # perturbations to pass in case weps=0
                         eps_cv=matrix(rnorm(n*n),n,n),
                         ess0=0.5,ess1=0.5,imin=0,imax=40,
                         ...){
  
  ## test parameters:
  if (!(weps_type %in% 0:1)) {
    warning('invalid weps_type value')
    stop()
  }
  if (!(adaptive %in% 0:5)) {
    warning('invalid adaptive value')
    stop()
  }
  if (weps_type==0){
    ## center eps:
    eps <- eps - rowMeans(eps)
  }
  
  
  W <- matrix(0,n,n)
  wabar <- numeric(n)
  
  Wmu <- matrix(0,n,n)
  wmubar <- numeric(n)
  Weps <- matrix(0,n,n)
  alphai <- numeric(n)
  ind_resample <- as.integer(numeric(n))
  
  output <- .Fortran('etkpf',
                     R_evc = as.double(R_evc),
                     R_evl = as.double(R_evl),
                     C     = as.double(C),
                     n     = as.integer(n),
                     W     = W,
                     wabar = wabar,
                     Wmu   = Wmu,
                     wmubar= wmubar,
                     Weps  = Weps,
                     alphai= alphai,
                     ind_resample   = ind_resample,
                     unif  = as.double(unif),
                     weps_type = as.integer(weps_type),
                     adaptive  = as.integer(adaptive),
                     center= as.integer(center),
                     gam   = gam,
                     ess0  = as.double(ess0),
                     ess1  = as.double(ess1),
                     imin  = as.double(imin),
                     imax  = as.double(imax),
                     eps   = as.double(eps),
                     eps_cv= as.double(eps_cv))
  
  if (center ==1){
    W <- output$W
    wm <- output$wabar
  } else{
    ## all in one such that xa = xbbar + Xb W
    W <- output$W + output$wabar %*% t(rep(1,n))
    wm <- rep(0,n)
  }
  
  
  return( list(W=W, wm=wm, Wmu=output$Wmu, wmubar=output$wmubar, Weps=output$Weps, ind=output$ind_resample,
               gam=output$gam, ess=1/sum(output$alphai^2)/n , alphai=output$alphai))
  
  
}

