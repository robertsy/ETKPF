# R implementation of main routines of etkpf



#' ETKPF analysis in ensemble space
#' compute the W matrices of the analysis of ETKPF
#'
#' @param R_evc eigenvectors of S=Yb'R^-1 Yb
#' @param R_evl eigenvalues of S
#' @param C Yb'R^-1(y - Hxbbar)
#' @param n ensemble size
#' @param weps_type Weps algorithm (0: stochastic, 1: Riccati)                     
#' @param adaptive scheme to choose gamma:
#'  0: not adaptive
#'  1: in ess bounds
#'  2: min MSE
#'  3: min ES  
#'  4: min CV(MSE)
#'  5: min CV(ES)
#' @param center (0: not centered, 1: centered such that  xabar = xbbar + Xb wm and xa = xabar + Xb W)
#' @param unif used in balanced sampling
#' @param gam default gam value if not adaptive
#' @param eps  perturbations to pass in case weps=0
#' @param eps_cv  perturbations to pass for CV scheme
#' @return a list with all the analysis elements (Wmu,Weps, etc)
etkpf_es <- function(R_evc, R_evl, C, n,
                     weps_type=0,  
                     # 0: stochastic (in ensemble space)
                     # 1: Riccati solver
                     adaptive=1,   # adaptive choice of gamma:
                     # 0: not adaptive
                     # 1: in ess bounds
                     # 2: min MSE
                     # 3: min ES  
                     # 4: min CV(MSE)
                     # 5: min CV(ES)
                     center=1,   
                     # parametrize such that (necessary for cov inflation):
                     # xabar = xbbar + Xb wm
                     # xa = xabar + Xb W
                     unif=runif(1),
                     gam=0.5,    # default gam value if not adaptive
                     eps=NULL,   # perturbations to pass in case weps=0
                     eps_cv=NULL,# perturbations to pass for CV scheme
                     ...){       # additional arguments passed to adaptive_gamma
  
  ## test parameters:
  if (!(weps_type %in% 0:1)) {
    warning('invalid weps_type value')
    stop()
  }
  if (!(adaptive %in% 0:5)) {
    warning('invalid adaptive value')
    stop()
  }
  if (is.null(eps) & weps_type==0) {
    eps <- matrix(rnorm(n*n),n,n)
  }
  if (is.null(eps_cv) & adaptive %in% 4:5) {
    eps_cv <- matrix(rnorm(n*n),n,n)
  }
  if (weps_type==0){
    ## center eps:
    eps <- eps - rowMeans(eps)
  }
  # browser()
  # print(eps_cv[1:4,1])
  
  #fix cases with negative eigenvalues
  #which should theoretically not happen
  rvl <- rep(0, length(R_evl))
  pos_mask <- R_evl > 1e-5 # .Machine$double.eps
  rvl[ pos_mask ] <- R_evl [ pos_mask ]
  
  ## adaptive gamma:
  if (adaptive==0){
    alphai <-  get_alpha    (R_evc, rvl, C, n, gam)$w
    ess    <-  1/sum(alphai^2)/n
  } else{
    if (adaptive==1) ad_fit <- adaptive_gamma_ess     (R_evc, rvl, C, n,...)
    if (adaptive==2) ad_fit <- adaptive_gamma_minSCORE(R_evc, rvl, C, n, score_type='MSE',unif=unif,...)
    if (adaptive==3) ad_fit <- adaptive_gamma_minSCORE(R_evc, rvl, C, n, score_type='ES',unif=unif,...)
    if (adaptive==4) ad_fit <- adaptive_gamma_minSCORE(R_evc, rvl, C, n, score_type='MSECV',unif=unif,eps_cv=eps_cv,...)
    if (adaptive==5) ad_fit <- adaptive_gamma_minSCORE(R_evc, rvl, C, n, score_type='ESCV',unif=unif,eps_cv=eps_cv,...)
    
    alphai <- ad_fit$w
    gam    <- ad_fit$gam
    ess    <- ad_fit$ess
  }
  
  # Wmu elements:
  Wmu    <- get_Wmu    (R_evc, rvl, n, gam)
  wmubar <- get_wmubar (R_evc, rvl, C, n, gam)
  
  ## resampling:
  ind_resample   <- bal_sample (alphai, n, unif)$index
  ind_resample   <- reorder_ind(ind_resample, n)
  
  ## Weps:
  if (weps_type==0) Weps <- get_Weps_stochastic   (R_evc, rvl, n, gam, eps=eps)
  if (weps_type==1) Weps <- get_Weps_riccati      (R_evc, rvl, n, gam, Wmu, ind_resample)
  
  
  ## center such that xbbar + Xb wm = xabar
  if (center ==1){
    shift  <- rowMeans(Wmu[,ind_resample] + Weps)
    wm <- wmubar + shift
    W <- Wmu[,ind_resample] + Weps - shift
  } else{
    ## all in one such that xa = xbbar + Xb W
    W = Wmu[,ind_resample] + Weps + wmubar %*% t(rep(1,n))
    wm <- rep(0,n)
  }
  
  
  return( list(W=W, wm=wm, Wmu=Wmu, wmubar=wmubar, Weps=Weps, ind=ind_resample,
               gam=gam, ess=ess, alphai=alphai))
  
  
}





#' Wrapper around etkpf_es
#' compute the ensemble space elements needed to call etkpf_es
#'
#' @param xb the background ensemble
#' @param y the observations
#' @param H observation operator
#' @param R observation covariance
#' @param etkpf_foo function to compute the ensemble space analysis
#' @param ... additional arguments to pass to etkpf_foo
#' @return a list with xa and mod=model returned by etkpf_foo (Wmu,Weps, etc)
etkpf_analysis <- function(xb, y, H, R,
                           etkpf_foo=etkpf_es, ## can be replaced with etkpf_es_f90 
                           ...){## all parameters to pass to etkpf_foo
  
  n <- ncol(xb)

  
  ## ensemble space elements:
  xbbar <- rowMeans(xb)
  Xb <- sweep(xb, 1, xbbar)
  Yb <- H%*%xb
  Yb <- sweep(Yb, 1, apply(Yb, 1, mean))
  
  S <- t(Yb)%*%solve(R, Yb)
  ee <- eigen(S, symmetric=TRUE)
  R_evc <- ee$vectors
  R_evl <- ee$values
  ## change order as in fortran code:
  R_evc <- R_evc[, order(R_evl)]
  R_evl <- R_evl[order(R_evl)]
  Cvec <- t(Yb)%*%solve(R, y - H%*%xbbar)
  
  mod <- etkpf_foo(R_evc, R_evl, Cvec, n, ...)
  
  
  xabar <- xbbar +  Xb %*% mod$wm 
  xa <- xabar%*%t(rep(1,n)) + Xb%*% mod$W
  
  return(list( xa=xa, mod=mod ))
}













