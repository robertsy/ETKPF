## R implementation of adaptive gamma functions


#' ESS target adaptive gamma criterion
#' Binary search such that ess0 <= ESS <= ess1
#'
#' @param ess0 ESS lower bound
#' @param ess1 ESS upper bound
#' @param imax, imin: gammas considered go from imin/(imax-imin+1) to imax/(imax-imin+1),
#' @param delta_gam gammas increments, can be changed
#' @inheritParams etkpf_es
#' @return list(gamma, w, ess) of chosen model
adaptive_gamma_ess <- function(R_evc, R_evl, C, n,
                               ess0=0.5, ess1=0.5,
                               imin=0, imax=40, delta_gam=1/imax){
  #find minimum gamma such that ESS is between ess0 and ess1
  best_model <- NA
  best_ess <- 2 # best ess (max=1)
  while( imax >= imin ){
    imid <- ceiling( imin + (imax-imin)/2)
    gam <- imid * delta_gam
    weights <- get_alpha(R_evc, R_evl, C, n, gam)
    ess <- weights$ess
    
    ## if ess is smaller than necessary, increase gam:
    if(ess < ess0){ #increase gam, no change in best_ess
      imin <- imid + 1
    } else{
      ## if ess.mid is bigger than best ess so far, decrease gam:
      if(ess > best_ess){
        imax <- imid - 1
        ## otherwise, update best_ess and best_model
      } else{ #update ess.0
        best_ess <- ess
        best_model <- list(gam=gam, w=weights$w, ess=ess)
        ## if current ess > e.1, decrease gam
        if (ess > ess1){
          imax <- imid - 1
        } else{ ## set imax such that it breaks out of the loop:
          imax <- -1
        }
      }
    }
  }
  
  return(best_model)
}



#' Score-based adaptive gamma criterion
#' Grid search to minimize criterion
#'
#' @param score_type describe which score to use (MSE, ES, MSECV, ESCV)
#' @param imax, imin: gammas considered go from imin/(imax-imin+1) to imax/(imax-imin+1),
#' @param delta_gam gammas increments, can be changed
#' @inheritParams etkpf_es
#' @return list(gam_best, w_best, ess_best, gammas, errors), where the last two are used to debug/plot
adaptive_gamma_minSCORE <- function(R_evc, R_evl, C, n,
                                    imin=0, imax=40, 
                                    delta_gam=1/imax,
                                    unif=runif(1), 
                                    score_type='MSE',
                                    # MSE, ES, MSECV, ESCV
                                    eps_cv=matrix(rnorm(n*n), n,n), ## for fake data in CV
                                    ...){
  ## test parameters:
  if (!(score_type %in% c('MSE', 'ES','MSECV', 'ESCV'))) {
    warning('invalid score_type value')
    stop()
  }
  gammas <- seq(imin, imax) * delta_gam
  ngams <- length(gammas)
  es_errors <-vector('numeric', ngams)
  for (i in 1:ngams){
    gam <- gammas[i]
    
    Wmu <- get_Wmu( R_evc, R_evl, n, gam )
    wmubar <- get_wmubar( R_evc, R_evl, C, n, gam )
    alphai <- get_alpha( R_evc, R_evl, C, n, gam )$w
    ind_resample <- bal_sample(alphai, n, unif)$index
    
    Wa <- Wmu[, ind_resample]
    shift <- apply(Wa, 1, mean)
    wm <- wmubar + shift
    
    S <- R_evc%*%diag(R_evl)%*%t(R_evc)
    
    ## component covariance:
    del_pag <- ( gam*R_evl )/ ( gam*R_evl^2 + 2*(n-1)*gam*R_evl + (n-1)^2 )
    

    if (score_type == 'MSE') {
    ## compute Err(gam) = wm'S wm - 2 wm'C:
    es_errors[i] <- t(wm) %*% S %*% wm -
      2*t(wm) %*% C
    }
    if (score_type == 'ES') {
      rvl_inv <- 1/R_evl
      rvl_inv[ rvl_inv >= 1e5 ] <- 0
      es_errors[i] <-  energy_score_es(Wmu, ind_resample, wmubar, R_evc, C, S, R_evl, rvl_inv, del_pag, n)
    }
    if (score_type == 'MSECV') {
      es_errors[i] <- cv_score(gam, R_evc, R_evl, S, unif,  eps_cv, 'MSE', n)
    }
    if (score_type == 'ESCV') {
      es_errors[i] <- cv_score(gam, R_evc, R_evl, S, unif,  eps_cv, 'ES', n)
    }
    
  }

  ind_opt <- min(which.min(es_errors))
  best_gam <- gammas[ind_opt]
  alphai <- get_alpha( R_evc, R_evl, C, n, best_gam )$w
  ess_opt <- 1/sum(alphai^2)/n
  best_model <- list(gam=best_gam, w=alphai, ess=ess_opt, gammas=gammas, errors=es_errors)
  

  return(best_model)
}










#' Energy score of yhat predicted by xahat, where hat is the projection onto orthogonal basis of ensemble space:
#'
#' @param Wmu as returned by get_Wmu
#' @param ind vector of resampling indices
#' @param wmubar as returned by get_wmubar
#' @param S S=Yb' R^-1 Yb'
#' @param rvl eigenvalues of S, already floored at 0
#' @param rvl_inv inverse of rvl such that 1/0=0
#' @param del_pag eigenvalues of Pag: Pag = R_evc diag(del_pag) R_evc'
#' @inheritParams etkpf_es
#' @return estimated energy score
energy_score_es <- function(Wmu, ind, wmubar, R_evc, Cvec, S, rvl, rvl_inv, del_pag, n){
  ycoord <- diag(sqrt(rvl_inv) ) %*% t(R_evc) %*% Cvec
  ## Wa in orthogonal basis:
  Worth <-  diag(sqrt(rvl)) %*% t(R_evc) %*% ( Wmu[,ind] + matrix(wmubar, n,n) )
  ## component-wise CRPS in orthogonal basis:
  all_crps_ens <- foreach(ii = 1:n, .combine = 'c')%do%
    f_crps_wgm( Worth[ii,], rep(rvl[ii]*del_pag[ii] + 1, n), ycoord[ii],  rep(1/n, n))
  
  mean(all_crps_ens)
}






#' CV score of ES, MSE or EMSE (in ensemble space)
#'
#' @param score_type MSE or ES
#' @inheritParams etkpf_es
#' @inheritParams energy_score_es
#' @return estimated CV score 
cv_score <- function(gam,
                     R_evc, rvl,
                     Sfull, 
                     unif,
                     eps_cv, ## perturbations for fake data
                     score_type, ## MSE or ES
                     n
){
  ## test parameters:
  if (!(score_type %in% c('MSE', 'ES'))) {
    warning('invalid score_type value in cv_score')
    stop()
  }
  
  ## compute eta_j :
  rvl_inv <- 1/rvl
  rvl_inv[ rvl_inv >= 1e5 ] <- 0
  Eta <- R_evc %*% diag( sqrt(rvl_inv)) %*% t(R_evc) %*% eps_cv
  
  onesj <- rep(1, n-1)
  
  cv_errors  <- vector('numeric', n)

  for (j in 1:n){
    ## S matrix:
    A_j <- diag(n)
    A_j[j,] <- 1/(n-1)
    A_j <- A_j[,-j]
    Stemp <- Sfull %*% A_j
    S_j <- t(A_j) %*% Sfull %*% A_j
    
    ## eigen S:
    ee_j <- eigen(S_j)
    R_evc_j <- ee_j$vectors
    R_evl_j <- ee_j$values
    R_evc_j <- R_evc_j[, order(R_evl_j)]
    R_evl_j <- R_evl_j[order(R_evl_j)]
    
    rvl_j <- rep(0, length(R_evl_j))
    pos_mask_j <- R_evl_j >  1e-5
    rvl_j[ pos_mask_j ] <- R_evl_j [ pos_mask_j ]
    
    
    rvl_inv_j <- 1/rvl_j
    rvl_inv_j[ rvl_inv_j >= 1e5 ] <- 0
    
    ## Cvec:
    ej <- rep(0, n)
    ej[j] <- 1
    C_j <- t(A_j) %*% Sfull %*% (  n/(n-1) * ej + Eta[,j] )
    
    ## fit:
    ## center=0 because everything should be with respect to the background mean
    ## eps=0 such that the W shows the mixture means
    fit_j <-  etkpf_es(R_evc_j, rvl_j, C_j, n-1, weps_type = 0, adaptive = 0, gam=gam,
                           center=0, unif = unif, eps=matrix(0,n-1,n-1)) ## if passing eps=0, the W is without random perturbations
    fit <- list(mod=fit_j) ## to get same structure
    
    
    del_pag_j <- ( gam*rvl_j )/ ( gam*rvl_j^2 + 2*(n-2)*gam*rvl_j + (n-2)^2 )
    

    if (score_type=='MSE'){
      ## MSE(yhat, xabarhat)
      wm <- fit$mod$wmubar + apply(fit$mod$Wmu[,fit$mod$ind], 1, mean)
      cv_errors[j] <- t(wm) %*% S_j %*% wm -
        2*t(wm) %*% C_j
    }
    if (score_type=='ES'){
      ## Energy score
      cv_errors[j] <-  energy_score_es(fit$mod$Wmu, fit$mod$ind, fit$mod$wmubar,
                                       R_evc_j, C_j, S_j, rvl_j, rvl_inv_j, del_pag_j, n-1)
      
      # if (gam==0.5 & j==2) print( cv_errors[j])
      
    }

  }
  
  return( mean(cv_errors) )
  
  
}




#' Compute CRPS score for a predictive ensemble (O(n*log(n)) Variante)
#' Author: Marco
#'
#' @param smple predictive ensemble (a vector)
#' @param obs observed value (a scalar)
#' @return the crps score
#' @examples
#' set.seed(1)
#' truth <- 0
#' n <- 100
#' f_crps( rnorm(n, truth, 1), truth)
#' f_crps( rnorm(n, truth, 0.1), truth)
#' f_crps( rnorm(n, truth+2, 1), truth)
#' ens.means <- rnorm(n, truth, 1)
#' plot(density(ens.means)); abline(v=truth)
#' @family forecasting_score
f_crps <- function(smple, obs)
{
  N <- length(smple)
  smple <- sort(smple) ## sort ensemble
  f <- 1/N*(1:(N-1))
  return((max(smple[1],obs)-obs) + (max(smple[N],obs)-smple[N])
         + sum(f^2*diff(smple)-(2*f-1)*diff(pmax(smple,obs))))
}


f_crps_wgm <- function(mu, sig2, obs, w)
{
  ## Purpose: Compute CRP score for a Gaussian mixture
  ## where sig2 is the same but mixtures have weights w
  ## Arguments: mu (mixture means), sig2 (mixture variance),
  ##                      obs (observations), w (mixture component probability)
  ## Author: Christian Kerkhoff, Date: 11 Nov 2013, 13:36
  if (length(sig2)==1) sig2 <- rep(sig2, length(mu))
  
  A <- function(mu, sig2){
    ifelse(  sig2 <= 1e-10, abs(mu),
             2*sqrt(sig2)*dnorm(mu/sqrt(sig2)) + mu*(2*pnorm(mu/sqrt(sig2))-1) )
  }
  
  return( weighted.mean( A(obs-mu, sig2), w) -
            0.5* sum(outer(w,w, FUN="*") * A( outer(as.numeric(mu),as.numeric(mu),FUN="-"),
                                              outer(as.numeric(sig2),as.numeric(sig2),FUN="+") ) ) )
}


