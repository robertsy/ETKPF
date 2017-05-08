## R implementation of basic etkpf functions



# Wmu part ---------------------------------------------------------------



#' Compute Wmu
#' 
#' @description 
#' Wmu is defined such that mu = xbbar + Xb wmubar + Xb Wmu
#'
#' @param gam fixed gamma value
#' @inheritParams etkpf_es
#' @return Wmu
get_Wmu <- function(R_evc, R_evl, n, gam){
  del_mu <- ( (n-1)*gam*R_evl + (n-1)^2 ) / ( gam*R_evl^2 + 2*(n-1)*gam*R_evl + (n-1)^2 )
  Wmu <- R_evc %*% diag(del_mu) %*% t(R_evc)
  return(Wmu)
}



#' Compute wmubar
#' 
#' @description 
#' wmubar is defined such that mu = xbbar + Xb wmubar + Xb Wmu
#'
#' @param gam fixed gamma value
#' @inheritParams etkpf_es
#' @return wmubar
get_wmubar <- function(R_evc, R_evl, C, n, gam){
  del_wmubar <- ( gam + ( (1-gam)*(n-1)*gam*R_evl /
                            (gam*R_evl^2 + 2*(n-1)*gam*R_evl + (n-1)^2)  ) ) / ( (n-1) + gam*R_evl )
  mat1 <- R_evc %*% diag(del_wmubar) %*% t(R_evc)
  wmubar <- mat1 %*% C
  return( as.numeric( wmubar ) )
}




# Walpha part --------------------------------------------------------------


#' Compute alphai weights
#'
#' @param gam fixed gamma value
#' @inheritParams etkpf_es
#' @return alphai
get_alpha <- function(R_evc, R_evl, C, n, gam){
  del_alpha <- ( (n-1)^2*(1-gam) ) / ( gam*R_evl^2 + 2*(n-1)*gam*R_evl + (n-1)^2 )
  mat1 <- R_evc %*% diag( del_alpha * R_evl) %*% t(R_evc)
  mat2 <- R_evc %*% diag( del_alpha )        %*% t(R_evc)
  alphai <- diag(mat1) - 2*mat2%*%C

  # to avoid over/under-flows:
  alphai <- exp( -.5 *(alphai- min(alphai)) )

  if (sum(alphai) < 1.e-12) { alphai <- rep(1, n) }
  alphai <- alphai/sum(alphai)


  return(list (w=as.numeric(alphai), ess=1/sum(alphai^2)/n) )
}



#' Purpose: Balanced sampling
#' 
#' @description 
#' Multiplicities N_1, ..., N_n, i.e. element j is sampled N_j times
#' index = numbers of sampled elements
#' Author: Hans-Rudolf Kuensch, Date: 21 Apr 2010.
#' Modified by Sylvain Robert (added unif as argument and some names)
#'
#' @param w = probabilities (must sum to one !)
#' @param R = sample size
#' @param unif = uniform (pass if deterministic behaviour wished)
#' @return list(N, index)
bal_sample <- function(w, R=length(w), unif=runif(1)){
  n <- length(w)
  M <- floor(R*cumsum(w) + unif)
  N <- M - c(0,M[-n])
  index <- rep(1:n,N)
  list(N=N,index=index)
}




#' Rearrange indices such that they match 1:n
#'
#' @param ind_resample vector of resampling indices
#' @param n ensemble size
#' @return permuted indices such that they match 1:n as much as possible
reorder_ind_fast <- function(ind_resample, n){
  ## permute indices such that
  ## there are a minimum of mismatches with
  ## 1:n

  ## Cost matrix:
  ## 1 if different, 0 if same particle:
  C <- outer(1:n,ind_resample,  FUN='!=')*1

  ## solve assignment problem:
  hungarian <- solve_LSAP(C)

  return( ind_resample[hungarian] )
}



#' Rearrange indices such that they match 1:n
#' 
#' @description 
#' Slow but same as in fortran90 to test code...
#'
#' @param ind vector of resampling indices
#' @param n ensemble size
#' @return permuted indices such that they match 1:n as much as possible
#' @example 
#' ind <- sample(1:10, 10, replace=TRUE) 
#' reorder_ind(ind, 10)
#' reorder_ind_fast(ind, 10)
#' ## timing:
#' library(microbenchmark)
#' nbench <- 40
#' microbenchmark(
#' reorder_ind(sample(1:nbench, nbench, replace=TRUE),10),
#' reorder_ind_fast(sample(1:nbench, nbench, replace=TRUE), 10)
#' )
reorder_ind <- function(ind, n){
  out <- rep(-1, n)
  #   out_filled <- logical(length(ind))
  #   ind_used   <- logical(length(ind))
  ## loop through the indices to fill
  for (i in 1:length(out)){
    ## look among the available indices:
    for (j in 1:length(ind)){
      ## forall construct in fortran: add mask (ind > 0)
      if (ind[j] > 0){
        if (ind[j]==i){
          out[i] <- i
          ind[j] <- -1
          break
        }
      }

    }
  }

  ## forall where out < 0
  for (i in 1:length(out)){
    if (out[i] < 0){
      ## minloc with mask ind > 0: (silly in R)
      minval <- min( ind[ind>0] )
      out[i] <- minval
      ind[which(ind==minval)[1] ] <- -1
    }
  }

  out
}




# Weps part: --------------------------------------------------------------


#' Compute Weps for stochastic EnKPF
#'
#' @param gam fixed gamma value
#' @inheritParams etkpf_es
#' @return Weps
get_Weps_stochastic <- function( R_evc, R_evl, n, gam, eps ){
  ## where eps is a nxn matrix of iid N(0,1)
  ## center eps:
  epsbar <- apply(eps, 1, mean)
  eps <- eps - epsbar %*% t(rep(1,n))

  del_pug <- ( gam*R_evl )/ ( gam*R_evl^2 + 2*(n-1)*gam*R_evl + (n-1)^2 )
  del_pug <- sqrt(del_pug)
  Pug_sqrt <- R_evc %*% diag(del_pug)  %*% t(R_evc)
  Weps <- Pug_sqrt %*% eps

  return(Weps)
}


#' Compute Weps for Riccati ETKPF
#'
#' @param gam fixed gamma value
#' @param Wmu Wmu as returned from get_Wmu
#' @param ind_resample vector of resampling indices
#' @param tol tolerance for stopping criterion
#' @param maxit maximum number of iteration of Newton's algorithm
#' @inheritParams etkpf_es
#' @return Weps
get_Weps_riccati <- function( R_evc, R_evl, n, gam, Wmu, ind_resample,
                              tol=10^-9, maxit=20){

  ## we solve for: A'X + XA + XX' = Pag
  ## Pag:
  del_pag <- ( gam*R_evl )/ ( gam*R_evl^2 + 2*(n-1)*gam*R_evl + (n-1)^2 )
  del_pag <- del_pag*(n-1)
  Pag <- R_evc %*% diag(del_pag)  %*% t(R_evc)

  if ( all(del_pag <= tol)) return(matrix(0,n,n))


  ## A centered:
  A <- Wmu[,ind_resample] - rowMeans(Wmu[,ind_resample])
  A <- t(A)

  ## special case where only one particle is selected: Weps=Pag^1/2
  if (all(A==0)) {
    X1 <- R_evc %*% diag(sqrt(del_pag))  %*% t(R_evc)
    return(X1)
  }


  ## Newton's solver:
  ## initialization:
  beta <- 1.1 * sqrt(sum(A^2))
  C = A + beta*diag(n)
  E = 2*diag(n)
  X1 <-  lyapunov(C,t(C), E, n)

  for (i in 1:maxit){
    X2 <- lyapunov(t(A) +  X1, A + t(X1), Pag + X1%*%t(X1), n)

    ff <- t(A)%*%X2 + X2%*%A + X2%*%X2 - Pag

    X1 <- X2

    if (1/n^2*sqrt(sum(ff^2)) <= tol) break
  }

  return(X1)
}





#' Lyapunov equation solver (used by get_Weps_riccati):
#' 
#' @description 
#' CX + XD = E
#' in Riccati Newton's step used to solve for X2 given X1:
#'  (A' + X1) X2  + X2  (A + X1) = (Pag + X1 X1')
#'  
#' @param C left multiplier (A' + X1)
#' @param D right multiplier (A + X1)
#' @param E rhs (Pag + X1 X1')
#' @param n ensemble size
#' @return X
lyapunov <- function(C,D,E,n){
  M <- matrix(0, n*n, n*n)
  b <- matrix(0, n*n, 1)
  l <- 0
  for (j in 1:n){
    for (i in 1:n){
      l <-  l+1
      b[l] <-  E[i,j]
      cind <-  ((j-1)*n + 1):(j*n)
      M[l, cind ] = C[i,]
      dind = seq(i, n*n, by=n)
      M[l, dind ] = M[l, dind ] +  t(D[,j])
    }
  }
  z <- solve(M, b)

  X <- matrix(z, n, n)
  return(X)
}




