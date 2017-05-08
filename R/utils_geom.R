

#' computes a symmetric matrix square root of a positive definite matrix
#' use eigenvalue decompostion
#'
#' @param A a symmetric matrix
#' @return S such that S'S = A
msqrt <- function(A){
  if(ncol(A)==1) return(sqrt(A))
  e <- eigen(A)
  V <- e$vectors
  return(V %*% diag(sqrt(e$values)) %*% t(V))
}


#' correlation function from Gaspari & Coh 1999, eq. 4.10
#' used in GC_taper
#'
#' @param z is a distance
#' @param c the support half-length
#' @family sweq.manip
#' @examples
#' GC_function(2, 5)
#' sapply(1:10, GC_function, c=5)
GC_function <- function(z, c){
  ## correlation function
  ##
  ##
  ## c is the support half-length
  zc <- z/c
  if (z >= 0 & z <= c){
    -1/4 * zc^5 + 1/2 * zc^4 + 5/8 * zc^3 - 5/3 * zc^2 + 1
  } else if (z > c & z <= 2*c) {
    1/12 * zc^5 - 1/2 * zc^4 + 5/8 * zc^3 + 5/3 * zc^2 - 5 * zc + 4 - 2/3 * c/z
  } else {
    0
  }
}


#' Compute the GC taper on a ring, with support half-length c
#'
#' @param ndim
#' @param c the support half-length
#' @return taper matrix
#' @examples
#' taper <- ring_GC(100, 10)
#' image_mat(taper)
ring_GC <- function(ndim, c){
  if (c*4 > ndim) warning("no zeroes in taper")
  dd <- ring_alldist(ndim)
  apply(dd, c(1,2), GC_function, c=c)
}



#' Compute minimal distances (discrete) on a ring...
#'
#' @param i from
#' @param j to
#' @param N size of the ring
#' @return distance on a ring between i and j
#' @examples
#' ring_dist(1, 9, 10)
ring_dist <- function(i, j, N){
  abs(N/2 - abs(abs(i-j) - N/2))
}


#' Compute all pairwise distances on a rings
#'
#' @param N the size of the ring
#' @return the matrix of pairwise distances on a ring of size N.
#' @examples
#' A <- ring_alldist(9)
ring_alldist <- function(N){
  ind <- 1:N
  A <- matrix(NA, nrow=N, ncol=N)
  for (i in ind){
    A[i,] <- ring_dist(i, ind, N)
  }
  return(A)
}


