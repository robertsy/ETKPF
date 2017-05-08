
!==========================================================================
subroutine etkpf (R_evc, R_evl, C, n, W, wabar, Wmu, wmubar, Weps, &
                  alphai, ind_resample, unif,   &
                  weps_type, adaptive, center, &
                  gam, ess0, ess1, imin, imax,  &
                  eps, eps_cv)
  use etkpf_util
!   use mo_matrix,       only: eigen          ! real sym. matrix eigenvec. dec.

  
  implicit none
  
  integer  ,intent(in)  :: n           ! ensemble size
  real(wp) ,intent(in)  :: R_evc (n,n) ! eigenvectors of Yt R^-1 Y
  real(wp) ,intent(in)  :: R_evl (n)   ! eigenvalues  of Yt R^-1 Y
  real(wp) ,intent(in)  :: C     (n)   ! Yt R^-1 (o-H(xb_mean))
  
  
  real(wp) ,intent(out) :: W (n,n)     ! xa = xabar + Xb W
  real(wp) ,intent(out) :: wabar(n)    ! xabar = xbbar + Xb wabar
  real(wp) ,intent(out) :: Wmu (n,n)   ! M = Xb Wmu !the centroids spread
  real(wp) ,intent(out) :: wmubar (n)  ! weights of each particle
  real(wp) ,intent(out) :: Weps(n,n)   ! pertubrations matrix 
  
  real(wp) ,intent(out) :: alphai (n)  ! weights of each particle
  integer  ,intent(out) :: ind_resample(n)! vector of indices
  real(wp), intent(in)  :: unif        ! uniform for balanced sampling
  
  integer  ,intent(in)  :: weps_type   ! 0: stochastic (in ensemble space)
                                       ! 1: Riccati solver
  integer  ,intent(in)  :: adaptive    ! adaptive choice of gamma:
                                       ! 0: not adaptive
                                       ! 1: in ess bounds
                                       ! 2: min MSE
                                       ! 3: min ES
                                       ! 4: min CV(MSE) 
                                       ! 5: min CV(ES)      
  integer  ,intent(in)  :: center      ! if 1 parametrize such that: (necessary for cov inflation)
                                       ! xabar = xbbar + Xb wm
                                       ! xa = xabar + Xb W                                                
  ! parameters for adaptive gamma:            
  real(wp),intent(inout):: gam         ! gamma parameter for ENKPF                          
  real(wp) ,intent(in)  :: ess0,ess1   ! ess bounds
  real(wp) ,intent(in)  :: imin, imax  ! default parameters for adaptive gamma increments

  ! parameters for Weps:
  real(wp) ,intent(in)  :: eps(n,n)    ! standard gaussian
  real(wp) ,intent(in)  :: eps_cv(n,n) ! standard gaussian
  
                                                              
  ! dummy variables:
  real(wp)              :: rvl(n)                       
  integer               :: ind_temp(n)! vector of indices
  real(wp)              :: ess
  real(wp)              :: shift(n)
  integer               :: i
                                       

  rvl = 0._wp; where (R_evl >  1.e-5) rvl = R_evl

  ! adaptive gamma:
  select case (adaptive)
    case (0) ! not adaptive
      call get_alpha( R_evc, rvl, C, n, alphai, gam, ess )
    case (1) ! ess in bounds
        call adaptive_gamma_ess(R_evc, rvl, C, n, alphai, gam, ess0, ess1, &
                    imin, imax)
    case (2) ! min MSE criterion
        call adaptive_gamma_minMSE(R_evc, rvl, C, n, alphai, gam, &
                    imin, imax, unif)
    case (3) ! min ES criterion:
        call adaptive_gamma_minES (R_evc, rvl, C, n, alphai, gam, &
                    imin, imax, unif)
    case (4) ! min CV(MSE)
        call adaptive_gamma_minMSECV(R_evc, rvl, C, n, alphai, gam, &
                    imin, imax, unif, eps_cv)        
    case (5) ! min CV(ES)
        call adaptive_gamma_minESCV(R_evc, rvl, C, n, alphai, gam, &
                    imin, imax, unif, eps_cv)                            
  end select 
                              

  !compute Wmu:
  call get_Wmu( R_evc, rvl, n, Wmu, gam)
  call get_wmubar (R_evc, rvl, C, n, wmubar, gam)


  !resample:
  call bal_sample (alphai, n, ind_temp, unif)
  ! reorder such that indices match postion whenever possible
  call reorder_ind(ind_temp, ind_resample,n)


  !compute Weps:
  select case (weps_type)
    case (0) ! stochastic weps:
      call get_Weps_stoch  (R_evc, rvl, n, Weps, gam, eps)
    case (1) ! Riccati:  
      call get_Weps_riccati(R_evc, rvl, n, Weps, gam, Wmu, ind_resample)
  end select


  !resampled W + Weps:
  W = Wmu(:, ind_resample) + Weps 

  !center such that xbbar + Xb wabar = xabar
  if (center==1) then
    shift = sum( W, DIM=2 )/n
    wabar = wmubar + shift
    do i=1,n
      W(:,i) = W(:,i) - shift
    end do
  else
    wabar = wmubar
  end if
  
end subroutine etkpf






!==========================================================================
subroutine get_Wmu (R_evc, R_evl, n, Wmu, gam)
  use etkpf_util
  implicit none
  
  real(wp) ,intent(in)  :: gam         ! gamma parameter for ENKPF
  integer  ,intent(in)  :: n           ! ensemble size
  real(wp) ,intent(in)  :: R_evc (n,n) ! eigenvectors of Yt R^-1 Y
  real(wp) ,intent(in)  :: R_evl (n)   ! eigenvalues  of Yt R^-1 Y
  real(wp) ,intent(out) :: Wmu (n,n)   ! M = Xb Wmu !the centroids spread
  real(wp)              :: del_mu (n)  ! eigenvalues of Wmu
  integer               :: i

  del_mu = ( (n-1)*gam*R_evl + (n-1)**2) / (gam*R_evl**2 + 2*(n-1)*gam*R_evl + (n-1)**2)
  !construct Wmu:
  do i=1,n
    Wmu (:,i) = R_evc(:,i) * del_mu(i)
  end do
  Wmu = matmul (Wmu, transpose (R_evc))
  
end subroutine get_Wmu




!==========================================================================
subroutine get_wmubar(R_evc, R_evl, C, n, wmubar, gam)
  use etkpf_util
  implicit none
  
  real(wp) ,intent(in)  :: gam         ! gamma parameter for ENKPF
  integer  ,intent(in)  :: n           ! ensemble size
  real(wp) ,intent(in)  :: R_evc (n,n) ! eigenvectors of Yt R^-1 Y
  real(wp) ,intent(in)  :: R_evl (n)   ! eigenvalues  of Yt R^-1 Y
  real(wp) ,intent(in)  :: C     (n)   ! Yt R^-1 (o-H(xb_mean))
  real(wp) ,intent(out) :: wmubar (n)  ! weights of each particle
  real(wp)              :: del_wmubar (n)
  real(wp)              :: mat1 (n, n)
  integer               :: i

  del_wmubar = ( gam + ( (n-1)*gam*(1-gam)*R_evl / & 
                         (gam*R_evl**2 + 2*(n-1)*gam*R_evl + (n-1)**2)  ) )/ &
               ( (n-1) + gam*R_evl )

  !construct mat1:
  do i=1,n
    mat1 (:,i) = R_evc(:,i) * del_wmubar(i)
  end do
  mat1 = matmul (mat1, transpose (R_evc))

  !construct wmubar:
  wmubar = matmul (mat1, C )


end subroutine get_wmubar








  !==========================================================================
  subroutine get_alpha(R_evc, R_evl, C, n, alphai, gam, ess)
    use etkpf_util
    implicit none
    
    real(wp) ,intent(in)  :: gam         ! gamma parameter for ENKPF
    integer  ,intent(in)  :: n           ! ensemble size
    real(wp) ,intent(in)  :: R_evc (n,n) ! eigenvectors of Yt R^-1 Y
    real(wp) ,intent(in)  :: R_evl (n)   ! eigenvalues  of Yt R^-1 Y
    real(wp) ,intent(in)  :: C     (n)   ! Yt R^-1 (o-H(xb_mean))
    real(wp) ,intent(out) :: alphai (n)  ! weights of each particle
    real(wp) ,intent(out) :: ess         ! return equivalent sample size
    real(wp)              :: del_alpha (n)
    real(wp)              :: mat1(n,n), mat2(n,n)
    integer               :: i

    del_alpha = ( (n-1._wp)**2 * (1-gam) ) / ( gam*R_evl**2 + 2*(n-1)*gam*R_evl + (n-1._wp)**2 )

    !construct mat1:
    do i=1,n
      mat1 (:,i) = R_evc(:,i) * del_alpha(i) * R_evl(i)
    end do
    mat1 = matmul (mat1, transpose (R_evc))
    
    !construct mat2:
    do i=1,n
      mat2 (:,i) = R_evc(:,i) * del_alpha(i)
    end do
    mat2 = matmul (mat2, transpose (R_evc))

    alphai = diag(mat1) - 2* matmul(mat2, C)
    ! catch some overflow cases:   
!     where (alphai >= 100)  alphai = 100._wp
!     where (alphai <= -100) alphai = -100._wp
!     alphai = exp( -.5_wp * alphai)
    alphai = exp( -.5_wp * (alphai - minval(alphai)))
  
    if (sum(alphai) <=  1.e-12 ) then
      if (gam > 0.001) then
       write(6,*) 'SYLVAIN: small weights, set all equal (gam=', gam,')'
       end if
      alphai = 1._wp
    end if
    
    alphai = alphai/sum(alphai)

    ess = 1._wp / sum(alphai**2)/ n

  end subroutine get_alpha



  !==========================================================================
  subroutine bal_sample (pa, n, ind, unif)
    use etkpf_util
    implicit none

    integer  ,intent(in)              :: n
    real(wp) ,intent(in)              :: pa (n) ! such that Pr(x = x^i) = pa(i)
    real(wp) ,intent(in)              :: unif   ! uniform number
    integer  ,intent(out)             :: ind (n)! return vector of indices

    !! dummy variables:
    real(wp)                          :: weights (n), cum_weights (n), temp (n)
    integer                           :: i,j,k

    !normalize the weights to probabilities:
    weights = pa/sum(pa)

    !compute the cumulative weights:
    cum_weights(1) = weights(1)
    do i=2,n
      cum_weights(i) = cum_weights(i-1) + weights(i)
    end do

    !number of samples by particle:
    temp = floor( n * cum_weights + unif )
    temp(2:n) = temp(2:n) - temp(1:n-1) 
    
    !produce the vector of indices
    !where each indice appears temp(i) times
    k=1
    do i=1,n
      do j=1, int( temp(i) )
        ind(k) = i
        k=k+1
      end do
    end do

  end subroutine bal_sample








  !==========================================================================
  ! Reorder indices in meaningful way
  ! brute force implementation
  ! 
  ! take a vector of ind and reorder it such that:
  ! 1. if possible indices match with their position
  ! 2. if not, keep indices in increasing order
  ! example:
  !  input: 1,1,2,3,4,4
  ! output: 1,2,3,4,1,4
  subroutine reorder_ind(ind_in,out,n)
    use etkpf_util
    implicit none
 
    integer, intent(in)     :: n
    integer, intent(in)     :: ind_in(n)
    integer, intent(out)    :: out(n)
    integer                 :: ind(n)
    integer                 :: i,j, minind
    
    ind = ind_in ! dummy copy of ind_in
    
    ! out = -1 means that it has not been filled yet
    out = -1
    ! ind = -1 means that it has been used already
    
    ! for each element in out
    do i = 1,n
      !look among the available indices:
      do j = 1,n
      ! if the current position in in ind, then fill in out
      ! mark that ind(j) was used by setting it to -1
        if ( ind(j) == i ) then 
          out(i) = i
          ind(j) = -1
          exit
        end if
      end do
    end do
    
    ! the remaining positions that need to be filled in out are -1
    ! the available indices are in ind (ind>0)
    ! for each missing position in out, take the minimum available ind left
    do i = 1,n
      if (out(i) < 0) then
        minind = minloc(ind, 1, ind > 0)
        out(i) = ind(minind)
        ind(minind) = -1
      end if
    end do

  end subroutine reorder_ind


 !==========================================================================
 subroutine adaptive_gamma_ess(R_evc, R_evl, C, n, alphai, gam, ess0, ess1, imin_p, imax_p)
    use etkpf_util
    implicit none

    !find minimum gamma such that ESS is between ess0 and ess1
    real(wp) ,intent(out) :: gam         ! gamma parameter for ENKPF
    integer  ,intent(in)  :: n           ! ensemble size
    real(wp) ,intent(in)  :: R_evc (n,n) ! eigenvectors of Yt R^-1 Y
    real(wp) ,intent(in)  :: R_evl (n)   ! eigenvalues  of Yt R^-1 Y
    real(wp) ,intent(in)  :: C     (n)   ! Yt R^-1 (o-H(xb_mean))
    real(wp) ,intent(out) :: alphai(n)   ! weights of each particle
    real(wp) ,intent(in)  :: ess0, ess1  ! bounds of ESS
    real(wp), intent(in)  :: imin_p, imax_p! fraction of gamma to test
    real(wp)              :: ess         ! chosen ESS
    real(wp)              :: imin, imax, imid
    real(wp)              :: delta_gam   ! increments to test
    real(wp)              :: best_alphai(n)
    real(wp)              :: best_ess
    real(wp)              :: best_gam
    
    
    ! search gamma for fractions:
    ! imin/ (imax-imin) to imax/ (imax-imin)
    !initialize imin and imax:
    imin = imin_p
    imax = imax_p

    !in original paper:
    !imin = 1._wp; imax = 15._wp
    !delta_gam = 1/15._wp
    delta_gam = 1._wp/imax
    
    best_ess = 2._wp ! best ess so far (max=1)
    
    
    do while ( imax >= imin )
      imid = ceiling( imin + (imax-imin)/2._wp )
      gam = imid * delta_gam
      call get_alpha(R_evc, R_evl, C, n, alphai, gam, ess)
      
      ! if ess is smaller than accepted, increase gam:
      if(ess < ess0) then !increase gam, no change in best_ess
        imin = imid + 1
      else
        ! if ess is bigger than best ess so far, decrease gam:
        if (ess > best_ess) then
          imax = imid - 1
          ! otherwise (ess < best_ess and ess > ess0), update best_ess and best_model
        else 
          best_ess = ess
          best_alphai = alphai
          best_gam = gam
          ! if best_ess > ess1, decrease gam
          if (ess > ess1) then
            imax = imid - 1
          else  ! Exit: set imax such that it breaks out of the loop:
            imax = -1
          end if
        end if
      end if

    end do
    
    ! set values to the best:
    gam = best_gam
    ess = best_ess
    alphai = best_alphai
    
    
  end subroutine adaptive_gamma_ess
  
  
  
  
 !==========================================================================
 subroutine adaptive_gamma_minMSE(R_evc, R_evl, C, n, alphai, gam, imin, imax, unif)
    use etkpf_util
    implicit none
    
    ! find gamma that minimizes Err(gam) = (y-Hxabar)'R^-1 (y-Hxabar) 
    !         equivalently that mimizes:    wm'S wm - 2 wm'C:  
    real(wp) ,intent(out) :: gam         ! gamma parameter for ENKPF
    integer  ,intent(in)  :: n           ! ensemble size
    real(wp) ,intent(in)  :: R_evc (n,n) ! eigenvectors of Yt R^-1 Y
    real(wp) ,intent(in)  :: R_evl (n)     ! eigenvalues  of Yt R^-1 Y
    real(wp) ,intent(in)  :: C     (n)   ! Yt R^-1 (o-H(xb_mean))
    real(wp) ,intent(out) :: alphai(n)   ! weights of each particle
    real(wp), intent(in)  :: imin, imax  ! fraction of gamma to test
    real(wp), intent(in)  :: unif
    real(wp)              :: delta_gam   ! increments to test
    real(wp)              :: current_alphai(n), current_gam
    real(wp)              :: current_error, best_error
    real(wp)              :: current_C(n)
    real(wp)              :: Wmu(n,n), wm(n)
    real(wp)              :: ess
    integer               :: ind_resample(n)
    integer               :: j

    delta_gam = 1._wp/imax
    
    ! we start at gamma=1 and go to 0 by increments of delta_gam
    current_gam = 1._wp
    
    ! loop through all gammas:
    gamloop: do while ( current_gam >= 0._wp)
    
        ! analysis for current_gam: returns wm with resampling correction:
        call get_alpha( R_evc, R_evl, C, n, current_alphai, current_gam, ess )
        call get_Wmu(   R_evc, R_evl, n, Wmu, current_gam)
        call get_wmubar (R_evc, R_evl, C, n, wm, current_gam)
        call bal_sample (current_alphai, n, ind_resample, unif)
        ! shift due to resampling:
        wm = wm + sum( Wmu(:, ind_resample), DIM=2 )/n ! == rowMeans(Wmu * Walpha)
      
        
        ! compute Err(gam) = wm'S wm - 2 wm'C:
        current_C = matmul( transpose (R_evc), wm)
        do j=1,n
          current_C(j) = current_C(j) * R_evl(j)
        end do
        current_C = matmul( R_evc, current_C)
        current_error = dot_product(wm, current_C)
        current_error = current_error - 2 * dot_product(wm,C)

        ! if improvement of error OR first iteration:
        if ( (current_error < best_error) .OR. current_gam == 1._wp ) then
          best_error = current_error
          gam = current_gam
          alphai = current_alphai
        end if
        
        ! increment gamma:
        current_gam =  current_gam - delta_gam
   
   
        ! to be sure we round down to 0
        if (abs( current_gam ) <= 0.1_wp * delta_gam ) then
          current_gam = 0._wp
        end if
      
      
      end do gamloop
    
   
  end subroutine adaptive_gamma_minMSE
  
  
  
  
   !==========================================================================
   subroutine adaptive_gamma_minES(R_evc, R_evl, C, n, alphai, gam, imin, imax, unif)
    use etkpf_util
    implicit none
    
    
    ! find gamma that minimizes ES(y) in ensemble space
    real(wp) ,intent(out) :: gam         ! gamma parameter for ENKPF
    integer  ,intent(in)  :: n           ! ensemble size
    real(wp) ,intent(in)  :: R_evc (n,n) ! eigenvectors of Yt R^-1 Y
    real(wp) ,intent(in)  :: R_evl (n)     ! eigenvalues  of Yt R^-1 Y
    real(wp) ,intent(in)  :: C     (n)   ! Yt R^-1 (o-H(xb_mean))
    real(wp) ,intent(out) :: alphai(n)   ! weights of each particle
    real(wp), intent(in)  :: imin, imax  ! fraction of gamma to test
    real(wp), intent(in)  :: unif
    real(wp)              :: delta_gam   ! increments to test
    real(wp)              :: current_alphai(n), current_gam, equal_alphai(n)
    real(wp)              :: current_error, best_error
    real(wp)              :: current_C(n)
    real(wp)              :: Wmu(n,n), wm(n)
    real(wp)              :: ess
    integer               :: ind_resample(n)
    integer               :: i,j
    real(wp)              :: rvl_inv(n)                 ! to compute Sinv:
    real(wp)              :: del_pag(n)                 ! to compute Pagens:
    real(wp)              :: crps_vec(n)                ! to compute CRPS:
    real(wp)              :: wy(n) ! y projection
    
    
    delta_gam = 1._wp/imax
    
    ! we start at gamma=1 and go to 0 by increments of delta_gam
    current_gam = 1._wp
!     current_gam = 0._wp
    
    
    ! compute Sinv:
    rvl_inv = 0._wp
    where (R_evl > 1.e-5) rvl_inv = 1._wp/R_evl
    
    equal_alphai = 1._wp / n
    
    ! y proj = Yb S^-1 C:
    wy = matmul(transpose(R_evc), C) * sqrt(rvl_inv)
    
   
    ! loop through all gammas:
    gamloop: do while ( current_gam >= 0._wp)
    
        ! analysis for current_gam: returns wm with resampling correction:
        call get_alpha( R_evc, R_evl, C, n, current_alphai, current_gam, ess )
        call get_Wmu(   R_evc, R_evl, n, Wmu, current_gam)
        call get_wmubar (R_evc, R_evl, C, n, wm, current_gam)
        call bal_sample (current_alphai, n, ind_resample, unif)
        ! shift due to resampling:
!         wm = wm + sum( Wmu(:, ind_resample), DIM=2 )/n ! == rowMeans(Wmu * Walpha)

        ! resampling:
        Wmu = Wmu(:,ind_resample)
      
        ! compute Pagens:
        del_pag = ( current_gam*R_evl ) / (current_gam*R_evl**2 + 2*(n)*current_gam*R_evl + (n)**2 )       
       
       

        
        ! in orthonormal basis:
        do i=1,n
          Wmu(:,i) = Wmu(:,i) + wm 
        end do

        Wmu = matmul(transpose(R_evc), Wmu) 
        
        do i=1,n
          Wmu(i,:) = Wmu(i,:) * sqrt(R_evl(i)) 
        end do
        

        
        ! compute CRPS componentwise:
        do i=1,n
          call f_crps( Wmu(i,:), R_evl(i) * del_pag(i) + 1._wp , equal_alphai, wy(i), n, crps_vec(i))
        end do
        
        current_error = sum(crps_vec) 
        
!         if (abs(current_gam - 0.5_wp ) < 0.001_wp) then
!         write(6,*) current_gam
!             write (6,*) current_error
! !             write (6,*) 'test', Wmu(2,:)
!         end if 
        
        ! if improvement of error OR first iteration:
        if ( current_error <= best_error .OR. current_gam == 1._wp ) then
          gam = current_gam
          best_error = current_error
          alphai = current_alphai
        end if
        
        

        
        ! increment gamma:
        current_gam =  current_gam - delta_gam
        
        ! to be sure we round down to 0
        if (abs( current_gam ) <= 0.1_wp * delta_gam ) then
          current_gam = 0._wp
        end if
        
      end do gamloop
    
   
  end subroutine adaptive_gamma_minES



!==========================================================================
 subroutine adaptive_gamma_minMSECV(R_evc, R_evl, C, n, alphai, gam, imin, imax, unif, eps)
    use etkpf_util
    use mo_matrix,               only: check_rs

    implicit none
    
    ! assimilate fake data yj = Hxbj + etaj
    real(wp) ,intent(out) :: gam         ! gamma parameter for ENKPF
    integer  ,intent(in)  :: n           ! ensemble size
    real(wp) ,intent(in)  :: R_evc (n,n) ! eigenvectors of Yt R^-1 Y
    real(wp) ,intent(in)  :: R_evl (n)     ! eigenvalues  of Yt R^-1 Y
    real(wp) ,intent(in)  :: C     (n)   ! Yt R^-1 (o-H(xb_mean))
    real(wp) ,intent(out) :: alphai(n)   ! weights of each particle
    real(wp), intent(in)  :: imin, imax  ! fraction of gamma to test
    real(wp), intent(in)  :: unif
    real(wp) ,intent(in)  :: eps(n,n)    ! standard gaussian

    ! elements for CV-loop (steps, best values, etc):
    real(wp)              :: delta_gam   ! increments to test
    real(wp)              :: current_gam
    real(wp)              :: current_error, best_error, cv_error(n)

    ! elements for analysis in CV-loop:
    real(wp)              :: Wmu(n-1,n-1), wmubar(n-1), wm(n-1), alphai_j(n-1)
    integer               :: ind_resample(n-1) ! to resample the n-1 particles

    ! elements for S and C in CV-loop:    
    real(wp)              :: Sfull(n,n), S_j(n-1,n-1)  ! S= Yb' R-1 Yb
    real(wp)              :: Stemp(n, n-1)
    real(wp)              :: R_evc_j (n-1,n-1)
    real(wp)              :: R_evl_j(n-1), rvl_j (n-1)
    real(wp)              :: C_j (n-1)
    real(wp)              :: e_xij(n)
    real(wp)              :: Ctemp (n)

    ! elements to computer error:    
    real(wp)              :: werr(n-1)


    ! elements for selecting subsets of variables:
    integer               :: ind_cv(n-1)               ! to select all but j
    real(wp)              :: A_j (n, n-1)              ! Xb-j = Xb A_j
    
    
    real(wp)              :: rvl_inv(n), xij(n,n) !for y_j = Hxbj + eta^j
                                                  !        = Hxbj + HXb xij + rj
    real(wp)              :: ess  
    integer               :: i,j,k,kk


      

    if (sum(R_evl) <= 1.e-5 ) then
      gam = 0._wp
      call get_alpha( R_evc, R_evl, C, n, alphai, gam, ess )
    else

      delta_gam = 1._wp/imax
    
      ! we start at gamma=1 and go to 0 by increments of delta_gam
      current_gam = 1._wp
    
      ! yj = Hxbj + etaj => etaj = HXb xij (+ residual)
      rvl_inv = 0._wp
      where (R_evl > 1.e-5) rvl_inv = 1._wp/R_evl
      do i=1,n
        if (rvl_inv(i) < 0._wp ) then ! to be extra sure
          rvl_inv(i) = 0._wp
        end if    
        xij (:,i) = R_evc(:,i) * sqrt( rvl_inv(i) )
      end do
      xij = matmul (xij, transpose (R_evc))
      xij = matmul (xij, eps)
    
      ! construct S:
      do i=1,n
        Sfull (:,i) = R_evc(:,i) * R_evl(i)
      end do
      Sfull = matmul (Sfull, transpose (R_evc))
      
      ! loop through all gammas:
      gamloop: do while ( current_gam >= 0._wp )
                            
          cvloop: do j=1,n 
            ! all indices except j:
            kk = 1
            do k=1,n
              if ( k /= j) then 
                ind_cv(kk) = k
                kk = kk + 1
              end if
            end do
            
            ! construct A_j:
            A_j = 0._wp
            kk = 1
            do i=1,n
              if ( i /= j ) then
                A_j (i,kk) = 1._wp
                kk = kk + 1
                else
                A_j (i,:) = 1._wp / (n - 1._wp)
              end if 
            end do
            

            
            ! construct S_j:
            Stemp = matmul( Sfull, A_j)
            S_j = matmul( transpose(A_j), Stemp)          
            
            ! eigenvalue decomposition of S_j:
!             call myeigen ( S_j, R_evc_j, R_evl_j, n-1, INFO )            
!             if (INFO /= 0) write (6,*) 'SYLVAIN:: fatal error in eigen ', INFO
           
!             call eigen ( S_j, R_evc_j, R_evl_j)
            call check_rs (S_j, evc=R_evc_j, evl=R_evl_j)

            
            rvl_j = 0._wp; where (R_evl_j >  1.e-5) rvl_j = R_evl_j


            ! construct C_j:
            e_xij = xij(:,j)
            e_xij(j) = e_xij(j) + n/(n-1._wp)
            Ctemp = matmul( Sfull, e_xij )
            C_j = matmul (transpose(A_j), Ctemp)
            

            
            ! compute analysis (Wmu and ind_resample):            
            call get_alpha( R_evc_j, rvl_j, C_j, n-1, alphai_j, current_gam, ess )
            call get_Wmu(   R_evc_j, rvl_j, n-1, Wmu, current_gam)
            call get_wmubar(R_evc_j, rvl_j, C_j, n-1, wmubar, current_gam)                        
            call bal_sample (alphai_j, n-1, ind_resample, unif)
            

            
            ! compute analysis mean wm:
            wm = wmubar + sum( Wmu(:, ind_resample), DIM=2 )/(n-1) ! == rowMeans(Wmu * Walpha)



            
            ! compute MSE:                
            werr = matmul(S_j, wm)
            cv_error(j) = dot_product(wm, werr) - 2 * dot_product(wm, C_j)          
                             
          
          end do cvloop

      
        current_error = sum(cv_error)


        ! if improvement of error OR first iteration:
        if ( (current_error < best_error) .OR. current_gam == 1._wp ) then
          best_error = current_error
          gam = current_gam
          call get_alpha( R_evc, R_evl, C, n, alphai, gam, ess )
        end if
        

        ! increment gamma:
        current_gam =  current_gam - delta_gam
      
        ! to be sure we round down to 0
        if (abs( current_gam ) <= 0.1_wp * delta_gam ) then
          current_gam = 0._wp
        end if
         
      end do gamloop
    
  end if   

    
!write(6,*) 'SYLVAIN best gam', gam, 'best error', best_error
    
  end subroutine adaptive_gamma_minMSECV
  
  
  
  !==========================================================================
 subroutine adaptive_gamma_minMSECVearly(R_evc, R_evl, C, n, alphai, gam, imin, imax, unif, eps)
    use etkpf_util
    use mo_matrix,               only: check_rs

    implicit none
    
    ! assimilate fake data yj = Hxbj + etaj
    real(wp) ,intent(out) :: gam         ! gamma parameter for ENKPF
    integer  ,intent(in)  :: n           ! ensemble size
    real(wp) ,intent(in)  :: R_evc (n,n) ! eigenvectors of Yt R^-1 Y
    real(wp) ,intent(in)  :: R_evl (n)     ! eigenvalues  of Yt R^-1 Y
    real(wp) ,intent(in)  :: C     (n)   ! Yt R^-1 (o-H(xb_mean))
    real(wp) ,intent(out) :: alphai(n)   ! weights of each particle
    real(wp), intent(in)  :: imin, imax  ! fraction of gamma to test
    real(wp), intent(in)  :: unif
    real(wp) ,intent(in)  :: eps(n,n)    ! standard gaussian

    ! elements for CV-loop (steps, best values, etc):
    real(wp)              :: delta_gam   ! increments to test
    real(wp)              :: current_gam
    real(wp)              :: current_error, best_error, cv_error(n)
real(wp)              :: early_max, early_jump ! for stopping criteria

    ! elements for analysis in CV-loop:
    real(wp)              :: Wmu(n-1,n-1), wmubar(n-1), wm(n-1), alphai_j(n-1)
    integer               :: ind_resample(n-1) ! to resample the n-1 particles

    ! elements for S and C in CV-loop:    
    real(wp)              :: Sfull(n,n), S_j(n-1,n-1)  ! S= Yb' R-1 Yb
    real(wp)              :: Stemp(n, n-1)
    real(wp)              :: R_evc_j (n-1,n-1)
    real(wp)              :: R_evl_j(n-1), rvl_j (n-1)
    real(wp)              :: C_j (n-1)
    real(wp)              :: e_xij(n)
    real(wp)              :: Ctemp (n)

    ! elements to computer error:    
    real(wp)              :: werr(n-1)


    ! elements for selecting subsets of variables:
    integer               :: ind_cv(n-1)               ! to select all but j
    real(wp)              :: A_j (n, n-1)              ! Xb-j = Xb A_j
    
    
    real(wp)              :: rvl_inv(n), xij(n,n) !for y_j = Hxbj + eta^j
                                                  !        = Hxbj + HXb xij + rj
    real(wp)              :: ess  
    integer               :: i,j,k,kk


!--------------------------------------------------------
! decrease gamma until error is > best_error + 2*early_jump
early_jump = 0._wp

      

    if (sum(R_evl) <= 1.e-5 ) then
      gam = 0._wp
      call get_alpha( R_evc, R_evl, C, n, alphai, gam, ess )
    else

      delta_gam = 1._wp/imax
    
      ! we start at gamma=1 and go to 0 by increments of delta_gam
      current_gam = 1._wp
    
      ! yj = Hxbj + etaj => etaj = HXb xij (+ residual)
      rvl_inv = 0._wp
      where (R_evl > 1.e-5) rvl_inv = 1._wp/R_evl
      do i=1,n
        if (rvl_inv(i) < 0._wp ) then ! to be extra sure
          rvl_inv(i) = 0._wp
        end if    
        xij (:,i) = R_evc(:,i) * sqrt( rvl_inv(i) )
      end do
      xij = matmul (xij, transpose (R_evc))
      xij = matmul (xij, eps)
    
      ! construct S:
      do i=1,n
        Sfull (:,i) = R_evc(:,i) * R_evl(i)
      end do
      Sfull = matmul (Sfull, transpose (R_evc))
      
      ! loop through all gammas:
      gamloop: do while ( current_gam >= 0._wp )
                            
          cvloop: do j=1,n 
            ! all indices except j:
            kk = 1
            do k=1,n
              if ( k /= j) then 
                ind_cv(kk) = k
                kk = kk + 1
              end if
            end do
            
            ! construct A_j:
            A_j = 0._wp
            kk = 1
            do i=1,n
              if ( i /= j ) then
                A_j (i,kk) = 1._wp
                kk = kk + 1
                else
                A_j (i,:) = 1._wp / (n - 1._wp)
              end if 
            end do
            

            
            ! construct S_j:
            Stemp = matmul( Sfull, A_j)
            S_j = matmul( transpose(A_j), Stemp)          
            
            ! eigenvalue decomposition of S_j:
!             call myeigen ( S_j, R_evc_j, R_evl_j, n-1, INFO )            
!             if (INFO /= 0) write (6,*) 'SYLVAIN:: fatal error in eigen ', INFO
           
!             call eigen ( S_j, R_evc_j, R_evl_j)
            call check_rs (S_j, evc=R_evc_j, evl=R_evl_j)

            
            rvl_j = 0._wp; where (R_evl_j >  1.e-5) rvl_j = R_evl_j


            ! construct C_j:
            e_xij = xij(:,j)
            e_xij(j) = e_xij(j) + n/(n-1._wp)
            Ctemp = matmul( Sfull, e_xij )
            C_j = matmul (transpose(A_j), Ctemp)
            

            
            ! compute analysis (Wmu and ind_resample):            
            call get_alpha( R_evc_j, rvl_j, C_j, n-1, alphai_j, current_gam, ess )
            call get_Wmu(   R_evc_j, rvl_j, n-1, Wmu, current_gam)
            call get_wmubar(R_evc_j, rvl_j, C_j, n-1, wmubar, current_gam)                        
            call bal_sample (alphai_j, n-1, ind_resample, unif)
            

            
            ! compute analysis mean wm:
            wm = wmubar + sum( Wmu(:, ind_resample), DIM=2 )/(n-1) ! == rowMeans(Wmu * Walpha)



            
            ! compute MSE:                
            werr = matmul(S_j, wm)
            cv_error(j) = dot_product(wm, werr) - 2 * dot_product(wm, C_j)          
                             
          
          end do cvloop

      
        current_error = sum(cv_error)


        ! if improvement of error OR first iteration:
        if ( current_error <= best_error + .5_wp * early_jump .OR. current_gam == 1._wp ) then
!         if ( (current_error < best_error) .OR. current_gam == 1._wp ) then
          best_error = current_error
          gam = current_gam
          call get_alpha( R_evc, R_evl, C, n, alphai, gam, ess )
        end if
        
        ! monitor variation in the first few steps
          if ( current_gam >= min( 0.8_wp, 1 - delta_gam ) ) then
            if ( current_error > early_max .OR. current_gam == 1._wp) then
              early_max = current_error
            end if
            early_jump = early_max - best_error
          end if


        ! increment gamma:
        current_gam =  current_gam - delta_gam
      
        ! to be sure we round down to 0
        if (abs( current_gam ) <= 0.1_wp * delta_gam ) then
          current_gam = 0._wp
        end if
         
      end do gamloop
    
  end if   

    
!write(6,*) 'SYLVAIN best gam', gam, 'best error', best_error
    
  end subroutine adaptive_gamma_minMSECVearly
  
  
  
  
  
  
!==========================================================================
 subroutine adaptive_gamma_minESCV(R_evc, R_evl, C, n, alphai, gam, imin, imax, unif, eps)
    use etkpf_util
    use mo_matrix,               only: check_rs

    implicit none
    
    ! assimilate fake data yj = Hxbj + etaj
    real(wp) ,intent(out) :: gam         ! gamma parameter for ENKPF
    integer  ,intent(in)  :: n           ! ensemble size
    real(wp) ,intent(in)  :: R_evc (n,n) ! eigenvectors of Yt R^-1 Y
    real(wp) ,intent(in)  :: R_evl (n)     ! eigenvalues  of Yt R^-1 Y
    real(wp) ,intent(in)  :: C     (n)   ! Yt R^-1 (o-H(xb_mean))
    real(wp) ,intent(out) :: alphai(n)   ! weights of each particle
    real(wp), intent(in)  :: imin, imax  ! fraction of gamma to test
    real(wp), intent(in)  :: unif
    real(wp) ,intent(in)  :: eps(n,n)    ! standard gaussian

    ! elements for CV-loop (steps, best values, etc):
    real(wp)              :: delta_gam   ! increments to test
    real(wp)              :: current_gam
    real(wp)              :: current_error, best_error, cv_error(n)
    
    ! elements for analysis in CV-loop:
    real(wp)              :: Wmu(n-1,n-1), wmubar(n-1), wm(n-1), alphai_j(n-1), equal_alphai(n-1)
    real(wp)              :: del_pag(n-1)
    integer               :: ind_resample(n-1) ! to resample the n-1 particles
    real(wp)              :: wy(n-1) ! y projection

    ! elements for S and C in CV-loop:    
    real(wp)              :: Sfull(n,n), S_j(n-1,n-1)  ! S= Yb' R-1 Yb
    real(wp)              :: Stemp(n, n-1)
    real(wp)              :: R_evc_j (n-1,n-1)
    real(wp)              :: R_evl_j(n-1), rvl_j (n-1), rvl_inv_j(n-1)
    real(wp)              :: C_j (n-1)
    real(wp)              :: e_xij(n)
    real(wp)              :: Ctemp (n)

    ! elements to computer error:    
    real(wp)              :: errtemp(n-1)


    ! elements for selecting subsets of variables:
    integer               :: ind_cv(n-1)               ! to select all but j
    real(wp)              :: A_j (n, n-1)              ! Xb-j = Xb A_j
    
    
    real(wp)              :: rvl_inv(n), xij(n,n) !for y_j = Hxbj + eta^j
                                                  !        = Hxbj + HXb xij + rj
    real(wp)              :: ess  
    integer               :: i,j,k,kk


    equal_alphai = 1._wp/(n-1._wp) ! for equal weight CRPS

    if (sum(R_evl) <= 1.e-5 ) then
      gam = 0._wp
      call get_alpha( R_evc, R_evl, C, n, alphai, gam, ess )
    else

      delta_gam = 1._wp/imax
    
      ! we start at gamma=1 and go to 0 by increments of delta_gam
      current_gam = 1._wp
    
      ! yj = Hxbj + etaj => etaj = HXb xij (+ residual)
      rvl_inv = 0._wp
      where (R_evl > 1.e-5) rvl_inv = 1._wp/R_evl
      do i=1,n
        if (rvl_inv(i) < 0._wp ) then ! to be extra sure
          rvl_inv(i) = 0._wp
        end if    
        xij (:,i) = R_evc(:,i) * sqrt( rvl_inv(i) )
      end do
      xij = matmul (xij, transpose (R_evc))
      xij = matmul (xij, eps)
    

    
      ! construct S:
      do i=1,n
        Sfull (:,i) = R_evc(:,i) * R_evl(i)
      end do
      Sfull = matmul (Sfull, transpose (R_evc))
      
      ! loop through all gammas:
      gamloop: do while ( current_gam >= 0._wp )
                            
          cvloop: do j=1,n 
            ! all indices except j:
            kk = 1
            do k=1,n
              if ( k /= j) then 
                ind_cv(kk) = k
                kk = kk + 1
              end if
            end do
            
            ! construct A_j:
            A_j = 0._wp
            kk = 1
            do i=1,n
              if ( i /= j ) then
                A_j (i,kk) = 1._wp
                kk = kk + 1
                else
                A_j (i,:) = 1._wp / (n - 1._wp)
              end if 
            end do
            

            
            ! construct S_j:
            Stemp = matmul( Sfull, A_j)
            S_j = matmul( transpose(A_j), Stemp)          
            
            ! eigenvalue decomposition of S_j:
!             call myeigen ( S_j, R_evc_j, R_evl_j, n-1, INFO )            
!             if (INFO /= 0) write (6,*) 'SYLVAIN:: fatal error in eigen ', INFO
!             call eigen ( S_j, R_evc_j, R_evl_j)
            call check_rs (S_j, evc=R_evc_j, evl=R_evl_j)


            
            rvl_j = 0._wp; where (R_evl_j >  1.e-5) rvl_j = R_evl_j
            
    
            ! compute Sinv:
            rvl_inv_j = 0._wp
            where (rvl_j > 1.e-5) rvl_inv_j = 1._wp/rvl_j
                 


            ! construct C_j:
            e_xij = xij(:,j)
            e_xij(j) = e_xij(j) + n/(n-1._wp)
            Ctemp = matmul( Sfull, e_xij )
            C_j = matmul (transpose(A_j), Ctemp)
            

            
                        ! y proj = Yb S^-1 C:
            wy = matmul(transpose(R_evc_j), C_j) * sqrt(rvl_inv_j)
            

            
            ! compute analysis (Wmu and ind_resample):            
            call get_alpha( R_evc_j, rvl_j, C_j, n-1, alphai_j, current_gam, ess )
            call get_Wmu(   R_evc_j, rvl_j, n-1, Wmu, current_gam)
            call get_wmubar(R_evc_j, rvl_j, C_j, n-1, wmubar, current_gam)                        
            call bal_sample (alphai_j, n-1, ind_resample, unif)
            
            
            
            ! compute analysis mean wm:
            wm = wmubar + sum( Wmu(:, ind_resample), DIM=2 )/(n-1) ! == rowMeans(Wmu * Walpha)


            ! in orthonormal basis:
 
            Wmu = Wmu(:,ind_resample)
            do i=1,n-1
              Wmu(:,i) = Wmu(:,i) + wmubar 
            end do

            Wmu = matmul(transpose(R_evc_j), Wmu) 
        
            do i=1,n-1
              Wmu(i,:) = Wmu(i,:) * sqrt(rvl_j(i)) 
            end do

            del_pag = ( current_gam*rvl_j ) /( current_gam*rvl_j**2 + 2*(n-2._wp)*current_gam*rvl_j + (n-2_wp)**2 )
                      
        
            ! compute ES:
            do i=1,n-1
              call f_crps( Wmu(i,:), rvl_j(i) * del_pag(i) + 1._wp, equal_alphai, wy(i), n-1, errtemp(i))
            end do
        
            cv_error(j) = sum(errtemp) /(n-1._wp)
!                               
!             if (abs(current_gam - 0.5_wp ) < 0.001_wp .AND. j==2) then
!               write (6,*) 'SYLVAIN:: ', cv_error(j)
!             end if                               
          
          end do cvloop

      
        current_error = sum(cv_error)


        ! if improvement of error OR first iteration:
        if ( (current_error < best_error) .OR. current_gam == 1._wp ) then
          best_error = current_error
          gam = current_gam
          call get_alpha( R_evc, R_evl, C, n, alphai, gam, ess )
        end if


        ! increment gamma:
        current_gam =  current_gam - delta_gam
      
        ! to be sure we round down to 0
        if (abs( current_gam ) <= 0.1_wp * delta_gam ) then
          current_gam = 0._wp
        end if
         
      end do gamloop
    
  end if   

    
!write(6,*) 'SYLVAIN best gam', gam, 'best error', best_error
    
  end subroutine adaptive_gamma_minESCV
  
  
  
  
  

  
  
  
  
  
  
  !==========================================================================
  subroutine get_Weps_stoch (R_evc, R_evl, n, Weps, gam, eps)
    use etkpf_util
    implicit none
    
    ! Stochastic version of the filter (in ensemble space)
    ! Weps = (Pag)^1/2 * eps,     eps ~ N(0,1)
    real(wp) ,intent(in)  :: gam         ! gamma parameter for ENKPF
    integer  ,intent(in)  :: n           ! ensemble size
    real(wp) ,intent(in)  :: R_evc (n,n) ! eigenvectors of Yt R^-1 Y
    real(wp) ,intent(in)  :: R_evl (n)   ! eigenvalues  of Yt R^-1 Y
    real(wp) ,intent(in)  :: eps(n,n)  ! standard gaussian
    real(wp) ,intent(out) :: Weps (n,n)  ! the square-root spread
    real(wp)              :: del_pag (n) ! eigenvalues of Pag
    integer               :: i
    real(wp)              :: eps_m (n)    ! mean of eps
    real(wp)              :: eps_c (n,n)  ! centered eps
    
    ! center eps:
    eps_c = eps
    eps_m = sum( eps, DIM=2 )/n
    do i=1,n
      eps_c(:,i) = eps(:,i) - eps_m
    end do
 
    ! compute sqrt:
    del_pag = ( gam*R_evl ) / ( gam*R_evl**2 + 2*(n-1._wp)*gam*R_evl + (n-1._wp)**2 )
    del_pag = sqrt(del_pag)

    !construct Revc * Del * Revc'
    do i=1,n
      Weps (:,i) = R_evc(:,i) *  del_pag(i)
    end do
    Weps = matmul( Weps, transpose(R_evc) )
    
    !Weps =  Revc * Del * Revc' eps
    Weps = matmul(Weps, eps_c)
    
  end subroutine get_Weps_stoch





 !==========================================================================
 subroutine get_Weps_riccati (R_evc, R_evl, n, Weps, gam, Wmu, ind_resample)
    use etkpf_util
    implicit none
 
   ! find Weps solution of: A(Weps)′+ Weps A′+ Weps(Weps)′ =(k−1)Pag
   ! using Newton method with a tolerance level
 
    real(wp) ,intent(in)  :: gam             ! gamma parameter for ENKPF
    integer  ,intent(in)  :: n               ! ensemble size
    real(wp) ,intent(in)  :: R_evc (n,n)     ! eigenvectors of Yt R^-1 Y
    real(wp) ,intent(in)  :: R_evl (n)       ! eigenvalues  of Yt R^-1 Y
    real(wp) ,intent(in)  :: Wmu(n,n)        ! Wmu matrix
    integer  ,intent(in)  :: ind_resample(n) !
    real(wp) ,intent(out) :: Weps (n,n)      ! the solution to Riccati equation
    real(wp)              :: del_pag (n)     ! eigenvalues of Pag
    real(wp)              :: Pag (n,n), A(n,n)
    real(wp)              :: X1(n,n), X2(n,n), resid(n,n)
    real(wp)              :: C(n,n), E(n,n)
    real(wp)              :: mean_resid, beta
    integer               :: i, maxit
    real(wp)              :: tol 
    integer               :: INFO
    
    ! newton stopping criterion:
    maxit = 20
    tol = 1.e-9
    
    
    
    !construct (n-1)Pag = R_evc * diag((n-1) del_pag) *R_evc'
    del_pag = ( gam*R_evl ) / ( gam*R_evl**2 + 2*(n-1._wp)*gam*R_evl + (n-1._wp)**2 )
    do i=1,n
      Pag (:,i) = R_evc(:,i) *  del_pag(i) * (n-1._wp)
    end do
    Pag = matmul( Pag, transpose(R_evc) )



    !-----------------------------------------------------------------
    !Rule out cases which are too close to pure PF:
    purePF: if (maxval(del_pag) < 1.e-9) then
      Weps = 0._wp
    else 
    
    
        !-----------------------------------------------------------------
        ! construct A = Wmu Walpha - 1/k Wmu Walpha 1:
        A = Wmu(:,ind_resample)
        ! equivalently:  A - reshape( sum(A,2)/real(n,wp), shape(A) )
        do i=1,n
          A(i,:) = A(i,:) - sum(A(i,:) )/real(n,wp)
        end do   
        
        
        !-----------------------------------------------------------------
        ! special case where only one particle is selected: Weps=Pag^1/2
        ! all A elements are zero because the colmean of A is the same as each
        ! column of A. 
        oneParticle: if ( sum( abs(A) ) < 1.e-12 ) then
          do i=1,n
            Weps (:,i) = R_evc(:,i) *  sqrt(del_pag(i) * (n-1._wp) )
          end do
          Weps = matmul( Weps, transpose(R_evc) )
        else 
              !-----------------------------------------------------------------
              ! Newton's solver:
              
              ! initialization X1:
              ! here: (A' + beta I)X + X (A + beta I) = 2*I
              beta = sqrt(sum(A**2)) * 1.1_wp
              C = transpose(A) + beta*ident(n)
              E = 2*ident(n)
              call lyap(X1, C, E, n, INFO) ! solves CX + XC' = E
    
              ! Loop until convergence:
              newton: do i =1,maxit
              
                  ! find X2 which solves:
                  ! (A + X1) X2 + X2 (A' + X1') = X1 X1' + Pag
                  C = A + X1
                  E = Pag + matmul(X1, X1)
                  call lyap(X2, C, E, n, INFO) ! solves CX + XC' = E

                  if (INFO /= 0) then
                    write(6,*) 'SYLVAIN:: warning in Newton step,',i,'skip'
                    exit
                  end if
      
                  ! update current solution:
                  X1 = X2
              
                  ! if accuracy level reached, exit the loop:
                  resid = matmul(A,X1) + matmul(X1,transpose(A)) + matmul(X1,X1) - Pag
                  mean_resid = 1._wp/(n**2) * sqrt(sum(resid**2))
!write(*,*) 'newton step ',i,'error',mean_resid                  
                  if ( mean_resid <= tol) then
                    exit
                  end if
  
              end do newton
      
            if (i==maxit .AND. mean_resid > tol) then
              write(6,*) 'SYLVAIN:: warning Newton failed to converged with mean_resid',mean_resid
            end if
        
            ! return Newton solution:
            Weps = X1
  
      ! special cases treated above:
      end if oneParticle
  end if purePF
    
  end subroutine get_Weps_riccati
  
  
  

  
  
  
  
  
  subroutine lyap (X,A,C,n,INFO)
    use etkpf_util
    implicit none
 
    ! solve for X (n,n) symmetric:
    ! AX + XA' = C 
    ! based on Bartels and Stewart (1972)
    real(wp) ,intent(out)     :: X(n,n)
    real(wp) ,intent(in)      :: A(n,n)
    real(wp) ,intent(in)      :: C(n,n)
    integer  ,intent(in)      :: n
    integer  ,intent(out)     :: INFO
    
    ! for the LAPACK solvers:
    integer                 :: LWORK
    real(wp)                :: SCALE
    real(wp) ,allocatable   :: TAU(:), WORK(:)
    real(wp) ,allocatable   :: H(:,:), Q(:,:)                !step 1
    real(wp) ,allocatable   :: Z(:,:), T(:,:), WR(:), WI(:)  !step 2
    real(wp) ,allocatable   :: Xsylv (:,:), CZ(:,:)          !step 3
    real(wp) ,allocatable   :: QZ (:,:)                      !step 4
    
    integer                 :: i,j
    
    ! 1. Compute the Hessenberg form of A:
    ! A = QHQ'
    LWORK = n*n !size of WORK (can check the minimum required when succeeded
    allocate ( WORK (LWORK) , TAU(n) )
    allocate ( H(n,n), Q(n,n) )
    
    ! compute H:
    H = A
    call dgehrd(n, 1, n, H, n, TAU, WORK, LWORK, INFO)
    if (INFO /= 0) write (6,*) 'SYLVAIN:: fatal error in dgehrd ', INFO

    ! compute Q:
    Q = H
    call dorghr(n, 1, n, Q, n, TAU, WORK, LWORK, INFO)
    if (INFO /= 0) write (6,*) 'SYLVAIN:: fatal error in dorghr ', INFO

    ! set lower triangle of H to zero:
    do i = 1,n - 2
       do j = i + 2, n
          H(j,i) = 0._wp
       end do
    end do

    
    ! 2. Compute the Shur factorization
    ! H = ZTZ'
    allocate( T(n,n), Z(n,n), WR(n), WI(n) )
    T = H
    !Z = Q ! initialization, compz='V' says that.
    call dhseqr('S', 'I', n, 1, n, T, n, WR, WI, Z, n, WORK, LWORK, INFO )
    
    if (INFO /= 0) write (6,*) 'SYLVAIN:: fatal error in dhseqr ', INFO
    
    
    ! 3. Solve sylvester equation:
    ! TX + XT’ = C  
    ! set CZ = Q'Z' C ZQ:
    allocate ( QZ (n,n), CZ(n,n) )
    QZ = matmul(Q,Z)
    CZ = matmul(transpose(QZ), C)
    CZ = matmul(CZ, QZ )
    
    allocate (Xsylv (n,n))
    Xsylv = CZ
    call dtrsyl( 'N', 'T', 1, n, n, T, n, T, n, Xsylv, n, SCALE, INFO )
    
    if (INFO /= 0) write (6,*) 'SYLVAIN:: fatal error in dtrsyl ', INFO


    ! 4. Construct solution:
    X = matmul(QZ, Xsylv)
    X = matmul(X, transpose(QZ) )

  end subroutine lyap


  subroutine f_crps( mu, sig2, w, obs, n, crps)
    use etkpf_util
    implicit none
  
    real(wp) ,intent(in)  :: mu(n)           ! mixture means
    real(wp) ,intent(in)  :: sig2            ! mixture variances
    real(wp) ,intent(in)  :: w(n)            ! mixture weights
    real(wp) ,intent(in)  :: obs             ! observation
    integer  ,intent(in)  :: n               ! number of mixture components
    real(wp) ,intent(out) :: crps            ! computed CRPS
    
    !dummy variables:  
    integer               :: i,j
    real(wp)              :: t1, t2
    
    t1=0._wp
    do i=1,n
      t1 = t1 + w(i) * A_crps( obs - mu(i), sig2)
    end do
    
    t2=0._wp
    do i=1,n
      do j=1,n
        t2 = t2 + w(i) * w(j) * A_crps( mu(i) - mu(j), 2*sig2 )
      end do
    end do
    
    crps = t1 - 0.5 * t2
  end subroutine f_crps


! unequal sig2:
  subroutine f_crps2( mu, sig2, w, obs, n, crps)
    use etkpf_util
    implicit none
  
    real(wp) ,intent(in)  :: mu(n)           ! mixture means
    real(wp) ,intent(in)  :: sig2(n)         ! mixture variances
    real(wp) ,intent(in)  :: w(n)            ! mixture weights
    real(wp) ,intent(in)  :: obs             ! observation
    integer  ,intent(in)  :: n               ! number of mixture components
    real(wp) ,intent(out) :: crps            ! computed CRPS
    
    !dummy variables:  
    integer               :: i,j
    real(wp)              :: t1, t2
    
    t1=0._wp
    do i=1,n
      t1 = t1 + w(i) * A_crps( obs - mu(i), sig2(i))
    end do
    
    t2=0._wp
    do i=1,n
      do j=1,n
        t2 = t2 + w(i) * w(j) * A_crps( mu(i) - mu(j), sig2(i) + sig2(j) )
      end do
    end do
    
    crps = t1 - 0.5 * t2

    
  end subroutine


