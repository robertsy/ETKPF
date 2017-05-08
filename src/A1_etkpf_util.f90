MODULE etkpf_util

  implicit none
  public
  !---------------------
  ! Real kind parameters
  !---------------------
  !
  !>  8 byte real kind parameter
  !
  integer, parameter :: dp = selected_real_kind(13)
  !
  !>  4 byte real kind parameter
  !
  integer, parameter :: sp = selected_real_kind(6)
  !
  !>  working precision used in the program
  !
  integer, parameter :: wp = dp
  


  contains
  
  !------------------------------------------------------------------------------
  ! subroutine A = evc d(evl) evc'
  subroutine myeigen( a, evc, evl, n, INFO)
    real(wp) ,intent(in)            :: a   (:,:)
    real(wp) ,intent(out)           :: evc     (:,:) ! eigenvectors of A
    real(wp) ,intent(out)           :: evl     (:)   ! eigenvalues  of A
    integer  ,intent(in)            :: n
    integer*4,intent(out)           :: INFO

    integer                         :: LWORK
    real(wp) ,allocatable           :: WORK(:)

    LWORK = n*n !size of WORK (can check the minimum required when succeeded
    allocate ( WORK (LWORK) )
    
    evc = a
    call dsyev( 'V', 'U', n, evc, n, evl, WORK, LWORK, INFO )

  end subroutine myeigen
!   

  subroutine printmat(a,n)
  real(wp) :: a(:,:)
  integer  :: n !nrow
  integer  :: i
  do i=1,n
    write(6,*) a(i,:)
  end do
  end subroutine  

  !Utils taken from mo_matrix in analysis/oo-model
  !------------------------------------------------------------------------------
  function ident (n) result (I)
  integer ,intent(in) :: n
  real(wp)            :: I (n,n)
    integer :: j
    I = 0._wp
    do j=1,n
      I(j,j) = 1._wp
    end do
  end function ident
  !------------------------------------------------------------------------------
  function diag_matrix (d) result (A)
  real(wp) ,intent(in) :: d(:)                ! diagonal
  real(wp)             :: A (size(d),size(d)) ! returned matrix
    integer :: j
    A = 0._wp
    do j=1,size(d)
      A(j,j) = d(j)
    end do
  end function diag_matrix
  !------------------------------------------------------------------------------
  function diag (A) result (y)
  real(wp) ,intent(in) :: A (:,:)
  real(wp)             :: y (size(A,1))
    integer :: i
    do i=1,size(y)
      y(i) = A(i,I)
    end do
  end function diag
  
  
  
  
  function alnorm ( x, upper )
! taken from: https://people.sc.fsu.edu/~jburkardt/f_src/asa066/asa066.html
!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none

  real ( kind = 8 ), parameter :: a1 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a2 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a3 = 5.92885724438D+00
  real ( kind = 8 ) alnorm
  real ( kind = 8 ), parameter :: b1 = -29.8213557807D+00
  real ( kind = 8 ), parameter :: b2 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: c1 = -0.000000038052D+00
  real ( kind = 8 ), parameter :: c2 = 0.000398064794D+00
  real ( kind = 8 ), parameter :: c3 = -0.151679116635D+00
  real ( kind = 8 ), parameter :: c4 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: c5 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: c6 = 3.99019417011D+00
  real ( kind = 8 ), parameter :: con = 1.28D+00
  real ( kind = 8 ), parameter :: d1 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: d2 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: d3 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: d4 = -15.1508972451D+00
  real ( kind = 8 ), parameter :: d5 = 30.789933034D+00
  real ( kind = 8 ), parameter :: ltone = 7.0D+00
  real ( kind = 8 ), parameter :: p = 0.398942280444D+00
  real ( kind = 8 ), parameter :: q = 0.39990348504D+00
  real ( kind = 8 ), parameter :: r = 0.398942280385D+00
  logical up
  logical upper
  real ( kind = 8 ), parameter :: utzero = 18.66D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  up = upper
  z = x

  if ( z < 0.0D+00 ) then
    up = .not. up
    z = - z
  end if

  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0D+00
    else
      alnorm = 1.0D+00
    end if

    return

  end if

  y = 0.5D+00 * z * z

  if ( z <= con ) then

    alnorm = 0.5D+00 - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 &
      / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0D+00 - alnorm
  end if

  return
end
  
  

  function dnorm ( x )
    implicit none
    real(wp) :: dnorm, x
    real(wp) :: pi = 3.14159265359
  
    dnorm = exp( - 0.5 * x**2 ) / sqrt(2*pi)
  
  end
  
  function A_crps ( mu, sig2 )
    implicit none
    real(wp) :: mu, sig2
    real(wp) :: A_crps
    if ( sig2 <= 1.e-10 ) then
      A_crps = abs(mu)
    else
      A_crps = 2 * sqrt(sig2) * dnorm( mu/ sqrt(sig2) ) + &
               mu * ( 2 * alnorm( mu/sqrt(sig2), .FALSE.) - 1._wp )
    end if
  
  end
 

end module etkpf_util