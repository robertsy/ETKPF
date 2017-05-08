!
!+ Fortran 90 wrapper for LAPACK singular value decomposition routines
!
! $Id$
!
!------------------------------------------------------------------------------
! #if defined (__SX__)
! ! This compiler directive must appear before the MODULE statement!
! !OPTION! -pvctl matmul
! #endif
!------------------------------------------------------------------------------
!
MODULE mo_matrix
!
! Description:
!   Fortran 90 wrapper for LAPACK singular value decomposition routines
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_8         2009/12/09 Harald Anlauf
!  optimization for NEC SX-9
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  new function diag_matrix: returns diagonal matrix
! V1_29        2014/04/02 Andreas Rhodin
!  check_rs: optionally return eigenvectors and eigenvalues
! V1_45        2015-12-15 Harald Anlauf
!  check_rs: work around crayftn performance problem: matmul(a,transpose(b))
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2002-2008
!================================================
!
! Diagnostics for vectorization
!
! #if defined (_FTRACE) && 0
! #define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
! #define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
! #else
! #define FTRACE_BEGIN(text)
! #define FTRACE_END(text)
! #endif
!============================================================================
  !-------------
  ! Modules used
  !-------------
  use etkpf_util,               only: wp       ! working precision
! use full_hessian_type_module, only: inverse  ! matrix inversion (SVD)
! use slatec_module,            only: svdc     ! singular value dec. (slatec)
! use mo_algorithms,            only: svdcmp, &! singular value dec. (wrapper)
!   use mo_exception,             only: finish
  implicit none
  !----------------
  ! Public entities
  !----------------
  private
!   public :: inverse        ! matrix inversion (singular value decomposition)
  public :: eigen       ! compute eigen decomposition of A
  public :: inverse_rs     ! matrix inversion (real symmetric)
  public :: sqrt_rs        ! square root (real symmetric)
  public :: pow_rs         ! A ** x (real symmetric)
  public :: check_rs       ! check eigenvalues of real symmetric matrix
  public :: cholesky       ! Cholesky factorization
  public :: solve_cholesky ! solve  A * x = b ,Cholesky factorized A
  public :: inv_cholesky   ! invert A         ,Cholesky factorized A
  public :: pack           ! convert matrix to   packed representation
  public :: unpack         ! convert matrix to unpacked representation
  public :: operator(.x.)  ! matrix-vector multiplication
  public :: operator(.o.)  ! outer product
  public :: operator(.i.)  ! inner product
!   public :: scale          ! scale a matrix
  public :: diag           ! returns diagonal
  public :: diag_matrix    ! returns diagonal matrix
  public :: ident          ! identity matrix
  !-----------
  ! Interfaces
  !-----------
  interface operator(.x.)
    module procedure matrix_times_vector
    module procedure vector_times_matrix
  end interface operator(.x.)

  interface operator(.o.)
    module procedure outer_product
  end interface operator(.o.)

  interface operator(.i.)
    module procedure matrix_times_matrix
  end interface operator(.i.)

!   interface inverse
!     module procedure inverse_svd
!     module procedure inverse_submatrix
!   end interface inverse

  interface cholesky
    module procedure cholesky_ap ! Cholesky factoriz., packed representation
    module procedure cholesky_a  ! Cholesky factoriz., full matrix repr.
  end interface cholesky

  interface pack
    module procedure pack_matrix
  end interface pack

  interface unpack
    module procedure unpack_matrix
  end interface unpack

!==================
! Module procedures
!==================
contains
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
!------------------------------------------------------------------------------
  !----------------
  ! matrix * matrix
  !----------------
  function matrix_times_matrix (A, B) result (Y)
  real(wp) ,intent(in) :: A (:,:)
  real(wp) ,intent(in) :: B (:,:)
  real(wp)             :: Y (size(A,1),size(B,2))
    Y = matmul (A, B)
  end function matrix_times_matrix
!------------------------------------------------------------------------------
  !----------------
  ! matrix * vector
  !----------------
  function matrix_times_vector (A, x) result (y)
  real(wp) ,intent(in) :: A (:,:)
  real(wp) ,intent(in) :: x (:)
  real(wp)             :: y (size(A,1))
    y = matmul (A, x)
  end function matrix_times_vector
!------------------------------------------------------------------------------
  !----------------
  ! vector * matrix
  !----------------
  function vector_times_matrix (x, A) result (y)
  real(wp) ,intent(in) :: x (:)
  real(wp) ,intent(in) :: A (:,:)
  real(wp)             :: y (size(A,2))
    y = matmul (x, A)
  end  function vector_times_matrix
!------------------------------------------------------------------------------
  !---------------------------
  ! matrix = vector .o. vector
  !---------------------------
  function outer_product (x1, x2) result (Y)
  real(wp) ,intent(in) :: x1 (:)
  real(wp) ,intent(in) :: x2 (:)
  real(wp)             :: Y  (size(x1),size(x2))
    integer :: i,j
    do j=1,size(x2)
      do i=1,size(x1)
        Y(i,j) = x1(i) * x2(j)
      end do
    end do
  end function outer_product

!------------------------------------------------------------------------------
  function inverse_rs (a, max_ev, min_ev) result (y)
  real(wp), intent(in)           :: a(:,:)
  real(wp), intent(in), optional :: max_ev
  real(wp), intent(in), optional :: min_ev
  real(wp)                       :: y (size(a,1),size(a,1))
  !--------------------------------------
  ! Inverts the (real symmetric) matrix a
  !--------------------------------------
    call check_rs (a, max_ev, min_ev, inv=y)
  end function inverse_rs
!------------------------------------------------------------------------------
  function sqrt_rs (a, min_ev) result (y)
  real(wp), intent(in)           :: a(:,:)
  real(wp), intent(in), optional :: min_ev
  real(wp)                       :: y (size(a,1),size(a,1))
  !---------------------------------------------
  ! Square root of the (real symmetric) matrix a
  !---------------------------------------------
    call check_rs (a, sqroot=y, min_ev=min_ev)
  end function sqrt_rs
!------------------------------------------------------------------------------
  function pow_rs (a, x, min_ev) result (y)
  real(wp), intent(in)           :: a(:,:)
  real(wp), intent(in)           :: x
  real(wp), intent(in), optional :: min_ev
  real(wp)                       :: y (size(a,1),size(a,1))
  !---------------------------------------
  ! A ** x for a (real symmetric) matrix A
  !---------------------------------------
    call check_rs (a, a_pow_x=y, x=x, min_ev=min_ev)
  end function pow_rs
  
!------------------------------------------------------------------------------
! subroutine A = evc d(evl) evc'
  subroutine eigen( a, evc, evl)
    real(wp) ,intent(in)            :: a   (:,:)
    real(wp) ,intent(out)           :: evc     (:,:) ! eigenvectors of A
    real(wp) ,intent(out)           :: evl     (:)   ! eigenvalues  of A
    
    call check_rs (a, evc=evc, evl=evl)
  end subroutine eigen
!------------------------------------------------------------------------------
  !--------------------------------
  ! inverse (real symmetric) matrix
  !--------------------------------
  subroutine check_rs (a, max_ev, min_ev, y, inv, hits, evmin, evmax, sqroot,&
                       a_pow_x, evc, evl, x)
  real(wp) ,intent(in)            :: a   (:,:)
  real(wp) ,intent(in)  ,optional :: max_ev
  real(wp) ,intent(in)  ,optional :: min_ev
  real(wp) ,intent(out) ,optional :: y   (:,:)
  real(wp) ,intent(out) ,optional :: inv (:,:)     ! inverse (A)
  integer  ,intent(out) ,optional :: hits
  real(wp) ,intent(out) ,optional :: evmin
  real(wp) ,intent(out) ,optional :: evmax
  real(wp) ,intent(out) ,optional :: sqroot  (:,:) ! sqrt (A)
  real(wp) ,intent(out) ,optional :: a_pow_x (:,:) ! A ** x
  real(wp) ,intent(out) ,optional :: evc     (:,:) ! eigenvectors of A
  real(wp) ,intent(out) ,optional :: evl     (:)   ! eigenvalues  of A
  real(wp) ,intent(in)  ,optional :: x
  !----------------------------------------------------
  ! checks eigenvalues of the (real symmetric) matrix a
  !----------------------------------------------------
    real(wp)              :: z    (size(a,1),size(a,1))
    real(wp)              :: zw   (size(a,1),size(a,1))
    real(wp)              :: w    (size(a,1))
    real(wp) ,allocatable :: work (:)
    integer               :: i, n, lwork, info, nb
    integer, external     :: ilaenv
! #if defined (__SX__) || defined (_CRAYFTN)
!     integer               :: j, k
! #endif

!     if (size(a,1)/=size(a,2)) call finish ('invert_rs','matrix is not square')
    if (size(a,1)==0) then
      if(present (hits)) hits   =  0
      if(present (evmin)) evmin =  huge(evmin)
      if(present (evmax)) evmax = -huge(evmax)
      return
    endif
    !-------------------------------------------------------
    ! call the Lapack driver routine (eigenvalues + vectors)
    !-------------------------------------------------------
    z     = a
    n     = size(a,1)
    nb    = ilaenv( 1, 'DSYTRD', 'U', n, -1, -1, -1 )
    lwork = (nb+2)*n
    allocate (work(lwork))
! FTRACE_BEGIN("check_rs:dsyev")
    call dsyev('V', 'U', n, z, n, w, work, lwork, info )
! FTRACE_END  ("check_rs:dsyev")
    deallocate (work)
    !-----------------
    ! check for errors
    !-----------------
    if (info<0) then
      write (0,*) 'inverse_rs: dsyev: ',info,      &
                  'th argument had an illegal value'
!       call finish('inverse_rs','dsyev: argument had an illegal value')
    else if (info>0) then
!       call finish('inverse_rs','dsyev: the algorithm failed to converge')
    endif
    !----------------------------------------------------------------
    ! bounds on eigenvalues, report if optional arguments are present
    !----------------------------------------------------------------
    if (present (hits)) then
      hits = 0
      if (present (min_ev)) hits = count (w < min_ev)
      if (present (max_ev)) hits = count (w > max_ev) + hits
    endif
    if (present (evmin)) evmin = minval (w)
    if (present (evmax)) evmax = maxval (w)
    !------------------------------
    ! bounds on eigenvalues, modify
    !------------------------------
    if (present (min_ev)) w = max (min_ev, w)
    if (present (max_ev)) w = min (max_ev, w)
    !--------------------------
    ! calculate modified matrix
    !--------------------------
    if (present(y)) then
! FTRACE_BEGIN("check_rs:mod")
      do i = 1, size (w)
        zw(:,i) = z(:,i) * w(i)
      end do
!     y = matmul ( zw, transpose (z)) !+++ segm.fault with NAG -mtrace
      call dgemm ('n','t', n, n, n, 1._wp, zw, n, z, n, 0._wp, &
                  y, size (y,dim=1))
! FTRACE_END  ("check_rs:mod")
    endif
    !------------------
    ! calculate inverse
    !------------------
    if (present(inv)) then
! FTRACE_BEGIN("check_rs:inv")
      do i = 1, size (w)
        zw(:,i) = z(:,i) * (1._wp / w(i))
      end do
! #if defined (__SX__) || defined (_CRAYFTN)
!       ! Work around poor performance on SX-9 and with Crayftn
!       if (n >= 4) then
! !        inv = matmul ( zw, transpose (z))
!          call dgemm ('n','t', n, n, n, 1._wp, zw, n, z, n, 0._wp, &
!                      inv, size (inv,dim=1))
!       else
!          ! Explicit coding for n <= 3
!          do j = 1, n
!             do i = 1, n
!                inv(i,j) = 0._wp
! !CDIR NOVECTOR LOOPCNT=3
!                do k = 1, n
!                   inv(i,j) = inv(i,j) + zw(i,k) * z(j,k)
!                end do
!             end do
!          end do
!       end if
! #else
      inv = matmul ( zw, transpose (z))
! #endif
! FTRACE_END  ("check_rs:inv")
    endif
    !----------------------
    ! calculate square root
    !----------------------
    if (present(sqroot)) then
! FTRACE_BEGIN("check_rs:sqroot")
      do i = 1, size (w)
        zw(:,i) = z(:,i) * sqrt (w(i))
      end do
      sqroot = matmul ( zw, transpose (z))
! FTRACE_END  ("check_rs:sqroot")
    endif
    !-----------------
    ! calculate a ** x
    !-----------------
    if (present(a_pow_x)) then
! FTRACE_BEGIN("check_rs:pow")
      do i = 1, size (w)
        zw(:,i) = z(:,i) * (w(i) ** x)
      end do
      a_pow_x = matmul ( zw, transpose (z))
! FTRACE_END  ("check_rs:pow")
    endif
    !------------------------------------
    ! return eigenvectors and eigenvalues
    !------------------------------------
    if (present (evc)) evc = z
    if (present (evl)) evl = w

  end subroutine check_rs

!------------------------------------------------------------------------------
 
!==============================================================================
  !-------------------------------
  ! Perform Cholesky factorization
  !-------------------------------
  subroutine cholesky_ap (ap, m, uplo, thr, info)
  real(wp) ,intent(inout)         :: ap (:) ! matrix in packed storage format
  integer  ,intent(in)            :: m      ! dimension
  character,intent(in)  ,optional :: uplo   ! 'U' or 'L', default: 'U'
  real(wp) ,intent(in)  ,optional :: thr    ! threshold
  integer  ,intent(out) ,optional :: info   ! error flag

    character :: ul
    integer   :: inf
    real(wp)  :: th

    ul = 'U';    if (present(uplo)) ul = uplo
    th = -1._wp; if (present(thr))  th = thr

    if (th < 0._wp) then
      call dpptrf    ( ul, m, ap, inf )
    else
!       if (ul /= 'U') call finish ('cholesky_ap',&
!                      'thr present, works only with uplo == U')
      if (th/=0) print *,'cholesky_ap: thr currently taken as 0.'
      call dpptrf_my ( ul, m, ap, inf)
    endif

    if (inf/=0) then
      write(0,*) 'cholesky_ap: info=',inf
      if(inf<0) write(0,*) 'the',-inf,'-th argument had an illegal value.'
      if(inf>0) then
        write(0,*)'the leading minor of order ',inf,'is not positive definite'
        write(0,*)'and the factorization could not be completed.'
      endif
!       if (.not.present(info)) call finish ('cholesky_ap','info/=0')
    endif
    if (present(info)) info = inf

  end subroutine cholesky_ap
!------------------------------------------------------------------------------
  subroutine cholesky_a (a, uplo, thr, info)
  real(wp) ,intent(inout)         :: a (:,:) ! matrix to factorize
  character,intent(in)  ,optional :: uplo    ! 'U' or 'L', default: 'U'
  real(wp) ,intent(in)  ,optional :: thr     ! threshold
  integer  ,intent(out) ,optional :: info    ! error flag

    real(wp)  :: ap (size(a,1)*(size(a,1)+1)/2)
!   character :: ul
    integer   :: m

    m  = size(a,1)
!   ul = 'U'; if (present(uplo)) ul = uplo

    call pack_matrix   (a, ap, uplo)
    call cholesky_ap   (ap, m, uplo, thr=thr, info=info)
    call unpack_matrix (ap, a, uplo)

  end subroutine cholesky_a
!------------------------------------------------------------------------------
  function solve_cholesky (ap, b, uplo, info) result (x)
  !---------------------------------------------------
  ! solve a * x = b using a cholesky factorized matrix
  !---------------------------------------------------
  real(wp), intent(in)            :: ap (:)     ! matrix, packed representation
  real(wp), intent(in)            :: b  (:)     ! rhs
  character,intent(in)  ,optional :: uplo       ! 'U' or 'L', default: 'U'
  integer,  intent(out) ,optional :: info       ! error message
  real(wp)                        :: x(size(b)) ! solution

    character :: ul
    integer   :: m
    integer   :: inf

    m  = size (b)
    ul = 'U'; if (present(uplo)) ul = uplo

    x = b
    call dpptrs (ul, m, 1, ap , x, m, inf)

    if (inf/=0) then
      write(0,*) 'solve_cholesky: info=',inf
!       if (.not.present(info)) call finish ('solve_cholesky','info/=0')
    endif
    if (present(info)) info = inf

  end function solve_cholesky
!------------------------------------------------------------------------------
  subroutine inv_cholesky (ap, uplo, info)
  !---------------------------------------------------
  ! derive the inverse of a cholesky factorized matrix
  !---------------------------------------------------
  real(wp), intent(inout)         :: ap (:)     ! matrix, packed representation
  character,intent(in)  ,optional :: uplo       ! 'U' or 'L', default: 'U'
  integer,  intent(out) ,optional :: info       ! error message

    character :: ul
    integer   :: m
    integer   :: inf

    m  = sqrt (float(size(ap)*2))
    ul = 'U'; if (present(uplo)) ul = uplo

    call dpptri(ul, m, ap, inf)

    if (inf/=0) then
      write(0,*) 'inv_cholesky: info=',inf
!       if (.not.present(info)) call finish ('inv_cholesky','info/=0')
    endif
    if (present(info)) info = inf

  end subroutine inv_cholesky
!------------------------------------------------------------------------------
  subroutine pack_matrix (x, y, uplo)
  real(wp) ,intent(in)            :: x (:,:) ! matrix to store in packed format
  real(wp) ,intent(out)           :: y (:)   ! matrix stored   in packed format
  character,intent(in)  ,optional :: uplo    ! 'U' or 'L', default: 'U'

    character :: ul
    integer   :: i, j, k
    integer   :: n

    ul = 'U'; if (present(uplo)) ul = uplo
    n = size (x,1)
    k = 0
    select case (ul)
    case ('U')
      do j = 1, n
        do i = 1, j
          k=k+1
          y (k) = x (i,j)
        end do
      end do
    case ('L')
      do j = 1, n
        do i = j, n
          k=k+1
          y (k) = x (i,j)
        end do
      end do
    case default
!       call finish ('pack_block','invalid triangle: '//ul)
    end select
  end subroutine pack_matrix
!------------------------------------------------------------------------------
  subroutine unpack_matrix (x, y, uplo, lsym)
  real(wp) ,intent(out)           :: y (:,:) ! matrix to store in full format
  real(wp) ,intent(in)            :: x (:)   ! matrix stored   in packed format
  character,intent(in)  ,optional :: uplo    ! 'U' or 'L', default: 'U'
  logical  ,intent(in)  ,optional :: lsym    ! unp.both triangles of sym.matrix

    character :: ul
    integer   :: i, j, k
    integer   :: n
    logical   :: l

    l  = .false. ; if (present(lsym))  l = lsym
    ul = 'U'     ; if (present(uplo)) ul = uplo
    n  = size (y,1)
    k  = 0
    if (.not.l) y = 0._wp
    select case (ul)
    case ('U')
      do j = 1, n
        do i = 1, j
          k=k+1
          y (i,j) = x (k)
  if(l)   y (j,i) = x (k)
        end do
      end do
    case ('L')
      do j = 1, n
        do i = j, n
          k=k+1
          y (i,j) = x (k)
  if(l)   y (j,i) = x (k)
        end do
      end do
    case default
!       call finish ('unpack_block','invalid triangle: '//ul)
    end select
  end subroutine unpack_matrix
!------------------------------------------------------------------------------
      SUBROUTINE DPPTRF_my( UPLO, N, AP, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * )
!     ..
!
!  Purpose
!  =======
!
!  DPPTRF computes the Cholesky factorization of a real symmetric
!  positive definite matrix A stored in packed format.
!
!  The factorization has the form
!     A = U**T * U,  if UPLO = 'U', or
!     A = L  * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the symmetric matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!          See below for further details.
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U**T*U or A = L*L**T, in the same
!          storage format as A.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  Further Details
!  ======= =======
!
!  The packed storage scheme is illustrated by the following example
!  when N = 4, UPLO = 'U':
!
!  Two-dimensional storage of the symmetric matrix A:
!
!     a11 a12 a13 a14
!         a22 a23 a24
!             a33 a34     (aij = aji)
!                 a44
!
!  Packed storage of the upper triangle of A:
!
!  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JC, JJ
      DOUBLE PRECISION   AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, DSPR, DTPSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPPTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
      IF( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U'*U.
!
         JJ = 0
         DO 10 J = 1, N
            JC = JJ + 1
            JJ = JJ + J
!
!           Compute elements 1:J-1 of column J.
!
            IF( J.GT.1 )                                             &
               CALL DTPSV( 'Upper', 'Transpose', 'Non-unit', J-1, AP,&
                           AP( JC ), 1 )
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            AJJ = AP( JJ ) - DDOT( J-1, AP( JC ), 1, AP( JC ), 1 )
            IF( AJJ.LE.ZERO ) THEN
!!!               AP( JJ ) = AJJ
!!!               GO TO 30
              AJJ = AP( JJ )
              AP (JC+1:JC+J-1) = 0.
            END IF
            AP( JJ ) = SQRT( AJJ )
   10    CONTINUE
      ELSE
!
!        Compute the Cholesky factorization A = L*L'.
!
         JJ = 1
         DO 20 J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            AJJ = AP( JJ )
            IF( AJJ.LE.ZERO ) THEN
               AP( JJ ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            AP( JJ ) = AJJ
!
!           Compute elements J+1:N of column J and update the trailing
!           submatrix.
!
            IF( J.LT.N ) THEN
               CALL DSCAL( N-J, ONE / AJJ, AP( JJ+1 ), 1 )
               CALL DSPR( 'Lower', N-J, -ONE, AP( JJ+1 ), 1, &
                          AP( JJ+N-J+1 ) )
               JJ = JJ + N - J + 1
            END IF
   20    CONTINUE
      END IF
      GO TO 40
!
   30 CONTINUE
      INFO = J
!
   40 CONTINUE
      RETURN
!
!     End of DPPTRF
!
      END SUBROUTINE DPPTRF_my
!==============================================================================
end module mo_matrix
