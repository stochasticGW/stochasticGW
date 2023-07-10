! ----------------------------------------------------------------------
! The file below is a  FORTRAN90 module to generate random numbers
! with uniform probability distribution
! Written June 2011, by Jean-Michel Brankart; 
! Updated address (2018): CNRS - Université Grenoble Alpes, IGE, Grenoble, France
! 
! A few tiny changes were inserted by DN for an easier mesh with stochasticGW; in addition,
! many of the original subroutine in this file were removed since only a uniform distribution is
! needed in stochasticGW.
!
! ----------------------------------------------------------------------
! The module is based on (and includes) the
! 64-bit KISS (Keep It Simple Stupid) random number generator
! distributed by George Marsaglia :
!
! http://groups.google.com/group/comp.lang.fortran/
!        browse_thread/thread/a85bf5f2a97f5a55
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! kiss          : 64-bit KISS random number generator (period ~ 2^250)
! kiss_seed     : Define seeds for KISS random number generator
! kiss_save     : Save current state of KISS (for future restart)
! kiss_load     : Load the saved state of KISS
! kiss_reset    : Reset the default seeds
! kiss_uniform  : Real random numbers with uniform distribution in [0,1]
! kiss_sample   : Select a random sample from a set of integers
! ----------------------------------------------------------------------
MODULE mod_kiss
  IMPLICIT NONE
  PRIVATE
! Public functions/subroutines
  PUBLIC :: kiss, kiss_seed,  kiss_seed_extract_DN
  PUBLIC :: kiss_uniform, kiss_sample
! Default/initial seeds
!  INTEGER(KIND=8), parameter :: x0=1234567890987654321_8  ! DN -- modified x to x0
!  INTEGER(KIND=8), parameter :: y0=362436362436362436_8
!  INTEGER(KIND=8), parameter :: z0=1066149217761810_8
!  INTEGER(KIND=8), parameter  :: w0=123456123456123456_8
  integer(kind=8), save :: x,y,z,w

! Parameters to generate real random variates
  REAL(KIND=8), PARAMETER :: huge64=9223372036854775808.0_8  ! +1
  REAL(KIND=8), PARAMETER :: zero=0.0_8, half=0.5_8, one=1.0_8, two=2.0_8
CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  FUNCTION kiss()
! --------------------------------------------------------------------
! The 64-bit KISS (Keep It Simple Stupid) random number generator.
! Components:
!  (1) Xorshift (XSH), period 2^64-1,
!  (2) Multiply-with-carry (MWC), period (2^121+2^63-1)
!  (3) Congruential generator (CNG), period 2^64.
! Overall period:
!  (2^250+2^192+2^64-2^186-2^129)/6 ~= 2^(247.42) or 10^(74.48)
! Set your own seeds with statement "CALL kiss_seed(ix,iy,iz,iw)".
! --------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(KIND=8) :: kiss, t

  t = ishft(x,58) + w
  IF (s(x).eq.s(t)) THEN
    w = ishft(x,-6) + s(x)
  ELSE
    w = ishft(x,-6) + 1 - s(x+t)
  ENDIF
  x = t + x
  y = m( m( m(y,13_8), -17_8 ), 43_8 )
  z = 6906969069_8 * z + 1234567_8

  kiss = x + y + z

  CONTAINS

    FUNCTION s(k)
      INTEGER(KIND=8) :: s, k
      s = ishft(k,-63)
    END FUNCTION s

    FUNCTION m(k, n)
      INTEGER(KIND=8) :: m, k, n
      m = ieor (k, ishft(k, n) )
    END FUNCTION m

  END FUNCTION kiss
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SUBROUTINE kiss_seed(ix, iy, iz, iw)
! --------------------------------------------------------------------
! Define seeds for KISS random number generator ! changed DN
! --------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(KIND=8) :: ix, iy, iz, iw

  x = ix
  y = iy
  z = iz
  w = iw
  
  END SUBROUTINE kiss_seed

  subroutine kiss_seed_extract_DN(ix,iy,iz,iw)
    implicit none
    integer(kind=8) ix,iy,iz,iw
    ix=x
    iy=y
    iz=z
    iw=w
  end subroutine kiss_seed_extract_DN


! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SUBROUTINE kiss_uniform(uran)
! --------------------------------------------------------------------
! Real random numbers with uniform distribution in [0,1]
! --------------------------------------------------------------------
  IMPLICIT NONE
  REAL(KIND=8) :: uran

  uran = half * ( one + REAL(kiss(),8) / huge64 )
  
  !write(77345,*)uran,x,y,z,w;call flush(77345)  !erase 
  END SUBROUTINE kiss_uniform
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SUBROUTINE kiss_sample(a,n,k)
!---------------------------------------------------------------------
! Select a random sample of size k from a set of n integers
! The sample is output in the first k elements of a
! Set k equal to n to obtain a random permutation
!   of the whole set of integers
!---------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(KIND=8), DIMENSION(:) :: a
  INTEGER(KIND=8) :: n, k, i, j, atmp
  REAL(KIND=8) :: uran

! Select the sample using the swapping method
! (see Devroye, Non-Uniform Random Variate Generation, p. 612)
  DO i=1,k
! Randomly select the swapping element between i and n (inclusive)
    CALL kiss_uniform(uran)
!   j = i + INT( REAL(n-i+1,8) * uran )
    j = i - 1 + CEILING( REAL(n-i+1,8) * uran )
! Swap elements i and j
    atmp = a(i) ; a(i) = a(j) ; a(j) = atmp
  ENDDO

  END SUBROUTINE kiss_sample
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE mod_kiss

