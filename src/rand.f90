!                     
!                     
!    The routine(s) in this file are a part of the  
!                     StochasticGW                  
!    suite, developed 2012-2018, and copyrighted    
!    to the authors of StochasticGW , see:          
!                                                   
!                 www.stochasticgw.com              
!                                                   
!   If you use or modify any part of this routine   
!   the header should be kept, unmodified.          
!                                                   
!                                                   
!                                                   
subroutine rand_r48(g,n,dv)
  implicit none
  integer               :: n
  integer               :: i,j

  real*4                :: g(n) ! note *4

  real*8            :: dv,r, b
  real*8, external  :: ran_ps

  b = dsqrt(3d0/dv)
  do i=1,20
     r = ran_ps()
  enddo

  do i=1,n;  
     g(i) = b* ( 2d0*ran_ps() - 1d0 ) 
  enddo
end subroutine rand_r48

subroutine rand_r(g,n,dv)
  implicit none
  integer               :: n
  integer               :: i,j
  real*8            :: g(n)
  real*8            :: dv,r, b
  real*8, external  :: ran_ps

  b = dsqrt(3d0/dv)
  do i=1,20
     r = ran_ps()
  enddo

  do i=1,n;  
     g(i) = b* ( 2d0*ran_ps() - 1d0 ) 
  enddo
end subroutine rand_r

subroutine rand_c(g,n,dv)
  implicit none
  integer                  :: n
  integer                  :: i
  complex*16               :: g(n)
  complex*16, parameter    :: ci=(0d0,1d0)
  real*8                   :: dv, b, pi, r
  real*8, external         :: ran_ps

  pi = dacos(-1d0)
  b = 1/dsqrt(dv)

  do i=1,20
     r = ran_ps()
  enddo

  do i=1,n;  
     g(i) = b * exp(2d0*pi*ci*ran_ps())
  enddo
end subroutine rand_c

subroutine rand_i(i,n)
  implicit none
  !
  ! i: random integer between 1 and n
  ! 
  integer i, n
  real*8, external :: ran_ps
  i = ceiling(ran_ps()*dble(n))
end subroutine rand_i
