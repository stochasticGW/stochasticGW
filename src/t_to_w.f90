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
subroutine ct_to_cw(ct, cw, nt, nw, dt)
  implicit none
  integer  i
  integer nt
  integer nw
  real*8 dt
  complex*16 ct(-nt:nt)
  complex*16 cw(-nw/2:nw/2-1)
  complex*16, allocatable :: c0(:)

  call check_le(nt*2+2, nw, ' nt2nw  ')
  allocate(c0(0:nw-1), stat=i); call check0(i,' nw  ')

  c0 = 0d0
  c0(0:nt)       = ct(0:nt)
  c0(nw-nt:nw-1) = ct(-nt:-1)

  c0             = conjg(c0)
  call fftw_good(c0,nw)
  c0             = conjg(c0)*dt

  cw = 0
  cw(0: nw/2-1)  = c0(0:nw/2-1)
  cw(-nw/2:-1 )  = c0(nw/2:nw-1)

  deallocate(c0)
end subroutine ct_to_cw
