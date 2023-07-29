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
subroutine hc(p1,chp,vks,n,ms)
  use gwm, only : nx,ny,nz,ekn,dv
  implicit none
  integer     n,ms
  real*8      vks(n)
  complex*16  p1(n),chp(n)

  if(ms/=1) stop "ERROR: Only single spin channel is supported right now!"
  chp = p1
  call fft3d_forward_many(nx, ny, nz, ms, chp)
  chp = chp * ekn
  call fft3d_backward_many(nx, ny, nz, ms, chp)
  chp = chp + vks(:)*p1  
  call addvnl_c(p1,chp,n,ms,1d0)   ! chp=chp+ vnlscaled p1
end subroutine hc

