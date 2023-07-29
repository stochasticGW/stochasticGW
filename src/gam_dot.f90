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
subroutine gam_dot( vtd, n, gam1, seg, j, dv, dot)
  use gwm, only : seg_sh, seg_w
  implicit none
  integer n, j, seg, i, sh
  real*8 vtd(n), dv, dot
  real*4 gam1(seg) 
  
  sh=seg_sh(j)
  dot = 0d0
  do i=1,seg_w(j)
     dot = dot + gam1(i)*vtd(i+sh)
  enddo
  dot = dot *dv
end subroutine gam_dot
