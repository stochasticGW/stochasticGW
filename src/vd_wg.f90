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
subroutine vd_wg(mb,mt,vr0,vr1,vd)
  use gwm
  implicit none
  integer mb,  mt, it
  real*4      vr0(mb:mt,0:nt)
  real*4      vr1(mb:mt,0:nt)
  complex*16  vd(mb:mt,0:nt)

  do it=0,nt
     vd(:,it) = wgf(it)*(dble(vr1(:,it))-dble(vr0(:,it)))/sm
  enddo
end subroutine vd_wg
  
  
