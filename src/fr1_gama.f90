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
subroutine fr1_gama(i,it)
  use gwm
  implicit none
  integer i, it
  real*8  wg
  wg = wgf(it)
  select case(i)
  case(0); ft(it,:) =             gamvr(:,sp0)
  case(1); ft(it,:) = wg/sm*(dble(gamvr(:,sp0)) - dble(ft(it,:))) ! for spin-dep extend
  case default; stop ' pfr '
  end select
end subroutine fr1_gama
