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
subroutine  prop_pert(it)
  use simple_mpi, only : rank
  use gwm
  implicit none
  integer it
  call makerho  
  if(it==0.and.flgdyn) call get_vxc_spn(rho, n, nsp, vxc0)
  call makevt(it)
  call propdt(pt,n,ns,map_sp_pt,1)
end subroutine prop_pert

