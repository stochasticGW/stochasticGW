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
subroutine calc_rho_det_ofpt(rho,n, nsp)
  use gwm, only : pt, ns, occ_det, map_sp_det, map_sp_pt, det_tddft
  implicit none
  integer i, n, nsp, sp
  real*8  rho(n, nsp)

  if(size(pt,1)/=n.or.size(pt,2)/=ns) stop ' size_pt_1,2 n,ns '
  if(.not.det_tddft) stop ' ERROR: stochastic tddft in calc_rho_det_of_pt '
  
  rho = 0d0
  do i=1,size(pt,2)
     sp = map_sp_pt(i)
     call check(sp, map_sp_det(i),' map_sp_det ')
     rho(:,sp) = rho(:,sp) + occ_det(i)*abs(pt(:,i))**2   ! occ goes from 0 to 2, includes spin
  enddo

end subroutine calc_rho_det_ofpt

