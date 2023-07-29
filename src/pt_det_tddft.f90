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
subroutine pt_det_tddft
  use gwm
  use simple_mpi, only : rank
  implicit none
  integer j
  if(ns>0) then
     call check(size(pt,1),size(phi_det,1),' sz_pt_wf_1 ')
     call check(size(pt,2),size(phi_det,2),' sz_pt_wf_2 ')
     pt = phi_det
     normpt = 1d0  ! occupation in occ_det
  endif
end subroutine pt_det_tddft
