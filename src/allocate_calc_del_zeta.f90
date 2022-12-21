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
subroutine allocate_calc_del_zeta
  use gwm
  use simple_mpi, only : color_bcast_r8, rank
  implicit none
  integer i
  allocate(del(n),stat=i);call check0(i,' dendel ')

  if(eta_flg) then
     allocate(zeta(n),stat=i);call check0(i,' zeta ')
     call rand_r(zeta,n,dv);
     allocate(den(n),stat=i);call check0(i,' den_zg ')
     den = zeta*g
     call vh_sub(den, del, n)
     deallocate(den)
  endif

  call color_bcast_r8(del,  size(del),  eta_rank)
  
end subroutine allocate_calc_del_zeta
