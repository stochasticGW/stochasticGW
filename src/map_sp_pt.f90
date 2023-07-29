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
 subroutine prep_map_sp_pt
  use gwm, only : det_tddft
  use gwm, only : map_sp_pt, map_sp_det
  use gwm, only : nsp, ns, is_map, ns_blk
  use simple_mpi, only : color_rank
  implicit none
  integer st, is, is_full

  allocate(map_sp_pt(ns), stat=st)
  call check0(st,' map_sp  ')

  if(nsp==1) then
     map_sp_pt = 1
     return
  endif

  if(det_tddft) then ! det. tddft
     map_sp_pt = map_sp_det
  else
     do is=1,ns
        is_full = is_map(is,color_rank)
        call check_lele(1,is_full,ns_blk,' one, is_full ns_blk ')
        map_sp_pt(is) = 1 + mod(is_full-1,nsp)  ! 1 2 3 4 5 ..--> 1 2 1 2 1
     enddo
  end if
end subroutine prep_map_sp_pt
