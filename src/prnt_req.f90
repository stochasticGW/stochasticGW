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
subroutine prnt_mem_or_scratch_req
  use gwm,        only : gamflg, n, nsp, flgdyn, ngam, nt
  use simple_mpi, only : rank, sync_mpi, nodes
  implicit none
  integer nspv
  real*8  memgb, wrttb

  nspv = 1
  if(flgdyn) nspv = nsp 

  call sync_mpi  
  if(rank==0) then
     write(6,*)
     select case(gamflg)
     case(1,2)    
        memgb = dble(nt+1)* dble(ngam) * 8d0/1e9
        write(6,*)' NOTE: memory requirements per core due to psi -- in GB ',real(memgb)
     case(3)
        wrttb = dble(nt+1)* dble(n) * dble(nspv)* dble(max(1,nodes)) *8d0/1e12
        write(6,*)' NOTE: scratch-disk requirements, total, due to storing vt(r,t) -- in TB ',real(wrttb)
     case default
        stop ' ERROR: gamflg wrong '
     end select
     write(6,*)
     call flush(6)
  endif
  call sync_mpi
end subroutine prnt_mem_or_scratch_req
