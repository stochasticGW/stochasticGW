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
subroutine debug_mpi(m)
  use simple_mpi, only: rank, nodes, sync_mpi
  implicit none
  integer m,ir,nds
  nds= max(nodes,1)

  call sync_mpi
  do ir=0, nds-1
     if(rank==ir) then
        write(60000+rank,*)' rank, stage ',rank,m
        call flush(60000+rank)
     endif
     call sync_mpi
  end do
end subroutine debug_mpi


subroutine debug_ns(m)
  use simple_mpi, only: rank, nodes, sync_mpi
  use gwm, only : ns
  implicit none
  integer m,ir,nds
  nds= max(nodes,1)

  call sync_mpi
  do ir=0, nds-1
     if(rank==ir) then
        write(60000+rank,*)' rank, ns, stage ',rank,ns, m
        call flush(60000+rank)
     endif
     call sync_mpi
  end do
end subroutine debug_ns
