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
subroutine prnt_time_mpi(a,tstart)
  use simple_mpi
  implicit none
  character(*) a
  real*8 tstart, tend

  call sync_mpi
  if(rank==0) then
    call cpu_time(tend); write(6,*) " Time: ",a,real(tend-tstart)
  end if
  call sync_mpi

end subroutine prnt_time_mpi

subroutine prnt_time_i_mpi(a,i,tstart)
  use simple_mpi, only : rank, sync_mpi
  implicit none
  character(*) a
  integer i
  real*8 tstart, tend

  call sync_mpi
  if(rank==0) then
     call cpu_time(tend)
     write( 6,*) " Time: ",a,i,"=> time: ",real(tend-tstart)
     write(17,*) " Time: ",a,i,"=> time: ",real(tend-tstart)
     call flush(6)
     call flush(17)
  end if
  call sync_mpi

end subroutine prnt_time_i_mpi
