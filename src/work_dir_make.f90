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
subroutine work_dir_make
  use simple_mpi, only : rank, sync_mpi
  use gwm,        only : work_dir
  implicit none
  logical fex
  call sync_mpi
  work_dir = 'GW_WORK'
  call work_dir_reset()
  if(rank==0) then 
     inquire(file=trim(work_dir),exist=fex)
     if(.not.fex) then
        call system("mkdir "//trim(work_dir)) 
     end if
  end if
  call sync_mpi
end subroutine work_dir_make
  
