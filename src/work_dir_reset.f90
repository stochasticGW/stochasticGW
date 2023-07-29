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
subroutine work_dir_reset
  use gwm, only : work_dir, icounter
  implicit none
  character*20 fr
  call cnvrt_number_to_string_20(icounter,fr)
  work_dir=trim(adjustl(work_dir))//'.'//trim(adjustl(fr))
  call work_dir_cleanup
end subroutine work_dir_reset

subroutine work_dir_cleanup
  use gwm, only : work_dir, icounter
  use simple_mpi, only : rank, sync_mpi
  implicit none
  if(rank==0) then
     call system("touch "//trim(work_dir)//"/tmp.txt")  ! ensures there's something in work_dir
     call system("rm "//trim(work_dir)//"/*")
  endif
  call sync_mpi
end subroutine work_dir_cleanup
