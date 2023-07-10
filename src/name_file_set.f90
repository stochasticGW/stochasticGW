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
subroutine set_details_file_name
  use simple_mpi, only : rank
  use gwm, only : work_dir
  implicit none
  character(200) filename
  if(rank/=0) return
  filename=trim(adjustl(work_dir))//'/details_output.txt'
  write(6,*)' Details on the results are in : ',trim(adjustl(filename))
  call flush(6)
  open(17,file=trim(adjustl(filename)),status='replace')
end subroutine set_details_file_name
  
