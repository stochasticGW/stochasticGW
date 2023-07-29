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
subroutine pp_prep
  use simple_mpi, only : rank
  implicit none
  call pp_allocate
  call pp_read
  if(rank==0) then; write(6,*)' Finished reading PP '; write(6,*); call flush(6); endif
  call pp_check
  if(rank==0) then; write(6,*)' Finished checking PP '; call flush(6); endif
  call pp_bcast
end subroutine pp_prep
