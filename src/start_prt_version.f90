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
subroutine prt_version
  use simple_mpi, only : rank
  implicit none
  
  if(rank==0) then
     write(6,*) "**** sGW CODE ****"
     write(6,*) " v1 - Mar/2017"
     write(6,*) ""
  end if

end subroutine prt_version
