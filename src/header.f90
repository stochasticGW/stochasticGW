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
subroutine write_header
  use simple_mpi, only : rank, nodes
  implicit none
  if(rank==0) then
     write(6,*)' ************************************************************** '
     write(6,*)
     write(6,*)' Stochastic GW, v3.0 (Oct/2022) '
     write(6,*)
     write(6,*)' ************************************************************** '
     write(6,*)
     write(6,*)' Program uses: ',nodes,' cores '
     call flush(6)
  end if
end subroutine write_header
  
