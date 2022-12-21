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
subroutine check_input_name
  use gwm, only : inputfname
  implicit none
  if(len(trim(inputfname)) .lt. 1) then 
     write(6,*) "    -- Input file not specified: stopping"
     stop
  endif
end subroutine check_input_name
