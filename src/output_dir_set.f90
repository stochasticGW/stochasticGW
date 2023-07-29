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
subroutine output_dir_set
  use gwm, only : output_dir, icounter
  implicit none
  character*20 fr
  call cnvrt_number_to_string_20(icounter,fr)
  output_dir='GW_OUTPUT.'//trim(adjustl(fr))
end subroutine output_dir_set
