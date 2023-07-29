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
subroutine scratch_dir_reset
  use gwm, only : scratch_dir, icounter
  implicit none
  character*20 fr
  call cnvrt_number_to_string_20(icounter,fr)
  scratch_dir=trim(adjustl(scratch_dir))//'.'//trim(adjustl(fr))
end subroutine scratch_dir_reset
