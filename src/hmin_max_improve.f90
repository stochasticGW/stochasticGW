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
subroutine improve_hmin_hmax
  use gwm, only : rddh, usegpu
  implicit none
  if (.not.rddh) then
#if GPU_ENABLED
    if (.not.usegpu) then

    else
!       call hmin_hmax_gpu
    endif
#endif
    call h_minmax
  endif

end subroutine improve_hmin_hmax
