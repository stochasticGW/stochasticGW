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
  use gwm, only : rddh, disable_gpu_hminmax
#if GPU_ENABLED
  use device_mem_module
#endif

  implicit none

  if (.not.rddh) then
#if GPU_ENABLED
    if (disable_gpu_hminmax) then
       call h_minmax
    else
       call init_device
       call hmin_hmax_gpu
       call flush_device
    endif
#else
    call h_minmax
#endif
  endif

  call proc_hminmax_info

end subroutine improve_hmin_hmax
