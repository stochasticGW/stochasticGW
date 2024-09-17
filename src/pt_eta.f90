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
subroutine pt_eta
  use gwm
  implicit none

  if(det_tddft) then
     call  pt_det_tddft
     call eta_det_tddft
     return
  end if

  call pt_stoch  ! prepare stochastic pt, to be projected.

  if(det_dft) then
     call pt_eta_effic
  else
#if GPU_ENABLED
     if (disable_gpu_filter) then
        call eta_fltr
        call pt_fltr
     else
        call pt_eta_fltr_gpu
     endif
#else
     call eta_fltr
     call pt_fltr
#endif
  endif
  
end subroutine pt_eta

