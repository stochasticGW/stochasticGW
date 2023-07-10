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
subroutine read_orb
  use gwm
  implicit none

  if (.not.use_unified_input) then
     if (orb_indx==-1.or.(.not.det_dft)) then
        call read_orb_1
     else
        call read_orb_fromwf
     endif

     if (nj>0) then
        if (read_jidx_input) then
           call read_orbj_fromwf
        else
           call read_orbj_1
        endif
     endif
  endif

end subroutine read_orb
