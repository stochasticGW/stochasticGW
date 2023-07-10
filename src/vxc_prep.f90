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
subroutine prep_vxc
  use gwm
  implicit none
  if(xc_type==1) then
     call read_vxc
  else
     call get_vxc_spn(dens0, n, nsp, vxc)
     vxc0 = vxc
  end if
end subroutine prep_vxc
