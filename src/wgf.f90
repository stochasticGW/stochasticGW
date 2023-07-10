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
function wgf(it)
  use gwm, only : dt, gama
  implicit none
  integer it
  real*8  wgf  
  wgf = exp(-(gama*dt*it)**2/2.0)
end function wgf
  
