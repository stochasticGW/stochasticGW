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
subroutine vk_possibly_reset
  use gwm, only : fuzzy_vk, vk, vk_exch
  implicit none
  call check(size(vk),size(vk_exch),' sz_vk_vk_exh ')
  if(fuzzy_vk) vk=vk_exch
end subroutine vk_possibly_reset
