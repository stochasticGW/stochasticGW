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
subroutine prep_rand
  use seed_gw_modu, only : line_seed, line_seed_general
  implicit none
  real*8 r
  real*8, external :: ran_ps
  line_seed = line_seed_general 
  r = ran_ps()
end subroutine prep_rand
