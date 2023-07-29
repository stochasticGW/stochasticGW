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
subroutine check_prnt
  use gwm, only : print_now, nmc, imc
  implicit none
  integer, external :: is_i_near_power_x
  real*8            :: r
  !
  ! prints for imc (monte carlo iteration on core) below 32 or multiples of 8 or final ones 
  !
  print_now = .false.
  if(imc.le.32.or.mod(imc,8)==0.or.imc==nmc) print_now = .true.
end subroutine check_prnt
