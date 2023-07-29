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
subroutine vo1(it)
  use gwm
  implicit none
  integer it
  select case(gamflg)
  case(1,2); call vo1_gam(it)
  case(3);   call vo1_drct(it)
  case default; stop ' gamflg-vo1 '
  end select
end subroutine vo1
