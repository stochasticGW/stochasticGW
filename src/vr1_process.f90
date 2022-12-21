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
subroutine vr1_process(i,it)
  use gwm
  implicit none
  integer  i, it
  select case(gamflg)
  case(1,2); 
     call gamvrmake
     call fr1_gama(i,it)
  case(3); call vr1_prnt(i,it)
  case default; stop ' gamflg '
  end select
end subroutine vr1_process
