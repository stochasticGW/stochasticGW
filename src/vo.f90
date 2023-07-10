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
subroutine vo_make     
  implicit none
  call vr_make 
  call vr_to_vo
end subroutine vo_make
  
subroutine vr_to_vo
  use gwm
  implicit none
  select case(gamflg)
  case(1,2);    call vr_to_vo_gam
  case(3);      call vr_to_vo_drct
  case default; stop ' gamflg '
  end select
end subroutine vr_to_vo
