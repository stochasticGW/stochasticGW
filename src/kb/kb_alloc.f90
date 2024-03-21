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
subroutine kb_alloc
  use kb_mod
  implicit none
  allocate(ngs( na),              stat=st); call check0(st,' ngs ')
  allocate(racut_a(       matop), stat=st); call check0(st,' racut_a ')
end subroutine kb_alloc

