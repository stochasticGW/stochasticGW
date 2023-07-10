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
  allocate(iub( na),              stat=st); call check0(st,' iub ')
  allocate(iut( na),              stat=st); call check0(st,' iut ')
  allocate(racut_a(       matop), stat=st); call check0(st,' racut_a ')
  allocate(pvp(    0:lmx, matop), stat=st); call check0(st,' vkb_nrmlz ')
end subroutine kb_alloc

subroutine kb_alloc_vp
  use kb_mod,     only : iumx, vp, st
  implicit none
  allocate(vp(iumx),  stat=st);  call check0(st,' vp       ')
end subroutine kb_alloc_vp

  
  
  
