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
subroutine kb_count
  use kb_mod, only : na                           ! inputs
  use kb_mod, only : iub, iut, iumx               ! outputs
  use kb_mod,  only : ch_a
  implicit none
  !
  ! counting the number of point for the non-local KB potential
  !
  integer iu,ia
  iu=0
  do ia=1,na
     iub(ia)=iu+1
     call kb_count_1a(ia,iu)
     iut(ia)=iu
  end do
  iumx = iu
end subroutine kb_count
