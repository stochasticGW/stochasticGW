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
subroutine gam_seg
  use gwm
  implicit none
  integer j, st

  if(.not.allocated(seg_sh))then; allocate(seg_sh(ngam),stat=st); call check0(st,' seg_sh ')  
  end if

  if(.not.allocated(seg_w)) then; allocate(seg_w(ngam),stat=st); call check0(st,' seg_w ')  
  end if

  if (block_gam_alg) then
     call gam_seg1_block
  else
     do j=1,ngam
        call gam_seg1(seg_sh(j), seg_w(j))
     end do
  endif 

end subroutine gam_seg
