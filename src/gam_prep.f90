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
subroutine gam_prep
  use gwm
  implicit none
  integer igam, st

  if(ngam<1) return
  if(gamflg==2) call gam_seg
  if(.not.allocated(gam)) then
     allocate(gam(seg,ngam),stat=st); call check0(st,' gam ')
  endif
  do igam=1,ngam
     call set_seed_gam( igam)
     call rand_r48(gam(:,igam),size(gam,1), dv)
  end do
end subroutine gam_prep
