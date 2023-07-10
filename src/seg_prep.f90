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
subroutine seg_prep
  use gwm
  use simple_mpi, only : rank
  implicit none

  if(gamflg/=2) return

  seg = ceiling(seg_frctn*n)
  seg = min(seg, n)
  seg = max(seg, 1)
  seg_fctr = dble(n+seg-1)/seg  ! see gam_seg1

  if(rank==0) then
!     write(6,*)
     write(6,*) ' segment width ',seg
     write(6,*)
  endif
end subroutine seg_prep
