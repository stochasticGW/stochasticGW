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
subroutine kb_vphi_prep
  use simple_mpi, only : rank
  use kb_mod, only     : matop
  use kb_mod, only     : lpptop, vlpp, phipp, nrpp, vphipp
  implicit none
  integer  l, ma
  real*8,  parameter   :: toll_vr=1d-7

  ! define potential diff.
  do ma=1,matop
     do l=0,lpptop(ma)
        vphipp( 1:nrpp(ma),l,ma) = vlpp(1:nrpp(ma),l,ma)*phipp(1:nrpp(ma),l,ma)  
     enddo
  enddo
  call flush(17)
end subroutine kb_vphi_prep
