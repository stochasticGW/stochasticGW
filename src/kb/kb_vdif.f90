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
subroutine kb_vdif
  use simple_mpi, only : rank
  use kb_mod, only     : matop
  use kb_mod, only     : lpptop, lpploc, vl=>vlpp, pl => phipp, nrpp, vpploc
  implicit none
  integer  l, ma
  real*8,  parameter   :: toll_vr=1d-7

  ! define potential diff.
  do ma=1,matop
     do l=0,lpptop(ma)
        vl( 1:nrpp(ma),l,ma) = vl(1:nrpp(ma),l,ma)-vpploc(1:nrpp(ma),ma) !dont exclude l=lloc 
     enddo
  enddo

  if(rank==0) then
     do ma=1,matop
        do l=0,lpptop(ma)
           write(17,*)' ma,l,vl_min vl_max ',ma,l,minval(vl(1:nrpp(ma),l,ma)),&
                                                  maxval(vl(1:nrpp(ma),l,ma))
        enddo
     enddo
     call flush(17)
  endif
end subroutine kb_vdif
