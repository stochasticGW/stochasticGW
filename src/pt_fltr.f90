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
subroutine pt_fltr
  use simple_mpi, only : rank
  use gwm
  implicit none
  integer   :: ib, isp, slc
  slc = 1
  if(det_tddft) stop ' ERROR: need stoch tddft in stoc_pt '

  if(filter_cheby) then
   do ib=1,ns
    isp = map_sp_pt(ib)
    call filter_chebyshev(isp, vks(:,isp),&
                          n,slc,pt(1,ib),normpt(ib),dh,havg,th_co,nchbc,nchbmx,dv,rank)
   enddo
  else
     stop ' problem in pt_fltr '
  endif
   
  call norm_check_lowest

contains
  subroutine norm_check_lowest
    implicit none
    !
    ! norm check for lowest orbital
    !
    
    if(ns>0) then
       if(normpt(1)>1d-10) then
          continue
       else
          write(6,*)' ERROR: normpt =',normpt,' too tiny.  mu=',mu,' is too small '
          stop
       endif
       
       if(rank==0) then
          write(17,  *)' filtered: rank=0, normpt ',ib,normpt
          call flush(17)
       endif
    end if
    call prnt(1,' pt_fltr ')
  end subroutine norm_check_lowest
end subroutine pt_fltr

