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
subroutine eta_fltr
  use simple_mpi, only : rank, color_rank
  use gwm
  implicit none
  integer                 :: st
  integer, save           :: n1=1
  complex*16, allocatable :: cz(:)
  real*8,     allocatable :: normz(:),rz(:)

  if(det_dft) stop           ' filter_eta only for non-(det_dft) '
  if(.not.filter_cheby) stop ' filter_eta needs filter_cheby     '
  if(color_rank==eta_rank) then
     allocate(cz(n), normz(1), stat=st); call check0(st,' cz ')
     cz = zeta
     call filter_chebyshev(sp0, vks(:,sp0),n,n1,cz,normz,dh,havg,th_co, nchbc, nchbmx, dv,rank)
     if(normz(1)<1d-10) then
        write(6,*)' ERROR: normz =',normz,' too tiny.  mu=',mu,' is too small '
        stop
     endif
     eta = cz * normz(1) !  
     if(abs(sum(abs(eta)**2*dv)-normz(1)**2)>1.d-6*normz(1)**2)stop ' tmp eta problems '
     deallocate(cz,normz)
  end if
end subroutine eta_fltr
