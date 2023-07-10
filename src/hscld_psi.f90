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
subroutine hc_scld_comp(isp, p1,p0,tmp,vks,n,ms,havg,dh)
! in p1,p0,tmp ; out p1,p0<=2*W*p1-p0,tmp_modified
  use gwm, only : nx,ny,nz,ekn,dv !,functional
  use simple_mpi, only : rank
  implicit none
  integer     n, ms, i, is, isp
  real*8      havg, dh, vks(n), two_ov_dh
  complex*16  p1(n,ms),p0(n,ms),tmp(n,ms)
  two_ov_dh = 2d0/dh

  if(ms/=1) then
    write(6,*) "ERROR: Only single spin channel is supported right now!"
    call flush(6)
    stop
  end if
  do is=1, ms
     p0(:,is) = two_ov_dh*(vks(:)-havg)*p1(:,is)-p0(:,is)   ! 2(v-havg)/dh
  enddo
  call addvnl_c(p1,p0,n,ms,two_ov_dh)   ! p0=p0+ vnlscaled p1
  tmp(:,1:ms) = p1(:,1:ms)
  call fft3d_forward_many(nx, ny, nz, ms, tmp)
  do is=1,ms
     tmp(:,is) = tmp(:,is)* ekn(:)
  enddo
  call fft3d_backward_many(nx, ny, nz, ms, tmp)
  p0(:,1:ms) = p0(:,1:ms)+two_ov_dh*tmp(:,1:ms)

!  if(functional=='bnl') then
!     tmp = 0d0
!     call add_xc_long_range_bnl(isp, p1, tmp, n, ms)
!     p0(:,1:ms) = p0(:,1:ms)+two_ov_dh*tmp(:,1:ms)
!  end if

end subroutine hc_scld_comp


