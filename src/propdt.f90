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
subroutine propdt(pt,n,ms,map_sp,ih)
  use simple_mpi, only : rank
  use gwm, only : nx,ny,nz,dv,expnk, nsp
  implicit none
  integer n,ms,it,i,ih
  integer  map_sp(ms)
  complex*16  :: pt(n,ms)
  complex*16, allocatable :: expvks(:,:)
  
  allocate(expvks(n,nsp),stat=i); if(i/=0) stop ' expvks '

  call make_expvks(expvks,n,nsp,ih);
  call propl( pt, expvks,map_sp, n, ms, nsp)
  call propnl(pt, n, ms);
  call propk( pt, expnk, nx, ny, nz, ms);
  call propnl(pt, n, ms)
  call propl( pt, expvks,map_sp, n, ms, nsp)
  call renrmlz(pt,n,ms,dv)
  deallocate(expvks)
end subroutine propdt

subroutine make_expvks(expvks,n,nsp,ih);              
  use gwm, only : vt,vks,dt,rho,flgdyn
  implicit none
  integer  n, nsp, ih, sp
  complex*16  expvks(n,nsp)
  complex*16, parameter :: ci = (0d0,1d0)
  do sp=1,nsp
     select case(ih)
     case(1)
        expvks(:,sp) = exp(-ci*dt/2d0* vt(:,sp)) 
     case(0)
        expvks(:,sp) = exp(-ci*dt/2d0* vks(:,sp))
     case default
        stop ' ih '
     end select
  enddo
end subroutine make_expvks

subroutine propl( pt, expvks,map_sp, n, ms, nsp)
  implicit none
  integer n,ms,map_sp(ms),nsp,is
  complex*16 pt(n,ms),expvks(n,nsp)
  do is=1,ms
     pt(:,is)= pt(:,is)*expvks(:,map_sp(is))
  enddo
end subroutine propl

subroutine propk(pt, expnk, nx,ny,nz,ms)
  implicit none
  integer nx,ny,nz,ms,is
  complex*16 pt(nx*ny*nz,ms), expnk(nx*ny*nz)
  call fft3d_forward_many(nx, ny, nz, ms, pt)
  do is=1,ms
     pt(:,is) = pt(:,is)* expnk(:)
  enddo
  call fft3d_backward_many(nx, ny, nz, ms, pt)
end subroutine propk

subroutine propnl(p,n,ms)
  use simple_mpi, only: rank
  use gwm,        only: na, pt
  implicit none
  integer n, ms, ia, st
  complex*16 p(n,ms)
  complex*16, allocatable :: q(:,:)

  allocate(q(n,ms), stat=st); call check0(st,' q-alloc ' )
  q = 0d0
!  do ia=1,na
!     call prop_nl_a(p,q, n,ms,ia)
!  enddo
  call prop_nl_hamann(p,q,n,ms)
  p=p+q
  deallocate(q)
end subroutine propnl

subroutine renrmlz(pt,n,ms,dv)   
  implicit none
  integer n,ms,is,i
  real*8  dv, nrm
  complex*16 pt(n,ms)
  do is=1,ms
     nrm = 0d0
     do i=1,n
        nrm = nrm + (dble(pt(i,is)))**2+(aimag(pt(i,is)))**2
     enddo
     nrm = 1d0/sqrt(nrm*dv)
     
     do i=1,n
        pt(i,is)=nrm*pt(i,is)
     enddo
  enddo
end subroutine renrmlz

