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
module denek_cut
  implicit none
  save
  real*8 ::  denekcut=8e10
end module denek_cut

subroutine grad(nx, dx, ny, dy, nz, dz, den, nrmgdn)
  use gwm, only: ekcut, n
  use denek_cut, only  : denekcut
  implicit none
  integer             :: nx, ny, nz, st
  real*8              :: dx, dy, dz
  real*8              :: den(nx*ny*nz), nrmgdn(nx*ny*nz)
  real*8, allocatable :: gden(:,:)

  call check(n,nx*ny*nz,' n, nxyz ')
  allocate(gden(n,3),stat=st); if(st/=0) stop ' grad_gden error '

  call grad_vec(nx, dx, ny, dy, nz, dz, den, gden)
  nrmgdn = sqrt(sum(gden**2,dim=2))

  deallocate(gden)
end subroutine grad

!
! check the signs, maybe need to put an overall minus.
!

subroutine grad_vec(nx, dx, ny, dy, nz, dz, den, gden)
  use gwm, only: ekcut, n
  use denek_cut, only  : denekcut
  implicit none
  integer     :: nx, ny, nz, ix, iy, iz, st, i
  real*8      :: dx, dy, dz, k2
  real*8      :: dkx, dky, dkz, kx, ky, kz
  real*8      :: den(nx*ny*nz)
  real*8      :: gden(nx*ny*nz,3)
  real*8      :: pi, tpi
  complex*16, allocatable, dimension(:) :: vk,vo1,vo2,vo3

  call check(n,nx*ny*nz,' n, nxyz ')
  allocate(vk(n), stat=st); if(st/=0) stop ' gradvk '
  allocate(vo1(n),stat=st); if(st/=0) stop ' gradv1 '
  allocate(vo2(n),stat=st); if(st/=0) stop ' gradv2 '
  allocate(vo3(n),stat=st); if(st/=0) stop ' gradv3 '

  pi = dacos(-1d0)
  tpi = 2d0*pi
  !later remove this whole subroutine and use only grad_cf

  dkx = tpi/(dble(nx)*dx)
  dky = tpi/(dble(ny)*dy)
  dkz = tpi/(dble(nz)*dz)

  vk = den
  call fft3d_forward(nx, ny, nz, vk, vk)

  i=0
  do iz = 1,nz
     do iy = 1,ny
        do ix = 1,nx
           i=i+1
           kx = dble(ix-1) * dkx
           ky = dble(iy-1) * dky
           kz = dble(iz-1) * dkz
           
           if(kx>pi/dx) kx = kx-tpi/dx
           if(ky>pi/dy) ky = ky-tpi/dy
           if(kz>pi/dz) kz = kz-tpi/dz
           
           k2 = kx**2.d0 + ky**2.d0 + kz**2.0
           if(k2.gt.denekcut*ekcut) then
              vk(i) = 0d0
           end if

           if(ix==nx/2+1) kx=0d0  !intentionally after the def. of k2.
           if(iy==ny/2+1) ky=0d0
           if(iz==nz/2+1) kz=0d0

           vo1(i) = vk(i)* kx
           vo2(i) = vk(i)* ky
           vo3(i) = vk(i)* kz
           
        end do
     end do
  end do
  
  call fft3d_backward(nx, ny, nz, vo1, vo1) 
  call fft3d_backward(nx, ny, nz, vo2, vo2) 
  call fft3d_backward(nx, ny, nz, vo3, vo3) 
  
  gden(:,1) = -aimag(vo1)/dble(n)
  gden(:,2) = -aimag(vo2)/dble(n)
  gden(:,3) = -aimag(vo3)/dble(n)
  
  deallocate(vo1,vo2,vo3)
end subroutine grad_vec

subroutine divrg(nx, dx, ny, dy, nz, dz, vec, div)
  use gwm, only: ekcut, n
  use denek_cut, only : denekcut
  implicit none
  integer     :: nx, ny, nz, ix, iy, iz, i, st
  real*8      :: dx, dy, dz, k2
  real*8      :: dkx, dky, dkz, kx, ky, kz
  real*8      :: vec(nx*ny*nz,3), div(nx*ny*nz)
  real*8      :: pi, tpi
  complex*16, allocatable :: vo1(:), vo2(:), vo3(:), vo(:)

  allocate(vo1(n), stat=st); if(st/=0) stop ' div_vo1 error '
  allocate(vo2(n), stat=st); if(st/=0) stop ' div_vo2 error '
  allocate(vo3(n), stat=st); if(st/=0) stop ' div_vo3 error '
  allocate( vo(n), stat=st); if(st/=0) stop ' div_vo  error '

  pi = dacos(-1d0)
  tpi = 2d0*pi

  dkx = tpi/(dble(nx)*dx)
  dky = tpi/(dble(ny)*dy)
  dkz = tpi/(dble(nz)*dz)

  vo1 = vec(:,1)
  call fft3d_forward(nx, ny, nz, vo1, vo1)

  vo2 = vec(:,2)
  call fft3d_forward(nx, ny, nz, vo2, vo2)

  vo3 = vec(:,3)
  call fft3d_forward(nx, ny, nz, vo3, vo3)

  i=0
  do iz = 1,nz
     do iy = 1,ny
        do ix = 1,nx
           i=i+1

           kx = dble(ix-1) * dkx
           ky = dble(iy-1) * dky
           kz = dble(iz-1) * dkz
           
           if(kx>pi/dx) kx = kx-tpi/dx
           if(ky>pi/dy) ky = ky-tpi/dy
           if(kz>pi/dz) kz = kz-tpi/dz
           
           k2 = kx**2.d0 + ky**2.d0 + kz**2.0
           if(k2.gt.denekcut*ekcut) then
              vo1(i) = 0d0
              vo2(i) = 0d0
              vo3(i) = 0d0
           end if
           
           if(ix==nx/2+1) kx=0d0  !intentionally after the def. of k2.
           if(iy==ny/2+1) ky=0d0
           if(iz==nz/2+1) kz=0d0

           vo1(i) = vo1(i) * kx
           vo2(i) = vo2(i) * ky
           vo3(i) = vo3(i) * kz
        end do
     end do
  end do
  
  vo = vo1 + vo2 + vo3
  call fft3d_backward(nx, ny, nz, vo, vo) 
  div(:) = -aimag(vo(:)) / dble(n)   ! check signs, they may be wrong. dont erase

  deallocate(vo1,vo2,vo3,vo)
end subroutine
  
