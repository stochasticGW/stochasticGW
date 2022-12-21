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
! Tuckerman Martyna
subroutine prep_vk_MT( nx,ny,nz,dx,dy,dz,vk)
  use simple_mpi, only : rank
  implicit none
  integer nx,ny,nz
  integer ix,iy,iz
  integer st
  real*8 vk(nx,ny,nz)
  real*8 dx,dy,dz
  real*8 kx,ky,kz
  real*8 dkx,dky,dkz
  real*8 r0(3),rdif(3)
  real*8 dv,rp,k2
  real*8, save :: aa
  real*8 pi
  real*8 k2s
  real*8, save :: cut=1d-8
  complex*16 cterm
  complex*16, allocatable, dimension(:,:,:)  :: cvk
  complex*16, parameter :: ci = (0.d0, 1.d0)
  pi = dacos(-1.d0)

  allocate(cvk(nx,ny,nz),stat=st); if(st/=0)stop 'cvk '  
  aa = 7d0/min(nx*dx,ny*dy,nz*dz)
  r0 = (/ -(nx)/2d0*dx, -(ny)/2d0*dy, -(nz)/2d0*dz /)

  if(rank==0)then
     write(17,*)' aa ( tuma ),r0 = ',aa,r0
     call flush(17)
  endif

  dv = 1d0
  if(dx>1e-8)dv=dv*dx
  if(dy>1e-8)dv=dv*dy
  if(dz>1e-8)dv=dv*dz

  dkx = 0d0 ; if(dx>1e-8) dkx = 2.d0*pi/(dble(Nx)*dx)
  dky = 0d0 ; if(dy>1e-8) dky = 2.d0*pi/(dble(Ny)*dy)
  dkz = 0d0 ; if(dz>1e-8) dkz = 2.d0*pi/(dble(Nz)*dz)
  
  !
  ! separate : 1/r = erf(aa*r)/r + (1-erf(ar))/r =  long-range-in-r + short-range in r.
  !
  !  FFT numerically the first part
  !

  cvk = 0d0
     
  do iz=1,nz
     do iy=1,ny
        do ix=1,nx
           rdif = r0 + (/(ix-1)*dx,(iy-1)*dy,(iz-1)*dz/)
           rp  = sqrt(max(sum(rdif**2),1e-26))
           cvk(ix,iy,iz) = erf(aa*rp)/rp * dv
        enddo
     enddo
  enddo
     
  call fft3d_general(cvk, Nx, Ny, Nz, 1)       
     
  kzdo : do iz=1,nz
     kydo : do iy=1,ny
        kxdo : do ix=1,nx
           
           kx = dble(ix-1)* dkx
           ky = dble(iy-1)* dky
           kz = dble(iz-1)* dkz
           
           if(kx>pi/dx) kx = kx-2d0*pi/dx  
           if(ky>pi/dy) ky = ky-2d0*pi/dy
           if(kz>pi/dz) kz = kz-2d0*pi/dz
           
           k2 = kx**2+ky**2+kz**2        
              
           cvk(ix,iy,iz) = &
           cvk(ix,iy,iz) * exp(ci*( r0(1)*kx+r0(2)*ky+r0(3)*kz ))   
              
           ! now to short-range in r.
            
           k2s = k2/(4d0*aa**2d0)
              
           if(k2s>cut) then
              cterm = 4.d0*pi/k2*(1.d0-exp(-k2s))
           else
              ! 1-exp(-k2s) = 1-(1-k2s+ks2^2/2 - k2s^3/6) = k2s-0.5d0*k2s**2 +k2s**3/6d0
              cterm = 4d0*pi/(4d0*aa**2d0)*(1d0-0.5d0*k2s +k2s**2/6d0)
           end if

           cvk(ix,iy,iz) =  cvk(ix,iy,iz) +  cterm   
           
       enddo kxdo
    enddo kydo
 enddo kzdo

 call check_c_real(cvk, size(cvk), ' cv_num ')
 vk =  cvk
 deallocate(cvk,stat=st); call check0(st,' cvkdeal ')
 
end subroutine prep_vk_MT

