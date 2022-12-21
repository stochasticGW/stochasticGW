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
subroutine prep_vk_3d( nx,ny,nz,dx,dy,dz,vk)
  implicit none
  integer nx,ny,nz
  integer ix,iy,iz
  integer st
  real*8 vk(nx,ny,nz)
  real*8 dx,dy,dz
  real*8 kx,ky,kz
  real*8 dkx,dky,dkz
  real*8 k2
  real*8 pi
  real*8, save :: cut=1d-8
  pi = dacos(-1.d0)

  dkx = 0d0 ; if(dx>1e-8) dkx = 2.d0*pi/(dble(nx)*dx)
  dky = 0d0 ; if(dy>1e-8) dky = 2.d0*pi/(dble(ny)*dy)
  dkz = 0d0 ; if(dz>1e-8) dkz = 2.d0*pi/(dble(nz)*dz)

  vk = 0d0
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

           if(k2>cut) then
              vk(ix,iy,iz) =  4d0*pi/k2
           end if

       enddo kxdo
    enddo kydo
 enddo kzdo

end subroutine prep_vk_3d

