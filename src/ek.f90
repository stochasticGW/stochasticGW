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
subroutine ek_prep
  use gwm
  implicit none
  integer i
  integer ix, iy, iz
  real*8  kx, ky, kz
  real*8  dkx, dky, dkz
  real*8  pi
  complex*16, parameter :: ci = (0d0,1d0)
  pi = dacos(-1d0)
  
  dkx = 2d0*pi/(nx*dx)
  dky = 2d0*pi/(ny*dy)
  dkz = 2d0*pi/(nz*dz)

  i= 0
  do iz=1,nz
     do iy=1,ny
        do ix=1,nx

           i=i+1

           kx = (ix-1)*dkx
           ky = (iy-1)*dky
           kz = (iz-1)*dkz
           
           if(kx>pi/dx) kx = kx - 2d0*pi/dx
           if(ky>pi/dy) ky = ky - 2d0*pi/dy
           if(kz>pi/dz) kz = kz - 2d0*pi/dz

           ekn(i) = min(ekcut, (kx**2+ky**2+kz**2)/2d0)

           expnk(i) = exp(-ci*dt*ekn(i))/(dble(nx)*dble(ny)*dble(nz))

           ekn(i)=ekn(i)/(dble(nx)*dble(ny)*dble(nz))
        enddo
     enddo
  enddo
end subroutine ek_prep
