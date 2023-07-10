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
!
! rconv:  subroutine to calculate fout = 1/N * ifft3d( bk * fft3d( fin))
!
! where "fin" and "fout" are real, and bk is real (and presumed to be bk(kvec) = bk(-kvec) 
! 
! The routine would also have worked with a tiny input change if bk was complex, presuming
! bk(-kvec)=bk*(kvec)
!
! Note that fout, fin, and bk are presumed to be (nx,ny,nz); the division to
! 2*(nx/2+1),ny,nz is done internally.
! 
!
subroutine rconv(gin,bk,gout,nx,ny,nz)
  implicit none
  integer nx, ny, nz
  integer ix, iy, iz, st
  real*8 gin( nx, ny, nz)
  real*8   bk(nx, ny, nz)
  real*8 gout(nx, ny, nz)
  real*8, allocatable :: ar(:,:,:)
  real*8 :: cnst
  cnst = 1d0/( dble(nx)*dble(ny)*dble(nz) )
  allocate(ar(2*(nx/2+1),ny,nz), stat=st); call check0(st,' ar_alloc')
  ar(1:nx,           :,:)= gin
  ar(nx+1:2*(nx/2+1),:,:)=0d0
  call fft3d_r2cinplace(nx, ny, nz, ar)
  call multip_ab(ar,bk)
  call fft3d_c2rinplace(nx, ny, nz, ar)
  gout = ar(1:nx,:,:)*cnst
  deallocate(ar)
contains
  subroutine multip_ab(ar, bk)
    implicit none
    real*8 ar( 2,  nx/2+1,ny,nz)
    real*8     bk(nx,    ny,nz)
    ar(1,:,:,:) = ar(1,:,:,:) * bk(1:nx/2+1,:,:)
    ar(2,:,:,:) = ar(2,:,:,:) * bk(1:nx/2+1,:,:)
  end subroutine multip_ab
end subroutine rconv

