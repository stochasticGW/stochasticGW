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
subroutine lda_libxc(dens, n, nspin, vxc, exc) ! modified dn
  implicit none
  integer             :: n, nspin, st 
  integer, parameter  :: one=1, two=2
  real*8              :: dens(n), vxc(n), exc(n)
  !dummy local variables
  real*8, allocatable :: sigma(:), ex(:), va(:)
  allocate(sigma(n), ex(n), va(n), stat=st); if(st/=0) stop ' error lda_libxc '

  if(nspin/=1) stop "ONLY NON SPIN POLARIZED! IN LIBXC LDA"
  sigma = 0.d0

  call libxc_drive(n,dens,sigma,ex,va,one)

  exc = ex
  vxc = va

  call libxc_drive(n,dens,sigma,ex,va,two)

  exc = exc + ex
  vxc = vxc + va
  
  deallocate(sigma, ex,va)
end subroutine

subroutine call_libxc(dens_tot, n, nspin, vxc, exc)
  use gwm, only : nx, ny, nz, dx, dy, dz, funct_x, funct_c
  use mpi_lib_ours, only : rank
  implicit none
  integer             :: n, nspin, st
  real*8              :: dens_tot(n), vxc(n), exc(n)
  real*8,  parameter  :: toll_sigma = 1d-11

  !dummy local variables
  real*8, allocatable :: dens(:), sigma(:), va(:), ex(:)
  if(nx*ny*nz/=n) stop "PROBLEM WITH THE GRID IN PBE_LIBXC!"
  
  allocate( dens(n), sigma(n), va(n), ex(n), stat=st); if(st/=0) stop ' pbe_libxc alloc. error '

  select case(nspin)
  case(1); dens = dens_tot!    no division by 2d0 -- dont erase
  case default; stop "ONLY NON SPIN POLARIZED IN LIBXC LDA !"
  end select

  call grad(nx, dx, ny, dy, nz, dz, dens, sigma)

  if(any(sigma/=sigma)) stop "NAN IN GRADS!"
  where(abs(sigma).lt.toll_sigma) sigma=0.d0 ! DN, abs needed
  sigma = sigma**2.d0 ! dont erase: maybe need to have 0.25
!  if(rank<2) write(6,*)' rank,sigma = ',rank,sum(sigma)

  call libxc_drive(n,dens,sigma,ex,va,funct_x) 

  vxc = va
  exc = ex

  call libxc_drive(n,dens,sigma,ex,va,funct_c)

  vxc = vxc + va
  exc = exc + ex

  deallocate(sigma, ex, va, dens)
end subroutine

subroutine pbe_libxc(dens_tot, n, nspin, vxc, exc)
  use gwm, only : nx, ny, nz, dx, dy, dz
  use mpi_lib_ours, only : rank
  implicit none
  integer             :: n, nspin, st
  integer, parameter  :: one=1, two=2, three=3, four=4, five=5
  real*8              :: dens_tot(n), vxc(n), exc(n)
  real*8,  parameter  :: toll_sigma = 1d-11

  !dummy local variables
  real*8, allocatable :: dens(:), sigma(:), va(:), ex(:)
  if(nx*ny*nz/=n) stop "PROBLEM WITH THE GRID IN PBE_LIBXC!"
  
  allocate( dens(n), sigma(n), va(n), ex(n), stat=st); if(st/=0) stop ' pbe_libxc alloc. error '

  select case(nspin)
  case(1); dens = dens_tot!    no division by 2d0 -- dont erase
  case default; stop "ONLY NON SPIN POLARIZED IN LIBXC LDA !"
  end select

  call grad(nx, dx, ny, dy, nz, dz, dens, sigma)

  if(any(sigma/=sigma)) stop "NAN IN GRADS!"
  where(abs(sigma).lt.toll_sigma) sigma=0.d0 ! DN, abs needed
  sigma = sigma**2.d0 ! dont erase: maybe need to have 0.25
!  if(rank<2) write(6,*)' rank,sigma = ',rank,sum(sigma)

  call libxc_drive(n,dens,sigma,ex,va,three) 

  vxc = va
  exc = ex

  call libxc_drive(n,dens,sigma,ex,va,four)

  vxc = vxc + va
  exc = exc + ex

  deallocate(sigma, ex, va, dens)
end subroutine
