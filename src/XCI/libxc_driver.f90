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
subroutine libxc_drive(ngrp,dens,sigma,exc,vxc,func)
  use xc_f03_lib_m
  use gwm, only : nx, ny, nz, dx, dy, dz
  use mpi_lib_ours, only : rank, sync_mpi
  implicit none

  !local dummy
  integer   :: ix, iy, iz, ii, st

  !from the DFT code
  integer :: ngrp
  integer*8 :: ngrp8
  real*8  :: dens(ngrp), sigma(ngrp), exc(ngrp), vxc(ngrp)
  ! modified -- allocatable

  real*8, allocatable :: v1mul(:), v2mul(:), vsigma(:), den_eps(:)
  integer :: func
  ! %% functional list allowed right now
  !    1 - LDAx
  !    2 - LDAc (Perdew & Wang) 
  !    3 - PBEx revised 
  !    4 - PBEc 
  !    5 - wPBEx



  !local vars for linking
!  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f03_func_t) :: xc_func
  TYPE(xc_f03_func_info_t) :: xc_info
!  TYPE(xc_f90_pointer_t) :: xc_info
  integer :: i, vmajor, vminor, vmicro, func_id 
  real*8 :: zero_r8=0.d0, one_r8=1.d0

  func_id=func
!  select case(func)
!  case(1)
!    func_id=1
!  case(2) 
!    func_id=12
!  case(3)
!    func_id=101
!  case(4)
!    func_id=130
!  case(5)
!    func_id=524 !478 !524=hse LR !478=LC_WPBE by Vydrov&Scuseria
!  case default
!    stop "WRONG SELECTION OF FUNCTIONAL IN LIBXC_DRIVER!"
!  end select

  allocate(v1mul(ngrp), v2mul(ngrp), vsigma(ngrp),den_eps(ngrp), stat=st)
  if(st/=0) stop ' libxc_driver allocation error '
  
  ngrp8=INT(ngrp,KIND(ngrp8))

!  call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)
  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
  xc_info = xc_f03_func_get_info(xc_func)

!  select case (xc_f90_info_family(xc_info))
  select case (xc_f03_func_info_get_family(xc_info))
  case(XC_FAMILY_LDA)
     call xc_f03_lda_exc_vxc(xc_func, ngrp8, dens(1), exc(1), vxc(1))
  case(XC_FAMILY_GGA)
     !
     ! note the structure:
     !  vmul is the multiplicative part of the potential
     !  the total potential has another term: -2 Div* (vsigma*Grad(dens))
     !
     
     call xc_f03_gga_exc_vxc(xc_func, ngrp8, dens(1), sigma(1), exc(1), v1mul(1), vsigma(1))
     
     
     where(sigma==0) vsigma=0.d0
     den_eps = dens
     where(den_eps.lt.1E-10) den_eps=0.d0
     
     call gga_pot2(ngrp, den_eps, vsigma,v2mul)
     
     vxc = v1mul + v2mul
     
  end select
  
  call xc_f03_func_end(xc_func)
  
  deallocate(v1mul, v2mul, vsigma, den_eps)

end subroutine

!
! removed denshlf
!

! may need need to use  1/2 for dens if 0.25 in sigma
subroutine gga_pot2(n, dens, vsigma,v2mul)
  use gwm,      only : nx, ny, nz, dx, dy, dz
  implicit none
  integer   :: idir, n, st
  real*8    :: dens(n), vsigma(n), v2mul(n)
  real*8, allocatable :: vec(:,:)

  allocate(vec(n,3), stat=st); if(st/=0) stop 'error gga_pot2  '
  call grad_vec(nx, dx, ny, dy, nz, dz, dens, vec)

  do idir = 1,3
    vec(:,idir) = vsigma(:) * vec(:,idir)  
  end do

  call divrg(nx, dx, ny, dy, nz, dz, vec, v2mul)

  v2mul = -2.d0*v2mul !factor of two due to derivative of |nabla n|^2
  deallocate(vec)

end subroutine gga_pot2

