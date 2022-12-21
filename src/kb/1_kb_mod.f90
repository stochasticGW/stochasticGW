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
module kb_mod
  use simple_mpi,    only : rnk => rank
  use atoms,         only : mapat, mapai, cnt, matop, ch_a
  use kb_top_module, only : ng=>n, nx, ny, nz, dx, dy, dz, dv
  use kb_top_module, only : na, lmx, scale_vh, dim_periodic
  use kb_top_module, only : vloc_tot
  use ppm,           only : nrpp, rrpp, lpptop, lpploc, vlpp, vpploc, phipp, vphipp
  implicit none
  save
  integer                  :: st
  integer                  :: ngsmall
  integer                  :: iumx
  integer, allocatable     :: ngs(:)
  integer, allocatable     :: iub(:)
  integer, allocatable     :: iut(:)
  integer, allocatable     :: mapkbg(:,:)
  real*8,  allocatable     :: rgn(:,:,:)
  real*8, allocatable      :: racut_a(:)
  real*8, allocatable      :: vp(:)
  real*8, allocatable      :: pvp(        :,:)
  complex*16, allocatable  :: cvn(:,:,:)       ! for formafactor
  complex*16, allocatable  :: cvkb(:,:,:)      ! for formafactor.  Really real
end module kb_mod

