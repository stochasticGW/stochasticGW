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
  use ppm,           only : nrpp, rrpp, lpptop, lpploc, vlpp, vpploc, phipp !, vphipp
  use ppm,           only : nproj_m, lpp_m, phipp_m, nsuper_ianl, indx_ianl, start_ianl
  use ppm,           only : vp_hamann, dij_diag_m
  implicit none
  save
  integer                  :: st
  integer                  :: ngsmall
  integer                  :: iumx
  integer, allocatable     :: ngs(:)
  integer, allocatable     :: mapkbg(:,:)
  real*8,  allocatable     :: rgn(:,:,:)
  real*8, allocatable      :: racut_a(:)
  complex*16, allocatable  :: cvn(:,:,:)       ! for formafactor
  complex*16, allocatable  :: cvkb(:,:,:)      ! for formafactor.  Really real
end module kb_mod

