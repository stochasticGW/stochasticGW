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
subroutine pp_bcast
  use simple_mpi, only : bcast_scalar_i
  use simple_mpi, only : bcast_i
  use simple_mpi, only : bcast_r8
  use simple_mpi, only : bcast_L
  use ppm, only: nrppmx
  use ppm, only: nrpp
  use ppm, only: lpptop
  use ppm, only: lpploc
  use ppm, only: rrpp
  use ppm, only: phipp
  use ppm, only: vphipp
  use ppm, only: vlpp
  use ppm, only: vpploc
  use ppm, only: rho_core_a
  use ppm, only: core_correction
  implicit none
  call bcast_scalar_i(nrppmx)
  call bcast_i(  nrpp,    size(nrpp),  0)
  call bcast_i(  lpptop,  size(lpptop),0)
  call bcast_i(  lpploc,  size(lpploc),0)
  call bcast_r8( rrpp,    size(rrpp),  0)
  call bcast_r8( vlpp,    size(vlpp),  0)
  call bcast_r8( vpploc,  size(vpploc),0)
  call bcast_r8( phipp,   size(phipp), 0)
  call bcast_r8(vphipp,   size(vphipp),0)
  call bcast_r8(rho_core_a, size(rho_core_a),0)
  call bcast_L(core_correction, size(core_correction),0)
end subroutine pp_bcast
