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
subroutine pp_allocate
  use ppm,   only : nrpp, nrppmx
  use ppm,   only : lpploc, lpptop
  use ppm,   only : rrpp, vpploc, vphipp, phipp, vlpp
  use ppm,   only : rho_core_a, core_correction
  use ppm,   only : matop
  use ppm,   only : lmx
  implicit none
  integer st
  call check_lele(50,nrppmx,10000, '50-nrppmx-1e4 ')
  
  allocate( lpploc(                   matop), stat=st);      call check0(st,' lpploc ')
  allocate( lpptop(                   matop), stat=st);      call check0(st,' lpploc ')
  allocate( nrpp(                     matop), stat=st);      call check0(st,' nrpp   ')
  allocate( rrpp(      1:nrppmx,      matop), stat=st);      call check0(st,' rrpp   ')
  allocate( rho_core_a(1:nrppmx,      matop), stat=st);      call check0(st,' rho_core_a ')
  allocate( vpploc(    1:nrppmx,      matop), stat=st);      call check0(st,' vpploc ')
  allocate( vlpp(      1:nrppmx,0:lmx,matop), stat=st);      call check0(st,' vlpp   ')
  allocate( phipp(     1:nrppmx,0:lmx,matop), stat=st);      call check0(st,' vpploc ')
  allocate( vphipp(    1:nrppmx,0:lmx,matop), stat=st);      call check0(st,' vpploc ')
  allocate( core_correction(          matop), stat=st);      call check0(st,' core_correction ')

  nrpp = 0 ! to be checked later that nrpp>0, i.e., all entries were read.
  core_correction = .false.
end subroutine pp_allocate
