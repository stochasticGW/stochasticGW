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
module ppm
  use gwm,   only : lmx
  use atoms, only : matop, mapat
  use atoms, only : atom_name, atom_upper_name
  implicit none
  save
  integer              :: nrppmx
  integer, allocatable :: nrpp(:)
  integer, allocatable :: lpptop(:)   
  integer, allocatable :: lpploc(:)
  real*8,  allocatable :: rrpp(:,:)
  real*8,  allocatable :: vpploc(:,:)
  real*8,  allocatable :: phipp(:,:,:)
  real*8,  allocatable :: vlpp(:,:,:)
  real*8,  allocatable :: vphipp(:,:,:)
  real*8,  allocatable :: rho_core_a(:,:)
  logical, allocatable :: core_correction(:)

end module ppm
