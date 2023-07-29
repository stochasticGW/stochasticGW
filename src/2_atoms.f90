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
module atoms
  use gwm, only : mapai
  implicit none
  save
  integer, parameter    :: n118 = 118
  integer               :: matop
  integer, allocatable  :: mapat(:)
  integer, allocatable  :: ch_a(:)
  integer, allocatable  :: valch_a(:)
  real*8,  allocatable  :: cnt(:,:)
  ! stuff below used in root only
  character*2           ::       atom_name(n118)
  character*2           :: atom_upper_name(n118)
  integer               :: valence_list(n118)
  logical               :: found_charge(n118)=.false.
end module atoms

