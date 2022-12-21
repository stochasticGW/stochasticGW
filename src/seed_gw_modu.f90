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
module seed_gw_modu
  use rand_seed_mod, only : nseed=>nlseed, line_seed
  use rand_seed_mod, only : seed_read, seed_array
  use rand_seed_mod, only : rank_dependent, mdiffseeds, ndiffseeds
  implicit none
  save
  integer, parameter :: line_seed_general = 1
  integer, parameter :: line_seed_pt      = 2
  integer, parameter :: line_seed_gam     = 3
  integer, parameter :: line_seed_g       = 4
  integer, parameter :: line_seed_zeta    = 5
  integer, parameter :: line_seed_exchange= 6
contains
  subroutine check_seed
    implicit none
    integer, save :: i1=1
    if(i1==1) then
       i1=-1
       if(any(rank_dependent(:)/=0)) stop ' rand_dependent not 0 '
       call check_le(6,ndiffseeds,' six, ndiffseeds ')
       call check_le(6,mdiffseeds,' six, mdiffseeds ')
    end if
  end subroutine check_seed
end module seed_gw_modu

