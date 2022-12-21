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
subroutine bcast_atoms
  use atoms
  use simple_mpi, only : bcast_scalar_i, bcast_i
  implicit none
  integer st
  call bcast_scalar_i(matop)
  call check_lele(1,matop,n118,' 1-matop-n118 ') 
  if(.not.allocated(mapat))   stop ' Error: mapat not allocated '
  if(.not.allocated(ch_a))    stop ' Error: charge array (ch_a) not allocated '
  call check_le(1,size(ch_a),' one vs. na from size_ch ')

  if(.not.allocated(valch_a)) then
     allocate(valch_a(size(ch_a)), stat=st); call check0(st,' valch_a ')
  endif  
  call bcast_i(valch_a,size(valch_a),  0)
end subroutine bcast_atoms
