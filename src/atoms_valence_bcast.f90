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
subroutine bcast_valence
  use atoms, only : valch_a, ch_a, valence_list
  use gwm,   only : na
  use simple_mpi, only : bcast_i
  implicit none
  integer i
  call bcast_i(valence_list,size(valence_list),0)
  do i=1,size(valch_a)
     valch_a(i) = valence_list(ch_a(i))
  enddo
end subroutine bcast_valence

subroutine bcast_valence_atoms
  use atoms, only : valch_a, ch_a, valence_list, matop, mapat, n118
  use simple_mpi, only : bcast_i, bcast_scalar_i
  implicit none
  integer :: i,st

  call bcast_scalar_i(matop)
  call check_lele(1,matop,n118,' 1-matop-n118 ')
  if(.not.allocated(mapat))   stop ' Error: mapat not allocated '
  if(.not.allocated(ch_a))    stop ' Error: charge array (ch_a) not allocated '
  call check_le(1,size(ch_a),' one vs. na from size_ch ')
  call bcast_i(valence_list,size(valence_list),0)
  allocate(valch_a(size(ch_a)), stat=st); call check0(st,' valch_a ')
  do i=1,size(valch_a)
     valch_a(i) = valence_list(ch_a(i))
  enddo

end subroutine bcast_valence_atoms

