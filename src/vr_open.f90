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
subroutine vr_open(i,ch)
  use gwm
  use simple_mpi, only : cs=>color_size, cr=>color_rank,rank
  implicit none
  character*3   :: ch
  integer       :: i, j, ip
  call checks

  do j=1,nvr/cs
     call ip_make(i,j,ip)
     call ip_open(i,ch,j,ip)
  end do
contains
  subroutine checks
    implicit none
    call check_le(1,nvr,    ' one-nvr       ')
    call check0(mod(nvr,cs),' nvr, cs       ')
    call check_le(nvr,nvrl, ' nvr vs. nvrl  ')
  end subroutine checks
end subroutine vr_open
