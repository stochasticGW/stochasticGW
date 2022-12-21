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
subroutine vr_close(i)
  use gwm
  use simple_mpi, only : cs=>color_size, cr=>color_rank,rank
  implicit none
  integer       :: i, j, ip
  do j=1,nvr/cs
     call ip_make(i,j,ip)
     close(ip)
  end do
end subroutine vr_close

subroutine vr_close_rmv(i)
  use gwm
  use simple_mpi, only : cs=>color_size, cr=>color_rank,rank
  implicit none
  integer       :: i, j, ip

  do j=1,nvr/cs
     call ip_make(i,j,ip)
     close(ip,status='delete')
  end do
end subroutine vr_close_rmv
