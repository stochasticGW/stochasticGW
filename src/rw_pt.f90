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
subroutine read_pt
  use gwm, only : pt
  implicit none
  call open_pt_file_old
  read( 320000)pt
  close(320000)
end subroutine read_pt

subroutine write_pt
  use gwm, only : pt
  implicit none
  call open_pt_file_replace
  write(320000)pt
  close(320000)
end subroutine write_pt
  
