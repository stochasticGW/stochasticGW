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
subroutine read_ge
  use gwm, only : ge
  implicit none
  call open_ge_file_old
  read( 310000)ge
  close(310000)
end subroutine read_ge

subroutine write_ge
  use gwm, only : ge
  implicit none
  call open_ge_file_replace
  write(310000)ge
  close(310000)
end subroutine write_ge
