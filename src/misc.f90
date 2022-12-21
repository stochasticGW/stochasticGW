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
subroutine check_real_1(c1,a)
  implicit none
  complex*16 c1
  character(*) a
  if(abs(aimag(c1))>1e-6) then
     write(6,*)' problem in ',a,' = ',c1
     stop
  end if
end subroutine check_real_1
