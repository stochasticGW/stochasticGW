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
!---------------------------------
!
!  Checks that a vector is real.
!
!--------------------------------

subroutine check_c_real(cvec, mm, aachar)
  implicit none
  
  integer mm
  complex*16    cvec(mm)
  character(*)  aachar
  real*8        temp
  
  temp = sum(abs(aimag(cvec)))/  &
       (sum(abs(cvec))+1.d-2)
  if(temp>1.d-6.and.sum(abs(cvec))>1.d-6) then
     write(6,*)'  problem in     ',aachar
     write(6,*)'  imag rel. part ', sum(abs(aimag(cvec))),sum(abs(cvec))
     stop
  end if
end subroutine check_c_real

