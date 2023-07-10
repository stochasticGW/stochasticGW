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
subroutine fftw_good(fw,nw)
  implicit none
  integer nw, many
  complex*16 fw(nw)
  many = 1
  call fftw_multiple(fw,nw,many)
end subroutine fftw_good
