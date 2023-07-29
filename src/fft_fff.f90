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
subroutine fft_fff(nx, ny, nz, ck, cr, rank) ! emulates cuda call for cpu using fftw
  implicit none
  integer nx, ny, nz, rank
  complex*16 cr(nx,ny,nz)
  complex*16 ck(nx,ny,nz)
  call fft3d_forward(nx, ny, nz, ck, cr)
end subroutine fft_fff
