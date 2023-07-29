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
subroutine vr_to_vo_gam
  use gwm, only : nw
  use gwm, only : ft
  use gwm, only : nt
  use gwm, only : ngam
  implicit none
  integer i, ih, many, iw, st
  integer, parameter :: cnk = 200
  complex*16, allocatable :: fw(:,:)
  
  call check_le(nt*2+2,  nw,' nt2nw        ')

  do i=1,ngam,cnk

     ih=min(i-1+cnk,ngam)
     many = ih+1-i

     allocate(fw(0:nw-1,i:ih), stat=st); call check0(st,' fw ')

     fw = 0d0;       
     fw(0:nt,:)   = ft(0:nt,i:ih)       
     fw(0,:)      = fw(0,:)/2d0                                             ! theta(t=0)= 1/2

     fw = conjg(fw);        
     call fftw_multiple(fw,nw,many);     
     fw = conjg(fw)/dble(nw)                                                ! fft with exp(iwt)

     fw(nw/2:nw-1,:) = conjg(fw( nw/2:nw-1,:))                              ! ret  to non ret
     call fftw_multiple(fw,nw,many)                                         ! fft with exp(-iwt)
     ft( 0:nt,i:ih) = fw(0:nt,:)                 

     deallocate(fw)
  enddo
end subroutine vr_to_vo_gam
