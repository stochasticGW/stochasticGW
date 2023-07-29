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
real*8 function filter_trace(ev, eorb, de0)  ! change if needed
  implicit none
  integer, save :: i1=1
  real*8,  save :: pi
  real*8        :: ev, eorb, de0
  real*8        :: x
  
  if(i1==1) then
     i1=-1
     pi = dacos(-1d0)
  end if

  x = (ev-eorb)/de0

  if(abs(x)>pi/2d0) then
     filter_trace = 0d0
  else
     filter_trace = cos(x*pi/2d0)
  end if
end function filter_trace
