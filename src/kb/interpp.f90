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
subroutine interpp(ra, pb, ma, np)
  use kb_mod, only : rr=>rrpp, nr=>nrpp
  use kb_mod, only : ff=>phipp
  implicit none
  integer j,k,jm,ma,np
  real*8 ra, ddr, dr
  real*8 pb(1:np)
  
  if(ra<rr(1,ma)) then
     pb(1:np) = ff(1,1:np,ma) 
     return
  elseif(ra>rr(nr(ma),ma)) then
     pb(1:np) = 0d0 
     return
  end if
  
  j=1
  k=nr(ma)
  searchi : do 
     if(k-j>1) then
        jm = (j+k)/2
        if(ra<rr(jm, ma)) then
           k=jm
        else
           j=jm
        end if
     else
        exit searchi
     end if
  end do searchi
  
  ! determined the point below r. 
  ! now linear interpolation
  call check(k,j+1,' k,j+1 ')
  
  ddr = ra      -rr(j,ma)
  dr  = rr(k,ma)-rr(j,ma)
  
  pb(1:np) = ff(j,1:np,ma)  * (dr-ddr)/dr + ff(k,1:np,ma) * ddr/dr
end subroutine interpp
