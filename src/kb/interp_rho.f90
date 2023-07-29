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
subroutine interp_rho_add(ri, vi, ra, va, nr)
  implicit none
  integer j, k, jm, l, nr
  real*8 vi, ra(nr), va(nr)
  real*8 ri, ddr, dr

  if(ri<ra(1)) then
     vi=vi+ va(1)
     return
  elseif(ri>ra(nr)) then
     vi=vi + va(nr)
     return
  end if
  
  j=1
  k=nr
  searchi : do 
     if(k-j>1) then
        jm = (j+k)/2
        if(ri<ra(jm)) then
           k=jm
        else
           j=jm
        end if
     else
        exit searchi
     end if
  end do searchi
  call check(k-1,j,' k-1, j ')
  
  ! determined the point below r. 
  ! now linear interpolation
  
  ddr = ri    -ra(j)
  dr  = ra(k) -ra(j)
  
  vi=vi + va(j)  * (dr-ddr)/dr + va(k) * ddr/dr
end subroutine interp_rho_add
