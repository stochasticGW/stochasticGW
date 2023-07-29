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
subroutine dr_make(ra,dr,mr)  ! for exponential grid.
  implicit none
  integer ir, mr
  real*8  ra(mr), dr(mr), ds
  !
  ! check that the grid is exponential
  !
  do ir=2,mr-1
     if(abs(ra(ir)/ra(ir-1) - ra(ir+1)/ra(ir))>1d-5) then
        write(6,*)' ra grid not logarithmic, ir,',ir, &
             ' ra(ir)/ra(ir-1) ',ra(ir)/ra(ir-1),&
             ' ra(ir+1)/ra(ir) ',ra(ir+1)/ra(ir)
        stop
     endif
  enddo
  !
  ! r = exp( s),  s linear variable
  ! dr = exp(s) ds = r ds. 
  ! s = log(r) --> ds = log(r(i))-log(r(i-1))
  
  ds = log(ra(2)/ra(1))
  dr = ds* ra
end subroutine dr_make
