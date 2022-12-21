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
subroutine addvnl_c(p,hp,n,ms,cnst)
  use gwm, only : na
  implicit none
  integer n, ms
  real*8 cnst
  complex*16 p(n,ms), hp(n,ms)

  integer i, st, ia, ilm, ma, igg, l,is
  complex*16 cs

  if(ms/=1) then
    write(6,*) "ERROR: ONLY SPIN UNPOLARIZED CALCULATION IMPLEMENTED!"
    call flush(6)
  end if
  do ia=1,na
    call addvnl_c_kb(p,hp,n,ms,cnst,ia)
  enddo
end subroutine addvnl_c
