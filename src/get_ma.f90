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
subroutine get_ma(znum,ma)
  use atoms, only : mapat, matop
  implicit none
  integer znum,ma,maf
  ma = 0
  maloop : do maf=1,matop
     if(mapat(maf)==znum) ma=maf
  end do maloop
end subroutine get_ma

