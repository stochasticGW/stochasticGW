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
subroutine vh_sub_exch(dens, vh, n)
  use gwm,        only : vk, vk_exch, vk_old
  use gwm,        only : dim_periodic
  implicit none
  integer n
  real*8 dens(n), vh(n)
  select case(dim_periodic)
  case(0) ! non periodic
     call vh_sub(dens,vh,n)
     return
  case(2,3)
     vk = vk_exch
     call vh_sub(dens,vh,n)
     vk = vk_old
  case default
     stop ' dimp problem '
  end select

end subroutine vh_sub_exch

subroutine cvh_sub_exch(dens, vh, n)
  use gwm,        only : vk, vk_exch, vk_old
  use gwm,        only : dim_periodic
  implicit none
  integer n
  complex*16 dens(n), vh(n)
  select case(dim_periodic)
  case(0) ! non periodic
     call cvh_sub(dens,vh,n)
     return
  case(2,3)
     vk = vk_exch
     call cvh_sub(dens, vh, n)
     vk = vk_old
  case default
     stop ' dimp problem '
  end select
end subroutine cvh_sub_exch

