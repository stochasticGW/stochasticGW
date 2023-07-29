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
subroutine prep_norm_etaxi
  use simple_mpi, only : rank, color_bcast_r8, color_bcast_c16
  use gwm
  implicit none
  integer, parameter :: n2=2
  real*8 :: norm(n2)

  if(eta_flg) then
     normeta  = sqrt(sum(conjg(eta)*eta)*dv)
     normxi   = sqrt(sum(conjg(xi )*xi )*dv)
     etaxi(:,1)= eta/normeta
     etaxi(:,2)= xi /normxi
     norm= (/ normeta, normxi /)
  endif
  
  call color_bcast_r8( norm,  size(norm), eta_rank)
  call color_bcast_c16(etaxi, size(etaxi), eta_rank)

  normeta=norm(1)
  normxi =norm(2)
end subroutine prep_norm_etaxi

