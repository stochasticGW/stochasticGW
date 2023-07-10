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
subroutine vh_spn(dens, vks, n, nsp)
  implicit none
  integer n, nsp, st
  real*8 dens(n, nsp)
  real*8 vks( n, nsp)
  select case(nsp)
  case(1)
     call  vh_sub( dens, vks, n)
  case(2)
     vks(:,1) = dens(:,1)+ dens(:,2)  ! trick to represnet the total dens. w/o add. storage.
     call  vh_sub( vks(:,1), vks(:,2), n)
     vks(:,1)=vks(:,2)
  case default
     stop ' nsp is neither 1 nor 2 '
  end select
end subroutine vh_spn
