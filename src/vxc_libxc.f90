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
subroutine vxc_libxc(dens, vxc, n, nspin) ! erase
  implicit none
  integer n, nspin, st
  real*8  dens(n), vxc(n)
  real*8, allocatable :: exc_pbe(:)
  allocate(exc_pbe(n), stat=st); call check0(st,' exc_pbe ')
  call call_libxc(dens, n, nspin, vxc, exc_pbe)
  deallocate(exc_pbe)
end subroutine vxc_libxc
