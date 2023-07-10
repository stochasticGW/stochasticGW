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
subroutine calc_eorb
  use gwm, only : n, dv, nj, eorb, vks, orb_rd, read_eorb
  use simple_mpi, only : rank, bcast_scalar_r8
  implicit none
  integer st
  real*8  eorb_calc
  complex*16,allocatable :: pa(:),p1(:), tmp(:)

  rnk0 : if(rank==0) then
     allocate(pa(n),p1(n),tmp(n),stat=st); call check0(st,' pa, p1 ')
     
     pa = orb_rd
     pa=pa/sqrt(sum(abs(pa)**2)*dv)
  
     p1=0d0; call hc(pa, p1, vks, n, 1)
     eorb_calc = ( sum(conjg(pa)*p1)*dv) / (sum(conjg(pa)*pa) * dv)
     
     write(6,*)' <phi|H|phi> computed from orbital from file: ',real(eorb_calc)
     write(6,*)'                     vs. eorb read from file: ',real(eorb)
     if(abs(eorb-eorb_calc)>0.02) &
          write(6,*)' WARNING: large difference between eorb read and <phi|H|phi> '
     if(.not.read_eorb) eorb = eorb_calc
     write(6,*)     
     call debug_analyze_h(pa,vks,n)
     call flush(6)

     deallocate(pa, p1, tmp)
  end if rnk0

  call bcast_scalar_r8(eorb)

  if(nj/=0) call calc_eorbj

end subroutine calc_eorb
