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
subroutine calc_eorbj
  use gwm, only : n, dv, eorbj, vks, rdorbj, nj, orbj_indx
  use simple_mpi, only : rank, bcast_r8
  implicit none
  integer st, ij, isp
  real*8  eorb_calc
  complex*16,allocatable :: pa(:),p1(:), tmp(:)

  allocate(eorbj(nj),stat=st)

  rnk0 : if(rank==0) then
     allocate(pa(n),p1(n),tmp(n),stat=st); call check0(st,' pa, p1 ')

     do ij=1,nj

       pa = rdorbj(:,ij)
       pa=pa/sqrt(sum(abs(pa)**2)*dv)

       p1=0d0; call hc(pa, p1, vks, n, 1)
       eorb_calc = ( sum(conjg(pa)*p1)*dv) / (sum(conjg(pa)*pa) * dv)

       write(6,'(X,A,I0,A,E14.7)')' KS energy computed for orbital ',&
             orbj_indx(ij),' (<phi|H|phi>): ',real(eorb_calc)
       call flush(6)

       eorbj(ij) = eorb_calc
     enddo

     deallocate(pa, p1, tmp)
  end if rnk0

  call bcast_r8(eorbj,size(eorbj),0)
end subroutine calc_eorbj
