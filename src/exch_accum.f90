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
subroutine exch_accum(ur, nin, isf, occ)
  use simple_mpi, only : rank
  use gwm
  implicit none
  integer nin, st, isf, ie, ij
  real*8  ur(nin), occ
  real*8, allocatable :: rh(:), rh1(:), pot(:)

  call check(n,nin,' n nin            ')
  call check(ne,1, ' ne should be one ')
  if(trace) stop   ' should modify xch_accumulate  if using trace '

  if(occ>tollocc.and.map_sa(isf)==sp0) then
     allocate(rh(n), pot(n), stat=st); call check0(st,' rh, pot ')
     ie=1
     rh(:) = ur(:)*orb_rd(:)*sqrt(nsp*occ/2d0)
     call vh_sub_exch( rh, pot, n)
     exce(ie) = exce(ie) - dv*sum(rh*pot)

     if(nj/=0) then
        allocate(rh1(n), stat=st); call check0(st,' rh1 ')
        do ij=1,nj
           rh1(:)=ur(:)*rdorbj(:,ij)*sqrt(nsp*occ/2d0)
           call vh_sub_exch( rh1, pot, n)
           excej(ie,ij) = excej(ie,ij) - dv*sum(rh*pot) !!! PT sum rh1 here ???
        enddo
        deallocate(rh1)
     endif

     deallocate(rh, pot)
  endif
end subroutine exch_accum
     
  
