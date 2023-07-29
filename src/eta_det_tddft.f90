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
subroutine eta_det_tddft
  use simple_mpi, only : rank, color_rank
  use gwm
  implicit none
  integer                 :: st
  complex*16, allocatable :: cz(:)

  allocate(cz(n), stat=st); call check0(st,' cz ')
  cz = 0d0
  if(color_rank==eta_rank)  cz = zeta

  call project_eta_det(cz,n,sp0)

  if(color_rank==eta_rank)  then
     eta = cz 
     call check_real(eta, size(eta))
  end if

  deallocate(cz)
end subroutine eta_det_tddft

subroutine project_eta_det(czeta, n, sp)
  use gwm, only : phi_det, ns, dv, eta_rank, map_sp_det, occ_det, tollocc, map_sp_det, nsp
  use simple_mpi, only : color_reduce_sum_c16, color_bcast_c16, rank, color_rank
  implicit none
  integer n, st, i, sp
  real*8     occ
  complex*16 czeta(n)
  complex*16 cov
  complex*16, allocatable :: cu(:)
  call color_bcast_c16(czeta, n, eta_rank)
  allocate(cu(n), stat=st); call check0(st,' cu ')
  cu = 0d0
  do i=1,ns
     occ = occ_det(i)
     if(sp==map_sp_det(i).and.occ>tollocc) then
        cov =   dv*sum(phi_det(:,i)*czeta) * sqrt(occ*nsp/2d0)
        cu  = cu + cov*phi_det(:,i)
     endif
  enddo
  call color_reduce_sum_c16(cu, n, eta_rank)
  if(color_rank==eta_rank) czeta = cu
  deallocate(cu)
end subroutine project_eta_det

