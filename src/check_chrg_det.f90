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
subroutine check_chrg
  use simple_mpi, only : cs=>color_size
  use simple_mpi, only : crnk => color_rank
  use simple_mpi, only : rank
  use simple_mpi, only : nodes
  use simple_mpi, only : color_reduce_sum_r8, color_gather_r8
  use simple_mpi, only : sync_mpi
  use gwm,        only : ptheader_flg, ptheader_rank, det_tddft, det_dft
  use gwm,        only : occall => occ_det_full, map_sp_det, map_sa
  use gwm,        only : nsp
  use gwm,        only : n
  use gwm,        only : ns
  use gwm,        only : dv
  use gwm,        only : rnel
  use gwm,        only : tollocc
  use gwm,        only : phi_det
  use gwm,        only : occ_det
  use gwm,        only : dens0
  implicit none
  integer  ir, i, sp, st, nds
  real*8 :: chrg(1)
  real*8, allocatable :: dens_b(:,:)

  nds = max(1,nodes)

  if(.not.det_dft) stop ' ERROR: check_chrg should only be used if deterministic dft is used'

  !
  ! check charge if det. tddft
  !
  ifdet_tddft : if(det_tddft) then 
     
     !
     ! diff. phi on separate colors.
     !
     allocate(dens_b(n, nsp), stat=st); if(st/=0) stop ' dens_b '
     
     dens_b = 0d0
     
     do i=1,ns
        sp = map_sp_det(i)
        if(occ_det(i)>tollocc) then
           dens_b(:,sp) = dens_b(:,sp) + occ_det(i)* phi_det(:,i)**2
        endif
     enddo
     
     ! checking charge 
     
     call color_reduce_sum_r8(dens_b, size(dens_b), ptheader_rank)
     chrg(1) = sum(dens_b)*dv

     if(crnk==ptheader_rank.and.abs(chrg(1)-rnel)>1d-6) then
        write(6,*)" ERROR: charge, rnel ",chrg,rnel
        stop
     endif
     
     deallocate(dens_b)
  end if ifdet_tddft
  
end subroutine check_chrg
