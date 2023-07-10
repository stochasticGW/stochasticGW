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
subroutine read_wf_det_calc_dens0(ip)
  use gwm
  use gwm, only : occall => occ_det_full
  implicit none
  integer :: st, ip

  call check_param
  call read_evlocc(ip)
  if(det_tddft) call distribute_occ_det
  if(det_tddft) call allocate_phi_det
  call read_wf
  call distribute_dens0
  call check_chrg 

contains
  subroutine check_param
    use simple_mpi, only : rank, bcast_scalar_i
    implicit none
    real*8 rnon_full
    rnk0 : if(rank==0) then
       det_w : if(det_tddft) then
          if(2d0/dble(nsp)*ns_blk<rnel) then
             write(6,*) " ERROR: At least ",ceiling(rnel*nsp/2d0-1d-8),&
                        " states are required, but #states=",ns_blk,"!"
             stop
          end if
       end if det_w
    end if rnk0
  end subroutine check_param
end subroutine read_wf_det_calc_dens0

subroutine allocate_phi_det
  use gwm
  implicit none
  integer st
  if(ns>0) then! deterministic tddft, w.f. distributed per cores.
     allocate(phi_det(n,ns),stat=st); if(st/=0) stop "p<hi_det allocation failed!"
  endif  
end subroutine allocate_phi_det

subroutine distribute_dens0
  use simple_mpi, only : bcast_r8
  use gwm, only : dens0
  implicit none
  call bcast_r8(dens0, size(dens0), 0)
end subroutine distribute_dens0

