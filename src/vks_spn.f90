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
subroutine dens_to_v_spn(dens,vks,n, nsp)
  use gwm, only : dv
  use gwm, only : rnel_orig
  use simple_mpi, only : rank
  implicit none
  integer n, nsp, st, sp
  real*8  dens(n, nsp)
  real*8   vks(n, nsp)
  real*8, allocatable :: dens_nrm(:, :)

  allocate(dens_nrm(n, nsp), stat=st);  call check0(st,' dens_nrm ')
  !  write(60000+rank,*)' rnel_orig ',rnel_orig

  dens_nrm = dens*(rnel_orig/sum(dens*dv))
  !do sp=1,nsp
  !   write(60000+rank,*)' sp, dens_nrm ',sp, sum(dens_nrm(:,sp))*dv; call flush(rank)
  !enddo
  call sub_vks_lda_spn(dens_nrm, vks, n , nsp)

  deallocate(dens_nrm)
end subroutine dens_to_v_spn

subroutine sub_vks_lda_spn(dens,vks,n, nsp)   ! vh+vx+vcn+vloc
  use gwm,        only : vloc_tot
  use simple_mpi, only : rank, nodes
  implicit none
  integer n,nsp,i,st,sp
  real*8 dens(n,nsp), vks(n,nsp)
  integer, save :: count=0
  real*8, allocatable, dimension(:) :: vxc(:,:)

  call vh_spn(dens, vks, n, nsp)

  allocate( vxc(n,nsp), stat=st);  call check0(st, ' vx vcn ')
  call get_vxc_spn( dens, n, nsp, vxc) ! note order

  ! hartree, xc, 1body_loc
  !do sp=1,nsp
  !  write(60000+rank,*)' sp dens, vh, vxc, vloc_tot ',sp, &
  !  sum(abs(dens(:,sp))),sum(abs(vks(:,sp))),sum(abs(vxc(:,sp))),sum(abs(vloc_tot))
  !enddo
  !call flush(60000+rank)

  vks = vks + vxc 
  do sp=1,nsp
     vks(:,sp) = vks(:,sp)+ vloc_tot
  enddo

  deallocate(vxc)
end subroutine sub_vks_lda_spn
