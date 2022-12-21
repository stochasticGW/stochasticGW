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
subroutine get_vxc_spn(dens, n, nsp, vxc)  ! include core_correction
  use gwm, only : xc_type, vxc0, rho_core, dim_periodic
  use simple_mpi, only : rank
  implicit none
  integer n,nsp,sp, st
  integer, save :: i1=1
  real*8  dens(n, nsp)
  real*8   vxc(n, nsp)
  real*8   vxc_avg
  real*8, allocatable :: dens_tot(:,:)
  if(xc_type==1) then
     vxc = vxc0
  else
     allocate(dens_tot(n, nsp),stat=st); if(st/=0) stop ' dens_tot  '
     do sp=1,nsp
        dens_tot(:,sp) = dens(:,sp) + rho_core(:)/dble(nsp)
     enddo
     call vxc_spn_choose(dens_tot, n, nsp, vxc);
     deallocate(dens_tot)
  end if

!  if(dim_periodic.gt.0) then
!     vxc_avg = sum(vxc)/dble(size(vxc))  
!     vxc = vxc-vxc_avg
!
!     if(i1==1.and.rank==0)then
!        write(6,*)' vxc_avg (removed) ',vxc_avg
!     end if
 ! end if
  i1=-1
end subroutine get_vxc_spn

subroutine vxc_spn_choose(dens, n, nspin, vxc)
!  use gwm, only : functional 
  implicit none
  integer n, nspin
  real*8 dens(n, nspin), vxc(n,nspin)

  call vxc_libxc(dens, vxc, n, nspin)
!  select case(functional)
!  case('pwlda','lda')
!     call vxc_spn_pwlda(dens, vxc, n, nspin)
!  case('pbe')
!     call check(nspin,1,' nspin-one for pbe ')
!     call vxc_pbe_new(dens, vxc, n, nspin) 
  !case (number)
  !call libxc with the code
!  case default
!     write(6,*)' ERROR, functional should be pwlda (or lda) or pbe, not  ',functional
!     stop
!  end select
end subroutine vxc_spn_choose

subroutine vxc_spn_pwlda(dens, vxc, n, nsp)  ! include core_correction
  use gwm, only : dv
  implicit none
  logical, parameter :: debug=.false. 
  integer n,nsp,sp, st
  real*8  dens(n, nsp)
  real*8   vxc(n, nsp)
  real*8, allocatable :: vx(:,:), vc(:,:)

  allocate(vx(n, nsp), vc(n, nsp), stat=st); if(st/=0) stop ' vx vc  '
  vx = 0d0
  vc = 0d0
  call get_vx_spn_lda(    dens, n, nsp, vx);
  call get_vcn_spn_lda(   dens, n, nsp, vc);

  vxc = vx + vc

  if(debug) then
     !vxc = vx
     write(6,*)' integral vx dr ',sum(vx)*dv
     write(6,*)' integral vc dr ',sum(vc)*dv
  endif

  deallocate(vx,vc)
end subroutine vxc_spn_pwlda

