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
subroutine gw_core
  use simple_mpi, only : color_reduce_sum_c16, rank
  use gwm
#if GPU_ENABLED
  use device_mem_module
#endif
  implicit none
  integer :: st
  integer, save :: n2=2
  real*8  :: t0, t1

  call alloc_ge_gge
  if(eta_flg) call set_seed_zeta()
  call set_ear
  call set_ge
  allocate(g(n),stat=st); call check0(st,' g allocation problem ')
  call reset_g
  call calc_gge
  call write_ge
  call allocate_calc_del_zeta  ! del, zeta: calc. in eta_rank, distributed.

  !
  ! all cores in a given color will project for psi representing W
  ! one core per color for eta(and xi) representing G(t<0).
  !  
  allocate(eta(n),  stat=st); call check0(st,' eta ')
  allocate(pt(n,ns),stat=st); call check0(st,' pt  ')

  call pt_eta
  call write_pt

  !
  ! finished the projection stage
  !
  !
  ! exchange from pt 
  !
  call xc_expect
  if(gamflg==1.or.gamflg==2) call gam_prep

  !
  ! now eta_xi
  !
  if(eta_flg) then
     allocate(xi(n), stat=st); call check0(st,' xi ')
     xi=zeta-eta
  end if
  if(allocated(zeta)) deallocate(zeta)

  allocate(etaxi(n, n2), stat=st); call check0(st,' etaxi ')
  call prep_norm_etaxi
  if(allocated(eta)) deallocate(eta)
  if(allocated(xi))  deallocate(xi)

  !
  ! now propagation
  !
#if GPU_ENABLED
  if (.not.usegpu) then
     call vo_make ! CPU version
  else
     call init_propdtpt_device
     call init_vomake_device
     call vo_make_gpu
  endif
#else
  call vo_make
#endif
  if(allocated(vh0))   deallocate(vh0)
  if(allocated(pt))    deallocate(pt)
  if(allocated(del))   deallocate(del)

  !
  !
  ! now ct=<phi eta_or_xi(t) | W(t) | phi zeta >
  ! Need to sum within color
  !
  call prnt(1,'DEBUG: pre G0 propagation')
  allocate(vo(n,nspv),stat=st); call check0(st,' vo ')
#if GPU_ENABLED
  if (.not.usegpu) then
     call make_ct ! CPU version
  else
     call make_ct_gpu
     call flush_vomake_device
     call flush_propdtpt_device
     call flush_device
  endif
#else
  call make_ct
#endif
  deallocate(vo)

  if(gamflg.eq.1.or.gamflg.eq.2) then ! no need in direct
    call color_reduce_sum_c16(ct, size(ct), ptheader_rank)
    if(nj>0) call color_reduce_sum_c16(ctj, size(ctj), ptheader_rank)
  endif

  call gwplt

  if(allocated(gam))    deallocate(gam)
  if(allocated(etaxi))  deallocate(etaxi)
  deallocate(gge, ge, g, stat=st); call check0(st,' gge dealloc. ')
  if(nj>0) deallocate(gej, stat=st)

contains

  subroutine make_ct
    implicit none
    integer :: it, ij
    integer :: map_sp_etaxi(2)
    real*8  :: wg

    map_sp_etaxi(:)=sp0

    ct=0d0
    if(nj>0) ctj=0d0
    do it=0,nt
       call vo1(it)
       wg = wgf(it); if(it==0)wg=wg/2d0
       call prepcv
       ct( it,:) = ct( it,:) + cvxi*normxi*wg
       ct(-it,:) = ct(-it,:) - cveta*normeta*wg

       if (nj>0) then
          do ij=1,nj
             ctj( it,:,ij) = ctj( it,:,ij) + cvxij(:,ij)*normxi*wg
             ctj(-it,:,ij) = ctj(-it,:,ij) - cvetaj(:,ij)*normeta*wg
          enddo
       endif

       call propdt(etaxi,n,n2,map_sp_etaxi,0)

       if(mod(it,100)==0) then
         call prnt(1,'DEBUG: G0 prop after : 0,100,200,... steps')
       end if
    enddo
  end subroutine make_ct

  subroutine alloc_ge_gge
    implicit none
    if(.not.allocated(ge)) then
       allocate(ge(n,ne), stat=st); if(st/=0) stop ' ge  allocation problem '
    endif
    if(.not.allocated(gge)) then
       allocate(gge( ne), stat=st); if(st/=0) stop ' gge allocation problem '
    endif
    if(nj>0) then
       if(.not.allocated(gej)) then
          allocate(gej(n,ne,nj), stat=st); if(st/=0) stop ' gej  allocation problem '
       endif
    endif
  end subroutine alloc_ge_gge

  subroutine prepcv
    use gwm
    implicit none
    integer ie, ij
    do ie=1,ne
       cveta(ie) = sum(conjg(etaxi(:,1))*vo(:,sp0)*ge(:,ie))*dv
       cvxi( ie) = sum(      etaxi(:,2) *vo(:,sp0)*ge(:,ie))*dv

       if(nj>0) then
          do ij=1,nj
            cvetaj(ie,ij) = sum(conjg(etaxi(:,1))*vo(:,sp0)*gej(:,ie,ij))*dv
            cvxij( ie,ij) = sum(      etaxi(:,2) *vo(:,sp0)*gej(:,ie,ij))*dv
          enddo
       endif
    enddo
  end subroutine prepcv
  
  subroutine reset_g
    implicit none
    integer ie
    g = 0d0
    do ie=1,ne
       g = g + ge(:,ie)
    enddo
  end subroutine reset_g

  subroutine calc_gge
    implicit none
    integer ie
    do ie=1,ne
       gge(ie) = sum(ge(:,ie)**2.d0)*dv
    enddo
  end subroutine calc_gge

  subroutine set_ear
    implicit none
    if(ne/=1) then
       write(6,*)' problem, not set for   ne ',ne,' stopping '
       stop
    end if
    ear(:)= eorb
  end subroutine set_ear

end subroutine gw_core

