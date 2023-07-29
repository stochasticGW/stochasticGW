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
subroutine pt_eta_effic
  use gwm
  use simple_mpi, only : rank, nodes, color_rank, allsum_i, bcast_i
  use simple_mpi, only : bcast_scalar_i, bcast_scalar_r8, allsum_c16, scatter_r8
  use simple_mpi, only : bcast_scalar_general_i
  implicit none
  !
  ! efficient projection of pt and eta
  !
  integer, save           :: ip=440
  integer                 :: nds,i
  integer                 :: nf
  integer                 :: isf 
  integer                 :: mpe(2)  ! 2=maxspin.
  integer                 :: mx(2)
  integer,    save        :: i0=0
  integer,    save        :: cnk=5
  real*8                  :: occ
  real*8,     allocatable :: uf(:)
  real*8,     allocatable :: ur(:)
  complex*16, allocatable :: pte1(:,:) 
  complex*16, allocatable :: pte2(:,:)

  if(eta_flg) then
     eta = zeta
  endif
  call check_flgs
  call set_nf
  call map_pt_eta
  call alloc_pteta
  call pt_eta_scatter
  call alloc_u
  call open_wf_file(ip)

  do isf=1,ntop
     call read_u     
     call proj(map_sa(isf))
  enddo
  call pt_eta_gather
  call dealloc
  call close_wf_file(ip)
contains

  subroutine check_flgs
    implicit none
    if(det_tddft)    stop ' det_tddft in pt_eta_effic    '
    if(filter_cheby) stop ' filter_cheby in pt_eta_effic '
    if(.not.det_dft) stop ' not det_dft in pt_eta_effic  '
  end subroutine check_flgs

  subroutine set_nf
    implicit none
    nds = max(1,nodes)
    nf  = (n/nds)
    if(nf*nds<n) nf=nf+1
    call check_le(n,nf*nds,' n nf*nds ')
  end subroutine set_nf
  
  subroutine map_pt_eta
    implicit none
    integer  ir, is, isp 
    mx(:) =0
    do ir=0,nds-1
       if(rank==ir) then
          do is=1,ns
             isp = map_sp_pt(is)
             mx(isp)=mx(isp)+1
          end do
          if(color_rank==eta_rank)mx(sp0)=mx(sp0)+1
       endif
    enddo
    call allsum_i(mx, size(mx))
  end subroutine map_pt_eta

  subroutine alloc_pteta
    implicit none
    integer st
    if(mx(1)>0) then; allocate(  pte1(nf, mx(1)), stat=st); call check0(st,' pte1 ')
    endif
    if(mx(2)>0) then; allocate(  pte2(nf, mx(2)), stat=st); call check0(st,' pte2 ')
    endif
  end subroutine alloc_pteta

  subroutine pt_eta_scatter
    use simple_mpi, only : scatter_c16
    implicit none
    integer ms, ir, ie, i1, i2, st, i
    integer,    allocatable :: map_s(:)
    complex*16, allocatable :: cz(:)
    allocate(cz(nf*nds), stat=st); call check0(st,' cz '); cz = 0d0
    call check_sz_pt_eta

    i1=0
    i2=0
    irloop : do ir=0,nds-1

       if(rank==ir) ms  = size(map_sp_pt)
       call bcast_scalar_general_i(ms, ir)
       msif : if(ms>0) then 
          allocate(map_s(ms), stat=st); call check0(st,' map_s ')
          if(rank==ir) map_s(:) = map_sp_pt(:)
          call bcast_i(map_s, size(map_s), ir)
       end if msif

       ptloop : do i=1,ms
          if(rank==ir) cz(1:n)=pt(:,i)
          select case(map_s(i))
          case(1);  i1=i1+1; call scatter_c16(cz, pte1(:,i1), nf, ir)
          case(2);  i2=i2+1; call scatter_c16(cz, pte2(:,i2), nf, ir)
          case default; stop ' map_s(i) '
          end select
       enddo ptloop

       ie = 0; if(rank==ir.and.color_rank==eta_rank) ie=1
       call bcast_scalar_general_i(ie, ir)
       eif : if(ie==1) then
          
          if(rank==ir) cz(1:n) = eta 
          select case(sp0)
          case(1); i1=i1+1; call scatter_c16(cz, pte1(:,i1), nf, ir)
          case(2); i2=i2+1; call scatter_c16(cz, pte2(:,i2), nf, ir)
          end select
       end if eif

       if(allocated(map_s)) deallocate(map_s)
    end do irloop
    if(allocated(cz)) deallocate(cz)
  end subroutine pt_eta_scatter

  subroutine check_sz_pt_eta
    implicit none
    call check(   n,size(eta),' sz_eta ')
    call check_le(n,nf*nds,   ' sz_cz  ')
    if(ns>0) then
       call check(size(pt,1),n, ' sz_pt_1 ')
       call check(size(pt,2),ns,' sz_pt_2 ')
    endif
  end subroutine check_sz_pt_eta

  subroutine alloc_u
    implicit none
    integer st
    allocate(uf(nf), stat=st); call check0(st,' wf_ch ')
    if(rank==i0) then; allocate(ur(nf*nds),stat=st); call check0(st,' wf_ur ')
    else;              allocate(ur(1),stat=st);      call check0(st,' wf_u1 ')
    endif
  end subroutine alloc_u

  subroutine read_u
    implicit none
    integer isp

    !ir=mod(isf-1,nds)
    occ = occ_det_full(isf)

    if(rank==i0) then
       ur = 0d0
       isp = map_sa(isf)
       call read_wf1(ur, n, ip, isf, isp)
       ur = ur*sqrt(occ*nsp/2d0)
    endif
    call bcast_scalar_i( isp)
    call bcast_scalar_r8(occ)
    
    if(occ>tollocc) then
       call scatter_r8(ur, uf, nf, i0)
    end if
  end subroutine read_u
     

  subroutine proj(isp)  ! sqrt(occ*nsp/2d0) absorbed already to u.
    implicit none
    integer isp, m, st
    integer jb, jt, j
    complex*16, allocatable :: ov(:)

    m = mx(isp)             ! mx= # states of pte of given spin.
    if(m==0.or.occ<tollocc) return    
     allocate(ov(m), stat=st); call check0(st,' ov ')       

    do jb=1,m,cnk
       jt=min(m, jb+cnk-1)
       if(isp==1) ov(jb:jt) = dv*matmul(uf,pte1(:,jb:jt))
       if(isp==2) ov(jb:jt) = dv*matmul(uf,pte2(:,jb:jt))
    enddo
    call allsum_c16(ov,size(ov))
    do j=1,m
       if(isp==1) pte1(:,j) = pte1(:,j)-ov(j)*uf
       if(isp==2) pte2(:,j) = pte2(:,j)-ov(j)*uf
    enddo

    deallocate(ov)
  end subroutine proj

  subroutine pt_eta_gather
    use simple_mpi, only : gather_c16
    implicit none
    integer ms, ir, ie, i1, i2, st, i
    integer,    allocatable :: map_s(:)
    complex*16, allocatable :: cz(:)
    allocate(cz(nf*nds), stat=st); call check0(st,' cz '); cz = 0d0
    call check_sz_pt_eta

    i1=0
    i2=0
    irloop : do ir=0,nds-1
       if(rank==ir) ms  = size(map_sp_pt)
       call bcast_scalar_general_i(ms, ir)
       msif : if(ms>0) then 
          allocate(map_s(ms),stat=st); call check0(st,' map_s ')
          if(rank==ir)map_s(:) = map_sp_pt(:)
          call bcast_i(map_s, size(map_s), ir)
       end if msif

       ptloop : do i=1,ms
          select case(map_s(i))
          case(1);  i1=i1+1; call gather_c16(pte1(:,i1), cz, nf, ir)
          case(2);  i2=i2+1; call gather_c16(pte2(:,i2), cz, nf, ir)
          case default; stop ' map_s(i) '
          end select
          if(rank==ir) then
             pt(:,i)=pt(:,i)-cz(1:n)
             call nrmlz_pt(pt(:,i),normpt(i))
          endif
       enddo ptloop

       ie = 0; if(rank==ir.and.color_rank==eta_rank) ie=1
       call bcast_scalar_general_i(ie, ir)
       eif : if(ie==1) then
          select case(sp0)
          case(1);  i1=i1+1; call gather_c16(pte1(:,i1), cz, nf, ir)
          case(2);  i2=i2+1; call gather_c16(pte2(:,i2), cz, nf, ir)
          end select
          if(rank==ir) eta=eta-cz(1:n)
       end if eif

       if(allocated(map_s)) deallocate(map_s)
    end do irloop
    if(allocated(cz))    deallocate(cz)
  end subroutine pt_eta_gather

  subroutine nrmlz_pt(pt1,normpt1)
    implicit none
    complex*16 pt1(n)
    real*8 normpt1(1)
    normpt1(1) = sqrt(sum(abs(pt1)**2)*dv)
    pt1 = pt1/normpt1(1)
  end subroutine nrmlz_pt

  subroutine dealloc
    implicit none
    if(allocated(uf))deallocate(uf)
    if(allocated(ur))deallocate(ur)
    if(allocated(pte1)) deallocate(pte1)
    if(allocated(pte2)) deallocate(pte2)
  end subroutine dealloc
end subroutine pt_eta_effic
