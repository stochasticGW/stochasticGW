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
subroutine read_wf
  !
  ! for density and for exchange
  !
  !
  use gwm
  use simple_mpi, only : rank, nodes, sr_mpi_r8, color_size, crnk=>color_rank
  use simple_mpi, only : send_receive_r8_in_mpi, color_scatter_r8, bcast_scalar_i
  use simple_mpi, only : sync_mpi
  implicit none
  integer isf,i,nds,ir,color_ir,tag,rank0,isp,ic
  integer isfmn, isfmx, ip
  character*9             :: ch
  real*8                  :: occ
  real*8,     allocatable :: ur(:), ur0(:)
  complex*16, allocatable :: uc(:)
  rank0 = 0
  nds = max(nodes, 1)
  ip  = 440

  call alloc_ur_uc
  call alloc_map
  call read_ch

  dens0=0d0
  icfull: do ic=0,(ntop-1)/nds
     isfmn=    ic*nds+1
     isfmx=    min(isfmn+nds-1, ntop)
     ur=0d0

     isstrp : do isf=isfmn, isfmx        
        occ = occ_det_full(isf)
        !
        ! reading 
        !
        rnk0: if(rank==0) then
           select case(nsp)
           case(1)
              if(binwf) then; read(ip,  end=99  )i
              else;           read(ip,*,end=99  )i
              endif
              isp = 1
           case(2)
              if(binwf) then; read(ip, end=99  )i, isp
              else;           read(ip,*,end=99 )i, isp
              endif
              call check_lele(1,isp,2,' 1-isp-2 ')
           case default; stop ' nsp not 1 nor 2 '
           end select
           call check(isf,i,' istate vs. istate_read ')           

           call read_ur
           call nrmlz(ur)
           dens0(:,isp) = dens0(:,isp) + occ*ur(:)**2
        end if rnk0
        
        call bcast_scalar_i(isp)  
        map_sa(isf)=isp
        if(det_tddft) call read_tddft
        if(.not.rdexch) call send_u_for_xch
     end do isstrp

     isf = isfmn + rank
     if(isf.le.isfmx.and..not.rdexch) then
       occ = occ_det_full(isf)
       call exch_accum(ur, n, isf, occ)
     endif
     ur = 0d0
  end do icfull
  
  call exch_sum

  if(allocated(ur0))deallocate(ur0)
  if(allocated(ur)) deallocate(ur)
  if(allocated(uc)) deallocate(uc)
  
  return
99 write(6,*)' problem, isf=',isf,' got state ',i,' but no spin direction '
  stop
contains

  subroutine read_tddft
    implicit none
    integer  st
    integer, save ::  is=0
    real*8, allocatable :: ub(:)
    allocate(ub(n), stat=st); call check0(st,' ub ')

    if(is_start(crnk).le.isf.and.isf.le.is_end(crnk)) is=is+1

    irloop : do ir=0,nds-1
       call sync_mpi
       color_ir = mod(ir, color_size)

       belng : if(is_start(color_ir).le.isf.and.isf.le.is_end(color_ir)) then

          if(rank==0)ub=ur

          call sr_mpi_r8(ub, size(ub), 0, ir)
          isp = map_sa(isf)
          call bcast_scalar_i(isp)
          hitrank : if(rank==ir) then  ! the present rank received ub
             call check_lele(1,is, ns,' 1-is-ns ')
             phi_det(:,is)=ub
             map_sp_det(is)=isp
          end if hitrank

       endif belng
    end do irloop
    call sync_mpi
    deallocate(ub)
  end subroutine read_tddft
  
  subroutine alloc_ur_uc
    implicit none
    integer st
    if(rank==0) then; allocate(ur0(n),stat=st); call check0(st,' ur0 ')
    endif
    allocate(ur(n),stat=st); call check0(st,' ur ')
    select case(orb_det_kind)  ! real or complex
    case(1)
    case(2)
       if(rank==0) allocate(uc(n),stat=st);call check0(st,' uc ')
    case default;  stop ' orb_det_kind is not 1 neither 2 '
    end select
    ur = 0d0
  end subroutine alloc_ur_uc

  subroutine alloc_map
    implicit none
    integer st
    allocate(map_sa(ntop),stat=st);  call check0(st,' map_sa ')
    if(.not.allocated(map_sp_det)) then
       allocate(map_sp_det(ns),stat=st); call check0(st,'map_sp_det ')
    else; call check(size(map_sp_det),ns,' size_map_sp_det, ns ')
    endif
  end subroutine alloc_map

  subroutine read_ch
    implicit none
    if(rank==0) then
       if(binwf) then;  read(ip  )ch
       else;            read(ip,*)ch
       endif
       if(ch/='orbitals') then; write(6,*)" ERROR reading wf.txt. need orbitals,got ",ch; stop
       end if
    endif
  end subroutine read_ch

  subroutine read_ur
    implicit none
    if(binwf) then
       select case(orb_det_kind)
       case(1);     read(ip  )ur(:)
       case(2);     read(ip  )uc(:); call check_real(uc,size(uc)); ur = uc
       end select
    else
       select case(orb_det_kind)
       case(1);  read(ip,*)ur(:)
       case(2);  read(ip,*)uc(:);  call check_real(uc,size(uc));  ur = uc
       end select
    endif
  end subroutine read_ur

  subroutine send_u_for_xch
    implicit none
    ir=mod(isf-1,nds)
    if(ir==0.and.rank==0) ur0 = ur

    if(occ>tollocc.and.map_sa(isf)==sp0) then
       call sync_mpi
       if(rank==0.or.rank==ir) then
          tag = ir; if(ir/=0) call send_receive_r8_in_mpi(ur, size(ur), rank0, ir, tag)
       endif
    end if

    if(rank==0) ur = ur0
  end subroutine send_u_for_xch

  subroutine nrmlz(ur)
    implicit none
    real*8 ur(n)
    real*8 norm2
    norm2 = sum(abs(ur)**2)*dv 
    if(abs(norm2-1d0)>0.01) then; write(6,*)' ERROR: orbital ',isf,' norm2= ',norm2; stop;
    endif
    ur = ur / sqrt(norm2)
  end subroutine nrmlz
   
end subroutine read_wf
