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
subroutine color_flgs
  use simple_mpi, only : rank, nodes, bcast_scalar_i
  use simple_mpi, only : cs=>color_size, color_rank
  use simple_mpi, only : color_allsum_scalar_i, sync_mpi
  use gwm,        only : cs_pt=> color_size_pt
  use gwm,        only : ns,      ns_blk,   is_map, nsp
  use gwm,        only : ngam,    ngam_blk, igam_map, ns_each
  use gwm,        only : ptheader_flg, ptheader_rank,  eta_flg, eta_rank
  use gwm,        only : det_tddft, rnel, det_dft
  use gwm, only : is_start, is_end
 
  implicit none
  integer nds, iold, is_top, igam_top, st
  integer ns_blk_p1, ns_blk_tmp, icr, is
  nds = max(nodes,1)

  !
  ! now assign ptheader and eta_rank
  ! 

  ptheader_rank = 0
  eta_rank      = cs-1
  
  if(color_rank==ptheader_rank) then
     ptheader_flg = .true.
  else
     ptheader_flg = .false.
  endif

  if(color_rank==eta_rank) then
     eta_flg = .true.
  else
     eta_flg = .false.
  end if

  ! eta is done in last color_rank, so put less states ("ns") there

  call check_le(1,cs,' one-cs')
  cs_pt = cs
  ns_blk_p1 = ns_blk + 1
  
  is_top = ns_blk_p1/cs_pt
  if(mod(ns_blk_p1, cs_pt)/=0) is_top = is_top+1

  if(ngam_blk>0) then
     igam_top = ngam_blk/cs_pt
     if(mod(ngam_blk, cs_pt)/=0) igam_top = igam_top + 1
     allocate(igam_map(igam_top, 0:cs_pt-1), stat=st); call check0(st,' igam_map ')
     call calc_map_n(igam_map, igam_top, ngam_blk, ngam)
  endif

  allocate(is_map(  is_top,   0:cs_pt-1), stat=st); call check0(st,' is_map ')
  call calc_map_n(is_map,   is_top,   ns_blk_p1,   ns)
  crloop : do icr=cs_pt-1,0,-1
     if(any(is_map(:,icr)>0)) then
        do is=is_top,1,-1
           if(is_map(is,icr)>0) then
              is_map(is,icr)=-99999
              if(icr==color_rank) ns = ns-1
              exit crloop
           endif
        enddo
     endif
  enddo crloop
           
  ns_blk_tmp = count(is_map>0)
  call check(ns_blk, ns_blk_tmp, ' ns, ns_blk_tmp ')
 
  call map_states

contains
  subroutine calc_map_n(ig_map, ig_top, ng_blk, ng)
    implicit none
    integer ig_top, ng_blk, ng, ig  ! ng: output
    integer ig_map(ig_top,0:cs_pt-1) ! output
    integer icnt, icr, ng_i, ng_s
    
    ig_map = -999999 ! intentional

    ng      = ng_blk/cs_pt
    if(color_rank<mod(ng_blk, cs_pt)) ng=ng+1

    ng_s = ng
    call color_allsum_scalar_i(ng_s)
    call check(ng_s, ng_blk,' ng_s, ng_blk ')
    
    icnt = 0
    ig_map = 0
    do icr=0,cs_pt-1
       ng_i = ng_blk/cs_pt
       if(icr<mod(ng_blk, cs_pt)) ng_i = ng_i+1
       do ig=1, ng_i
          icnt = icnt+1
          ig_map(ig, icr) = icnt
       enddo
    enddo
    call check(icnt, ng_blk, ' icnt, ng_blk ')
  end subroutine calc_map_n

  subroutine map_states
    implicit none
    integer ir, i, i_strt, ns_sum, i0, ntop_in

    ! 
    ! rechecking
    !
    ns_sum = ns
    call color_allsum_scalar_i(ns_sum)
    call check(ns_sum,ns_blk,' ns_sum, ns_blk ')

    if(det_dft) then ! deterministic deterministic dft
       ntop_in = ceiling(rnel*dble(nsp)/2d0-1d-8)
       if(det_tddft) call check(ns_sum, ntop_in ,' ns_sum, ntop ')
    endif

    allocate(is_start(0:cs-1), stat=st); call check0(st,' is_start ')
    allocate(is_end(  0:cs-1), stat=st); call check0(st,' is_end   ')
    allocate(ns_each( 0:cs-1), stat=st); call check0(st,' ns_each  ')
 
    i0=0
    do ir=0,cs-1
       if(color_rank<ir) then
          i=ns
       else
          i=0
       endif
       call color_allsum_scalar_i(i) 

       if(i+1.le.ns_blk) then
          is_start(ir)= i+1
       else
          is_start(ir) = 0
       endif
    enddo

    do ir=0,cs-1
       if(color_rank.le.ir) then
          i=ns
       else
          i=0
       endif
       call color_allsum_scalar_i(i) 
       if(i.le.ns_blk.and.is_start(ir)/=0) then
          is_end(ir)= i
       else
          is_end(ir) = -1
       end if
    enddo

    do ir=0,cs-1
       if(is_start(ir)>is_end(ir).and.is_start(ir)>0)then
          write(6,*)' problem, rank,ir, is_start, is_end ',rank,ir,is_start(ir),is_end(ir)
          stop
       end if
    end do

    do ir=0,cs-1
       ns_each(ir) = max(0,is_end(ir)+1-is_start(ir))
    enddo 

    if(rank==0) then
       write(17,*)' is_start(0:color_size-1)= ',is_start
       write(17,*)'   is_end(0:color_size-1)= ',is_end
       write(17,*)'  ns_each(0:color_size-1)= ',ns_each
       call flush(17)
    end if

    call check(ns_blk, sum(ns_each),' ns_blk, sum_ns_each ')

    !write(60000+rank,*)' is_start ',is_start
    !write(60000+rank,*)' is_end   ',is_end
    !write(60000+rank,*)' ns_each  ',ns_each
    !write(60000+rank,*)' ns       ',ns
    !call flush(60000+rank)

  end subroutine map_states
end subroutine color_flgs
