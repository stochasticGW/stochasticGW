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
!
! NOTE: need to check sum over g, wherther g is defined and what to divide by
!

module gwpltm
  use simple_mpi, only : allsum_r8, allsum_scalar_r8, allsum_c16, rank
  implicit none
  save
  real*8, allocatable     :: sge(:), sr2(:,:), si2(:,:), dra(:,:), dia(:,:)
  real*8, allocatable     :: swg(:)
  real*8, allocatable     :: sx(:), sx2(:), dex(:)
  real*8, allocatable     :: sv(:), sv2(:), dev(:)
  real*8, allocatable     :: fx(:),  fv(:)
  complex*16, allocatable :: cw(:),  scw(:,:)
  real*8, allocatable     :: sgej(:,:), sr2j(:,:,:), si2j(:,:,:), draj(:,:,:), diaj(:,:,:)!wf
  real*8, allocatable     :: swgj(:,:), svj(:,:),    sv2j(:,:),   devj(:,:)!wf
  real*8, allocatable     :: sxrj(:,:), sx2rj(:,:),  dexrj(:,:)!wf
  real*8, allocatable     :: sxij(:,:), sx2ij(:,:),  dexij(:,:)!wf
  real*8, allocatable     :: fxrj(:,:), fxij(:,:), fvj(:,:)!wf
  complex*16, allocatable :: cwj(:,:),  scwj(:,:,:)!wf 

end module gwpltm

subroutine gwplt
  use gwpltm
  use gwm,        only : ne, imc, ct,  nw, dw, nt, dt, work_dir
  use gwm,        only : output_dir, ear, war, print_now
  use gwm,        only : imctot, stoch_x, det_tddft
  use gwm,        only : ctj, nj, orbj_indx !wf
!  use gwm,        only : xlrnew, functional !!! PT pre-bnl merge
  use simple_mpi, only : color_size
  implicit none
  integer, save :: i1=1
  integer ie,st, iw, iwmn, iwmx, i, nds
  real*8     :: sigr_ks, sigi_ks, dsigr_ks, dsigi_ks, z, elin
  real*8     :: sgm, eq, eiq, deq, deiq, dsgm !wf
  integer ij!wf
  real*8     :: sigr_ksj(nj), sigi_ksj(nj), dsigr_ksj(nj), dsigi_ksj(nj)!wf
  real*8     :: sgmrj(nj), sgmij(nj), dsgmrj(nj), dsgmij(nj)!wf

  call prnt(1,'GW Plotting results')

  call prep_plot
  call prep_sge
  call accum_ex_vks
  call make_sigw
  call dev_sigw
  call plot_sigw
  call plot_qp
  !call make_plot_sigt
contains
  subroutine prep_plot
    use gwm, only : nmc
    use simple_mpi, only : nodes
    implicit none
    
    nds = max(nodes,1)
    imctot = imc*(nds/color_size)
    
    if(mod(nw,2)/=0) stop ' ERROR: nw should be even '
    iwmn = -nw/2
    iwmx =  nw/2-1
    
    firstime : if(i1==1) then
       allocate( sge(ne), stat=i); if(i/=0) stop ' sge problems '
       sge = 0d0

       if(nj>0) then!wf
         allocate( sgej(ne,nj), stat=i); if(i/=0) stop ' sgej problems '!wf
         sgej = 0d0!wf
       endif!wf

       if(rank==0) then
          allocate(scw(iwmn:iwmx,ne),stat=i); if(i/=0) stop ' scw problems '
          allocate(sr2(iwmn:iwmx,ne),stat=i); if(i/=0) stop ' sr2 problems '
          allocate(si2(iwmn:iwmx,ne),stat=i); if(i/=0) stop ' si2 problems '
          allocate(dra(iwmn:iwmx,ne),stat=i); if(i/=0) stop ' dra problems '
          allocate(dia(iwmn:iwmx,ne),stat=i); if(i/=0) stop ' dia problems '
          allocate(swg(ne),          stat=i); if(i/=0) stop ' swg problems '
          allocate(sx(ne),           stat=i); if(i/=0) stop ' sx  problems '
          allocate(sx2(ne),          stat=i); if(i/=0) stop ' sx2 problems '
          allocate(dex(ne),          stat=i); if(i/=0) stop ' dex problems '
          allocate(fx(ne),           stat=i); if(i/=0) stop ' fx  problems '
          allocate(sv(ne),           stat=i); if(i/=0) stop ' sv  problems '
          allocate(sv2(ne),          stat=i); if(i/=0) stop ' sv2 problems '
          allocate(dev(ne),          stat=i); if(i/=0) stop ' dev problems '
          allocate(fv(ne),           stat=i); if(i/=0) stop ' fv  problems '

          scw = 0d0
          sr2 = 0d0
          si2 = 0d0
          swg = 0d0
          sx  = 0d0
          sx2 = 0d0
          sv  = 0d0
          sv2 = 0d0

          if(nj>0) then!wf
             allocate(scwj(iwmn:iwmx,ne,nj),stat=i); if(i/=0) stop ' scwj problems '!wf
             allocate(sr2j(iwmn:iwmx,ne,nj),stat=i); if(i/=0) stop ' sr2j problems '!wf
             allocate(si2j(iwmn:iwmx,ne,nj),stat=i); if(i/=0) stop ' si2j problems '!wf
             allocate(draj(iwmn:iwmx,ne,nj),stat=i); if(i/=0) stop ' draj problems '!wf
             allocate(diaj(iwmn:iwmx,ne,nj),stat=i); if(i/=0) stop ' diaj problems '!wf
             allocate(swgj(ne,nj),          stat=i); if(i/=0) stop ' swgj problems '!wf
             allocate(svj(ne,nj),           stat=i); if(i/=0) stop ' sv j problems '!wf
             allocate(sv2j(ne,nj),          stat=i); if(i/=0) stop ' sv2j problems '!wf
             allocate(devj(ne,nj),          stat=i); if(i/=0) stop ' devj problems '!wf
             allocate(sxrj(ne,nj),          stat=i); if(i/=0) stop ' sxrj problems '!wf
             allocate(sx2rj(ne,nj),         stat=i); if(i/=0) stop 'sx2rj problems '!wf
             allocate(dexrj(ne,nj),         stat=i); if(i/=0) stop 'dexrj problems '!wf
             allocate(sxij(ne,nj),          stat=i); if(i/=0) stop ' sxij problems '!wf
             allocate(sx2ij(ne,nj),         stat=i); if(i/=0) stop 'sx2ij problems '!wf
             allocate(dexij(ne,nj),         stat=i); if(i/=0) stop 'dexij problems '!wf
             allocate(fxrj(ne,nj),          stat=i); if(i/=0) stop ' fxrj problems '!wf
             allocate(fxij(ne,nj),          stat=i); if(i/=0) stop ' fxij problems '!wf
             allocate(fvj (ne,nj),          stat=i); if(i/=0) stop ' fvj  problems '!wf
   
             scwj = 0d0!wf
             sr2j = 0d0!wf
             si2j = 0d0!wf
             swgj = 0d0!wf
             svj  = 0d0!wf
             sv2j = 0d0!wf
             sxrj = 0d0!wf
             sx2rj= 0d0!wf
             sxij = 0d0!wf
             sx2ij= 0d0!wf
          endif!wf
  
       end if
       i1=-1
    end if firstime

  end subroutine prep_plot
  
  subroutine prep_sge
    use gwm, only : ge, g, dv, ptheader_flg
    use gwm, only : gej !wf
    implicit none
    real*8 geg
    real*8 gegj(nj)
    do ie=1,ne
       if(ptheader_flg) then
          geg= sum(ge(:,ie)*g(:))*dv
          
          if(nj>0) then!wf
             do ij=1, nj!wf
               gegj(ij)= sum(gej(:,ie,ij)**2)*dv!wf
             enddo!wf
          endif!wf
       
       else
          geg = 0d0
          if (nj>0) gegj=0d0 !wf
       end if
       call allsum_scalar_r8(geg)
       if(nj>0) call allsum_r8(gegj,size(gegj)) !wf
       sge(ie) = sge(ie) + geg
       if(nj>0) sgej(ie,:) = sgej(ie,:) + gegj !wf
    enddo
  end subroutine prep_sge

  subroutine accum_ex_vks
    use simple_mpi, only : rank, nodes
    use gwm, only        : wge, exce, vxce, ptheader_flg, trace, wgej, nj, excej, vxcej !wf
    implicit none
    integer :: ie
    real*8, allocatable :: aw(:),ax(:), av(:), ax2(:), av2(:)
    
    real*8, allocatable :: awj(:,:), avj(:,:), av2j(:,:)!wf
    real*8, allocatable :: axrj(:,:), ax2rj(:,:)!wf
    real*8, allocatable :: axij(:,:), ax2ij(:,:)!wf
    
    allocate(aw(ne),ax(ne),av(ne),ax2(ne),av2(ne),stat=st); call check0(st,' ax-w ')
    if(nj/=0) allocate(awj(ne,nj),axrj(ne,nj),axij(ne,nj),avj(ne,nj),&!wf
        ax2rj(ne,nj),ax2ij(ne,nj),av2j(ne,nj),stat=st); call check0(st,' ax-wj ') !wf
   
    if(ptheader_flg) then
       aw = wge 
       ax = exce
       av = vxce
!       if(functional=='bnl') av = av + wge*xlrnew !!! PT pre-bnl merge
       ax2 = (ax/wge)**2*wge
       av2 = (av/wge)**2*wge

       if(nj/=0) then!wf
          awj = wgej !wf
          axrj= dble(excej)!wf
          axij= aimag(excej)!wf
          avj = vxcej!wf
        
          do ie=1,ne!wf
            do ij=1,nj!wf
              av2j(ie,ij)  = &!wf
                  vxcej(ie,ij)**2/(sqrt(wge(ie))*sqrt(wgej(ie,ij)))!wf
              ax2rj(ie,ij) = &!wf
                  dble(excej(ie,ij))**2/(sqrt(wge(ie))*sqrt(wgej(ie,ij)))!wf
              ax2ij(ie,ij) = &!wf
                  aimag(excej(ie,ij))**2/(sqrt(wge(ie))*sqrt(wgej(ie,ij)))!wf
            enddo!wf
          enddo!wf
       endif!wf
      
    else
       aw = 0d0
       ax = 0d0
       av = 0d0
       ax2 = 0d0
       av2 = 0d0

       if(nj/=0) then!wf
          awj = 0d0!wf
          axrj = 0d0!wf
          axij = 0d0!wf
          avj = 0d0!wf
          av2j  = 0d0!wf
          ax2rj = 0d0!wf
          ax2ij = 0d0!wf
       endif  !wf
       
    end if

    call allsum_r8(aw, size(aw))
    call allsum_r8(ax, size(ax))
    call allsum_r8(av, size(av))
    call allsum_r8(ax2,size(ax2))
    call allsum_r8(av2,size(av2))
    
    if(nj/=0) then!wf
       call allsum_r8(awj, size(awj))!wf
       call allsum_r8(avj, size(avj))!wf
       call allsum_r8(av2j, size(av2j))!wf
       call allsum_r8(axrj, size(axrj))!wf
       call allsum_r8(axij, size(axij))!wf
       call allsum_r8(ax2rj,size(ax2rj))!wf
       call allsum_r8(ax2ij,size(ax2ij))!wf
    endif!wf

    rnk0 : if(rank==0) then
       swg = swg + aw
       sv  = sv  + av
       sv2 = sv2 + av2
       
       if(nj/=0) then!wf
         swgj = swgj + awj!wf
         svj  = svj  + avj!wf
         sv2j = sv2j + av2j!wf
       endif!wf

       if(stoch_x) then
          sx  = sx  + ax
          sx2 = sx2 + ax2
          
          if(nj/=0) then!wf
             sxrj  = sxrj  + axrj!wf
             sxij  = sxij  + axij!wf
             sx2rj = sx2rj + ax2rj!wf
             sx2ij = sx2ij + ax2ij!wf
          endif!wf

       endif

       neloop : do ie=1,ne          

          fv(ie) = sv(ie)/swg(ie)
          call stat_err(sv2(ie),sv(ie),swg(ie),imctot,dev(ie),' dev ')
  
          if(nj/=0) then!wf
             do ij=1,nj!wf
                fvj(ie,ij) = svj(ie,ij)/sqrt(swg(ie)*swgj(ie,ij))!wf
                call stat_err(sv2j(ie,ij),svj(ie,ij),sqrt(swgj(ie,ij)*swg(ie)),imctot,devj(ie,ij),' devj ') !wf
             enddo!wf
          endif!wf  

          if(stoch_x)  then;  
             fx(ie) = sx(ie)/swg(ie)
             call stat_err(sx2(ie),sx(ie),swg(ie),imctot,dex(ie),' dex ')
     
             if(nj/=0) then!wf
                do ij=1,nj!wf
                   fxrj(ie,ij)=sxrj(ie,ij)/sqrt(swg(ie)*swgj(ie,ij))!wf
                   fxij(ie,ij)=sxij(ie,ij)/sqrt(swg(ie)*swgj(ie,ij))!wf
                   call stat_err(sx2rj(ie,ij),sxrj(ie,ij),sqrt(swgj(ie,ij)*swg(ie)),imctot,dexrj(ie,ij),' dexrj ') !wf
                   call stat_err(sx2ij(ie,ij),sxij(ie,ij),sqrt(swgj(ie,ij)*swg(ie)),imctot,dexij(ie,ij),' dexij ') !wf
                enddo!wf
             endif!wf
     
          else;               
             fx(ie) = exce(ie)
             dex(ie) = 0d0 
     
             if(nj/=0) then!wf
                do ij=1,nj!wf
                   fxrj(ie,ij)=dble(excej(ie,ij))!wf
                   fxij(ie,ij)=aimag(excej(ie,ij))!wf
                   dexrj(ie,ij)=0d0!wf
                   dexij(ie,ij)=0d0!wf
                enddo!wf
             endif!wf
     
          endif          
       enddo neloop

    end if rnk0
    deallocate(aw,ax,av,ax2,av2)
  end subroutine accum_ex_vks

  subroutine make_sigw
    use gwm, only : g, ge, dv, gej, nj!wf
    use simple_mpi, only : color_rank
    implicit none
    real*8                  :: geg
    real*8,     allocatable :: br2(:), bi2(:)
    complex*16, allocatable :: cw(:), ct1(:)
    
    real*8                  :: gegj(nj)!wf
    real*8,     allocatable :: br2j(:,:), bi2j(:,:)!wf
    complex*16, allocatable :: cwj(:,:)!wf

    allocate( cw( iwmn:iwmx),stat=i); if(i/=0) stop '  cw problems '
    allocate( br2(iwmn:iwmx),stat=i); if(i/=0) stop ' br2 problems '
    allocate( bi2(iwmn:iwmx),stat=i); if(i/=0) stop ' bi2 problems '
    allocate( ct1(size(ct,1)), stat=i); if(i/=0) stop ' ERROR: ct1 problems '

    if(nj/=0) then!wf
       allocate( cwj( iwmn:iwmx,nj),stat=i); if(i/=0) stop '  cwj problems '!wf
       allocate( br2j(iwmn:iwmx,nj),stat=i); if(i/=0) stop ' br2j problems '!wf
       allocate( bi2j(iwmn:iwmx,nj),stat=i); if(i/=0) stop ' bi2j problems '!wf
    endif!wf
    
    !
    ! checking
    !
    neloop : do ie=1,ne
       if(color_rank==0) then
          ct1=ct(:,ie)
          call ct_to_cw(ct1, cw, nt, nw, dt) 
          geg= sum(ge(:,ie)*g(:))*dv

          br2 = (dble( cw/geg))**2*geg
          bi2 = (aimag(cw/geg))**2*geg
  
          if(nj/=0) then!wf
             do ij=1,nj!wf
                ct1=ctj(:,ie,ij)!wf
                call ct_to_cw(ct1,cwj(:,ij),nt,nw,dt)!wf
                gegj(ij) = sum(gej(:,ie,ij)**2)*dv!wf
              
                br2j(:,ij)=dble( cwj(:,ij))**2/sqrt(geg)/sqrt(gegj(ij))!wf
                bi2j(:,ij)=aimag(cwj(:,ij))**2/sqrt(geg)/sqrt(gegj(ij))!wf
             enddo!wf
         endif!wf
          
       else
          cw=0d0
          geg = 0d0
          br2 = 0d0
          bi2 = 0d0
  
          if(nj/=0) then!wf
             cwj=0d0!wf
             gegj=0d0!wf
             br2j=0d0!wf
             bi2j=0d0!wf
          endif!wf
  
       end if

       call allsum_c16(cw, size(cw))
       call allsum_r8( br2, size(br2))
       call allsum_r8( bi2, size(bi2))
       
       if(nj/=0) then!wf
          call allsum_c16(cwj, size(cwj))!wf
          call allsum_r8( br2j, size(br2j))!wf
          call allsum_r8( bi2j, size(bi2j))!wf
       endif!wf
       
       if(rank==0) then
          scw(:,ie)  = scw(:,ie)  + cw(:)
          sr2(:,ie)  = sr2(:,ie)  + br2(:) 
          si2(:,ie)  = si2(:,ie)  + bi2(:) 
          
          if(nj/=0) then!wf
             scwj(:,ie,:)  = scwj(:,ie,:)  + cwj(:,:)!wf
             sr2j(:,ie,:)  = sr2j(:,ie,:)  + br2j(:,:) !wf
             si2j(:,ie,:)  = si2j(:,ie,:)  + bi2j(:,:) !wf
          endif!wf

       end if
    enddo neloop

    deallocate( cw)
    deallocate( br2)
    deallocate( bi2)
    deallocate( ct1) 
    
    if(nj/=0) then!wf
      deallocate( cwj)!wf
      deallocate( br2j)!wf
      deallocate( bi2j)!wf
    endif!wf

  end subroutine make_sigw

  subroutine dev_sigw
    implicit none
    integer ie, iw
    if(rank==0) then
       do ie=1,ne
          do iw=iwmn,iwmx
             call stat_err(sr2(iw,ie), dble(scw(iw,ie)),sge(ie),imctot,dra(iw,ie),'dra')
             call stat_err(si2(iw,ie),aimag(scw(iw,ie)),sge(ie),imctot,dia(iw,ie),'dia')

             if(nj/=0) then!wf
                do ij=1,nj!wf
                   call stat_err(sr2j(iw,ie,ij), dble(scwj(iw,ie,ij)),sqrt(sge(ie))*sqrt(sgej(ie,ij)),imctot,draj(iw,ie,ij),'draj')!wf
                   call stat_err(si2j(iw,ie,ij),aimag(scwj(iw,ie,ij)),sqrt(sge(ie))*sqrt(sgej(ie,ij)),imctot,diaj(iw,ie,ij),'diaj')!wf
                enddo!wf
             endif!wf
     
          end do
       end do
    end if
  end subroutine dev_sigw


  subroutine plot_sigw
    implicit none
    character(len=1024) :: filename
    character(len=1024) :: fm,fe
    integer ip, ie
    if(rank/=0) return
    if(.not.print_now)return

    eloop : do ie=1,ne
       
       if(nj==0) call find_linear(ie,sigr_ks, sigi_ks, dsigr_ks, dsigi_ks, z, elin)!wf
       if(nj==0) call find_qp(ie, z, sgm, eq, deq, eiq, deiq)  !wf
       if(nj/=0) call find_linear_j(ie,sigr_ks, sigi_ks, dsigr_ks, dsigi_ks, z, elin,&!wf
                      sigr_ksj, sigi_ksj, dsigr_ksj, dsigi_ksj)!wf
       if(nj/=0) call find_qp_j(ie, z, sgm, dsgm, eq, deq, eiq, deiq,&!wf
                      sgmrj, sgmij, dsgmrj, dsgmij)  !wf

       select case(ie)
       case(0  :9 );       write(fe,"('',I1)")ie
       case(10 :99);       write(fe,"('',I2)")ie
       case default; write(6,*)' problem with ie printing ',ie; stop
       end select

       ip=18
       if(ne>1) then
          open(ip,file=trim(adjustl(output_dir))//'/SigW.'//trim(adjustl(fe)),status='replace')
       else
          open(ip,file=trim(adjustl(output_dir))//'/SigW.txt',status='replace')
       endif
       call plot_sigw_to(ip,ie)
       close(ip)

       ip=18
       fm=' '
       select case(imctot)
       case(0  :9 );         write(fm,"('000000',I1)")imctot
       case(10 :99);         write(fm,"('00000',I2)")imctot
       case(100:999);        write(fm,"('0000',I3)")imctot
       case(1000:9999);      write(fm,"('000',I4)")imctot
       case(10000:99999);    write(fm,"('00',I5)")imctot
       case(100000:999999);  write(fm,"('0',I6)")imctot
       case(1000000:9999999);write(fm,"('',I7)")imctot
       case default; write(6,*)' problem with imctot ',imctot; stop
       end select 
       
       filename = &
       trim(adjustl(work_dir))//'/SigW.'//trim(adjustl(fe))//'.'//trim(adjustl(fm))

       open(ip,file=trim(adjustl(filename)),status='replace')
       call plot_sigw_to(ip,ie)
       close(ip)
    end do eloop

  end subroutine plot_sigw

  subroutine plot_sigw_to(ip,ie)
    implicit none
    integer ip,ie,iw
    write(ip,*)'# '   
    write(ip,*)'# SigW output from stochastic GW '   
    write(ip,*)'# ',iwmx-iwmn+1,' plotting frequencies '
    write(ip,*)'# '   
    write(ip,*)'# MC_STEP (Monte Carlo steps on each core) '
    write(ip,*)'# ',imc
    write(ip,*)'# '
    write(ip,*)'# N_cores '
    write(ip,*)'# ',nds
    write(ip,*)'# '
    write(ip,*)'# Buffer_size (# cores per indep. run) '
    write(ip,*)'# ',color_size
    write(ip,*)'# '
    write(ip,*)'# MC (accumulative # Monte Carlo runs) =MC_STEP*(N_cores/buffer_size)'
    write(ip,*)'# ',imctot
    write(ip,*)'# '
    write(ip,*)'# KS energy '
    write(ip,*)'# ',real(ear(ie))
    write(ip,*)'# '
    write(ip,*)'#  <X>              d<X> '
    write(ip,*)'# ',real(fx(ie)), real(dex(ie))
    write(ip,*)'# '
    write(ip,*)'# <vxc>                  '
    write(ip,*)'# ',real(fv(ie))
    write(ip,*)'#'
    write(ip,*)'# e_qp_real,       de_qp_real  '
    write(ip,*)'# ',real(eq), real(deq)
    write(ip,*)'#'
    write(ip,*)'# e_qp_imag,       de_qp_imag  '
    write(ip,*)'# ',real(eiq), real(deiq)
    write(ip,*)'#'
    
    if(nj/=0) then!wf
      do ij=1,nj!wf
        write(ip,*)'# j=', ij!wf
        write(ip,*)'# <X>ji_R           d<X>ji_R  '!wf
        write(ip,*)'# ',real(fxrj(ie,ij)),real(dexrj(ie,ij))!wf
        write(ip,*)'# '!wf
        write(ip,*)'# <X>ji_I           d<X>ji_I  '!wf
        write(ip,*)'# ',real(fxij(ie,ij)),real(dexij(ie,ij))!wf
        write(ip,*)'# '!wf
        write(ip,*)'# <vxc>ji                     '!wf
        write(ip,*)'# ',real(fvj(ie,ij))!wf
        write(ip,*)'#'!wf
      enddo!wf
    endif!wf
      
    write(ip,*)'#     w             SigW_real       SigW_imag        dSigW_real      dSigw_imag '
    write(ip,*)'# '
    do iw=iwmn,iwmx
       write(ip,"(' ',5e16.6)")iw*dw, scw(iw,ie)/sge(ie), dra(iw,ie), dia(iw,ie)
    enddo

    if(nj/=0) then!wf
      do ij=1,nj!wf
        write(ip,*) '# off diagonal with j=',ij!wf
        do iw=iwmn,iwmx!wf
          write(ip,"(' ',5e16.6)")iw*dw, scwj(iw,ie,ij)/sqrt(sge(ie)*sgej(ie,ij)), draj(iw,ie,ij), diaj(iw,ie,ij)!wf
        enddo!wf
      enddo!wf
    endif!wf
    
  end subroutine plot_sigw_to

  subroutine stat_err(sa2,sa1,g,imctot,da,ach)
    implicit none
    integer imctot
    real*8 sa2,sa1,g,da, d2
    character(*) ach
    
    d2 = sa2/g-(sa1/g)**2

    if(d2>0d0) then
       da = sqrt(d2/dble(imctot))
    else
       if(abs(d2)<1d-10) then
          da = 0d0
       else
          write(6,*)ach,' is negative, stopping ',d2
          stop
       endif
    end if
  end subroutine stat_err

  subroutine plot_qp
    use simple_mpi, only : nodes
    implicit none
    character(len=10) :: se='+-stat.err'
    rnk0prnt: if(rank==0.and.print_now) then

       eloop : do ie=1,ne
          
          if(nj==0) call find_linear(ie,sigr_ks, sigi_ks, dsigr_ks, dsigi_ks, z, elin)!wf
          if(nj==0) call find_qp(ie, z, sgm, eq, deq, eiq, deiq)  !wf
          if(nj/=0) call find_linear_j(ie,sigr_ks, sigi_ks, dsigr_ks, dsigi_ks, z, elin, &!wf
             sigr_ksj, sigi_ksj, dsigr_ksj, dsigi_ksj)!wf
          if(nj/=0) call find_qp_j(ie, z, sgm, dsgm, eq, deq, eiq, deiq,&!wf
             sgmrj, sgmij, dsgmrj, dsgmij)  !wf
  
          write(6,*)
          write(6,1) ear(ie),                      'E                 '
          write(6,1) fv(ie),                       '<V>               '
          write(6,2) imctot,fx(ie),dex(ie),        'MC, <X>           ',se
          
          if(nj==0) write(6,2)imctot,sgm, deq,     'MC, Sig_R(e_qp)   ',se!wf
          if(nj/=0) write(6,2)imctot,sgm,dsgm,     'MC, Sig_R(e_qp)   ',se!wf

          write(6,2) imctot,sigr_ks, dsigr_ks,     'MC, Sig_R(e_ks)   ',se
          write(6,2) imctot,sigi_ks, dsigi_ks,     'MC, Sig_I(e_ks)   ',se
          write(6,3) imctot,Z,                     'MC, Z             '
          write(6,3) imctot,Elin,                  'MC, Elin          '
          write(6,2) imctot,eq,deq,                'MC, E_qp          ',se
          write(6,2) imctot,eiq,deiq,              'MC, Ei_qp         ',se
          
          if(nj/=0) then!wf
!             write(6,*)!wf
             do ij=1,nj!wf
                write(6,*)!wf
                write(6,4) fvj(ie,ij),                            '<vxc>ji; j =',orbj_indx(ij)!wf
                write(6,2) imctot,fxrj(ie,ij),dexrj(ie,ij),       'MC, <X>ji real    ',se!wf
                write(6,2) imctot,fxij(ie,ij),dexij(ie,ij),       'MC, <X>ji imag    ',se!wf
                write(6,2) imctot,sigr_ksj(ij),dsigr_ksj(ij),     'MC, Sig_R(e_ks)ji ',se!wf
                write(6,2) imctot,sigi_ksj(ij),dsigi_ksj(ij),     'MC, Sig_I(e_ks)ji ',se!wf
                write(6,2) imctot,sgmrj(ij),dsgmrj(ij),           'MC, Sig_R(e_qpi)ji',se!wf
                write(6,2) imctot,sgmij(ij),dsgmij(ij),           'MC, Sig_I(e_qpi)ji',se!wf     
             enddo!wf
          endif!wf
 
1         format(x,13x,1f15.6,15x,10x,A,X)
2         format(x,i13,2f15.6,10x,A,X,A)
3        format(x,i13,1f15.6,15x,10x,A)
4        format(x,13x,1f15.6,15x,10x,A,X,I0)

       enddo eloop
       write(6,*)
       call flush(6)
    end if rnk0prnt
  end subroutine plot_qp

  subroutine find_linear(ie,sigr_ks, sigi_ks, dsigr_ks, dsigi_ks, z, elin)
    implicit none
    integer    ie, iw
    real*8     sigr_ks, sigi_ks, dsigr_ks, dsigi_ks, z, elin
    real*8     r0, r1
    complex*16 sg, sq
    iw=floor(ear(ie)/dw);  
    if(iw>iwmx-2.or.iw<iwmn+3) stop ' w problems '
    if(war(iw).gt.ear(ie).or.war(iw+1).lt.ear(ie)) stop ' we problems '
    r0 = (war(iw+1)-ear(ie))/dw
    r1 = ( ear(ie) -war(iw))/dw     
    if(abs(r0+r1-1d0)>1d-4.or.min(r0,r1)<-0.001.or.max(r0,r1)>1.001) &
         stop ' r0,r problems, stopping '

    sg =  scw(iw,  ie)*r0+ scw(iw+1,ie )*r1
    sq = (scw(iw+1,ie)   - scw(iw  ,ie))/dw 

    sg = sg/sge(ie)
    sq = sq/sge(ie)

    sigr_ks= dble(sg)
    sigi_ks= aimag(sg)

    dsigr_ks = dra(iw, ie)
    dsigi_ks = dia(iw, ie)
    z = 1d0/(1d0- dble(sq))  ! xxx check if not wrong, need to get it from dividing by weight

    !
    ! E(w=E)=  eks + Sig(w=E) + X-vks
    ! Assume
    !       Sig(w)= Sig(eks)+ dS/dw * (w-eks)
    ! So E = eks + Sig(eks) + dS/dw * (E-eks) +X-Vks
    ! i.e., (E-eks)*(1-dS/dw) = (Sig(eks) + X-Vks)
    ! i.e., E = eks + Z( Sig(eks)+X-Vks)

    elin = ear(ie) + Z*(sigr_ks+fx(ie)-fv(ie))

  end subroutine find_linear
    
  subroutine find_qp(ie, z, sgm, eq, deq, eiq, deiq)  
    implicit none
    integer iw, ie, iw0, iwu, iwd
    real*8 sgm, eq, deq, eiq, deiq, r0, r1, z
    real*8, allocatable :: f(:)
    logical zerou, zerod
    allocate(f(iwmn:iwmx), stat=st); call check0(st,' f in qp ')

    f = -war + ear(ie) + scw(:,ie)/sge(ie) + fx(ie)-fv(ie) ! (sx(ie)-sv(ie))/swg(ie)

    iw=floor(ear(ie)/dw);  
    if(iw>iwmx-2.or.iw<iwmn+3) stop ' w problems '
    if(war(iw).gt.ear(ie).or.war(iw+1).lt.ear(ie)) stop ' we problems '

    iw0 = iw
    ! first go up 
    zerou = .false.
    do iw=iw0,iwmx-2
       if(f(iw)*f(iw+1)<0d0) then
          iwu = iw
          zerou = .true.
       endif
    enddo
    zerod = .false.
    do iw=iw0-1,iwmn+2,-1
       if(f(iw)*f(iw+1)<0d0) then
          iwd = iw
          zerod = .true.
       endif
    enddo

    if((.not.zerou).and.(.not.zerod)) stop ' no fit found '
    if(.not.zerou) iw = iwd
    if(.not.zerod) iw = iwu
    if(zerou.and.zerod) then
       iw = iwu
       if(abs(iwd-iw0)<abs(iwu-iw0)) iw =  iwd
    endif

    if(f(iw)*f(iw+1).ge.0d0) stop ' problem in fit '
    
    r1 = (0d0-f(iw))/(f(iw+1)-f(iw))
    r0 = 1d0-r1
    if(r1<0d0.or.r1>1d0) stop ' r1 problems '

    eq   = r0*war(iw)           +r1*war(iw+1)
    sgm  =  dble(r0*scw(iw,ie)+r1*scw(iw+1,ie))/sge(ie)
    eiq  = aimag(r0*scw(iw,ie)+r1*scw(iw+1,ie))/sge(ie)
    deq  = r0*dra(iw,ie)        +r1*dra(iw+1,ie)
    deiq = r0*dia(iw,ie)        +r1*dia(iw+1,ie)

    deq  = sqrt(deq**2 + z**2*dex(ie)**2)
    
    deallocate(f)
    
  end subroutine find_qp

!wf it's the entire subroutine
  subroutine find_linear_j(ie,sigr_ks, sigi_ks, dsigr_ks, dsigi_ks, z, elin, &
             sigr_ksj, sigi_ksj, dsigr_ksj, dsigi_ksj)
    implicit none
    integer    ie, iw, ij
    real*8     sigr_ks, sigi_ks, dsigr_ks, dsigi_ks, z, elin
    real*8     sigr_ksj(nj), sigi_ksj(nj), dsigr_ksj(nj), dsigi_ksj(nj)
    real*8     r0, r1
    complex*16 sg, sq
    iw=floor(ear(ie)/dw);  
    if(iw>iwmx-2.or.iw<iwmn+3) stop ' w problems '
    if(war(iw).gt.ear(ie).or.war(iw+1).lt.ear(ie)) stop ' we problems '
    r0 = (war(iw+1)-ear(ie))/dw
    r1 = ( ear(ie) -war(iw))/dw     
    if(abs(r0+r1-1d0)>1d-4.or.min(r0,r1)<-0.001.or.max(r0,r1)>1.001) &
         stop ' r0,r problems, stopping '

    sg =  scw(iw,  ie)*r0+ scw(iw+1,ie )*r1
    sq = (scw(iw+1,ie)   - scw(iw  ,ie))/dw 

    sg = sg/sge(ie)
    sq = sq/sge(ie)

    sigr_ks= dble(sg)
    sigi_ks= aimag(sg)

    dsigr_ks = dra(iw, ie)
    dsigi_ks = dia(iw, ie)
    z = 1d0/(1d0- dble(sq))  ! xxx check if not wrong, need to get it from dividing by weight

    !
    ! E(w=E)=  eks + Sig(w=E) + X-vks
    ! Assume
    !       Sig(w)= Sig(eks)+ dS/dw * (w-eks)
    ! So E = eks + Sig(eks) + dS/dw * (E-eks) +X-Vks
    ! i.e., (E-eks)*(1-dS/dw) = (Sig(eks) + X-Vks)
    ! i.e., E = eks + Z( Sig(eks)+X-Vks)

    elin = ear(ie) + Z*(sigr_ks+fx(ie)-fv(ie))

    do ij=1,nj
      sg = scwj(iw, ie, ij)*r0+ scwj(iw+1, ie, ij)*r1
      sg = sg/sqrt(sge(ie)*sgej(ie,ij))

      sigr_ksj(ij) = dble(sg)
      sigi_ksj(ij) = aimag(sg)

      dsigr_ksj(ij) = draj(iw,ie,ij)
      dsigi_ksj(ij) = diaj(iw,ie,ij)
    enddo
  
  end subroutine find_linear_j

  subroutine find_qp_j(ie, z, sgm, dsgm, eq, deq, eiq, deiq, &
             sgmrj, sgmij, dsgmrj, dsgmij)  
    implicit none
    integer iw, ie, iw0, iwu, iwd, ij
    real*8 sgm, dsgm, eq, deq, eiq, deiq, r0, r1, z
    real*8 sgmrj(nj), sgmij(nj), dsgmrj(nj), dsgmij(nj)
    real*8, allocatable :: f(:)
    logical zerou, zerod
    allocate(f(iwmn:iwmx), stat=st); call check0(st,' f in qp ')
    f = -war + ear(ie) + scw(:,ie)/sge(ie) + fx(ie)-fv(ie) ! (sx(ie)-sv(ie))/swg(ie)

    iw=floor(ear(ie)/dw);  
    if(iw>iwmx-2.or.iw<iwmn+3) stop ' w problems '
    if(war(iw).gt.ear(ie).or.war(iw+1).lt.ear(ie)) stop ' we problems '

    iw0 = iw
    ! first go up 
    zerou = .false.
    do iw=iw0,iwmx-2
       if(f(iw)*f(iw+1)<0d0) then
          iwu = iw
          zerou = .true.
       endif
    enddo
    zerod = .false.
    do iw=iw0-1,iwmn+2,-1
       if(f(iw)*f(iw+1)<0d0) then
          iwd = iw
          zerod = .true.
       endif
    enddo

    if((.not.zerou).and.(.not.zerod)) stop ' no fit found '
    if(.not.zerou) iw = iwd
    if(.not.zerod) iw = iwu
    if(zerou.and.zerod) then
       iw = iwu
       if(abs(iwd-iw0)<abs(iwu-iw0)) iw =  iwd
    endif

    if(f(iw)*f(iw+1).ge.0d0) stop ' problem in fit '
    
    r1 = (0d0-f(iw))/(f(iw+1)-f(iw))
    r0 = 1d0-r1
    if(r1<0d0.or.r1>1d0) stop ' r1 problems '

    eq   = r0*war(iw)           +r1*war(iw+1)
    sgm  =  dble(r0*scw(iw,ie)+r1*scw(iw+1,ie))/sge(ie)
    eiq  = aimag(r0*scw(iw,ie)+r1*scw(iw+1,ie))/sge(ie)
    deq  = r0*dra(iw,ie)        +r1*dra(iw+1,ie)
    deiq = r0*dia(iw,ie)        +r1*dia(iw+1,ie)

    deq  = sqrt(deq**2 + z**2*dex(ie)**2)
    
    do ij=1,nj
      sgmrj(ij) = dble(r0*scwj(iw,ie,ij)+r1*scwj(iw+1,ie,ij))/sqrt(sge(ie)*sgej(ie,ij))
      sgmij(ij) =aimag(r0*scwj(iw,ie,ij)+r1*scwj(iw+1,ie,ij))/sqrt(sge(ie)*sgej(ie,ij))
      dsgmrj = r0*draj(iw,ie,ij)        +r1*draj(iw+1,ie,ij)
      dsgmij = r0*diaj(iw,ie,ij)        +r1*diaj(iw+1,ie,ij)
    enddo
    
    deallocate(f)
    
  end subroutine find_qp_j

end subroutine gwplt
