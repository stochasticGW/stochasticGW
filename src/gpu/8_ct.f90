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
subroutine make_ct_gpu
        
  use gwm
  use simple_mpi, only : rank
  use device_mem_module

  implicit none
  type (cudaGraph) :: graph
  type (cudaGraphExec) :: graph_exec
  type (cudaGraphNode) :: error_node
  character(c_char) :: buffer
  integer(c_size_t) :: buffer_len
  integer :: it,ie,ij,st,j,ncolcrci
  real*8, allocatable :: crveta(:,:),civeta(:,:),crvxi(:,:),civxi(:,:),qri(:,:)

! Error handling
  if (gamflg.ne.2) stop ' prop_etaxi_gpu(): gamflg !=2 '
  if (.not.device_setup) stop ' make_ct_gpu(): device not setup'
  if (.not.ffts_setup) stop ' make_ct_gpu(): ffts not setup'
  if (.not.propdtpt_setup) stop ' make_ct_gpu(): propdtpt not setup'
  if (.not.vomake_setup) stop ' make_ct_gpu(): vomake not setup'
  if (.not.allocated(etaxi)) stop ' etaxi not allocated'
  if (.not.allocated(ge)) stop ' ge not allocated'
  if (.not.allocated(vo)) stop ' vo not allocated'
  if (.not.allocated(cveta)) stop ' cveta not allocated'
  if (.not.allocated(cvxi)) stop '  cvxi not allocated'
  if (.not.allocated(ct)) stop ' ct not allocated'
  if (SIZE(etaxi,1).ne.n .or. SIZE(etaxi,2).ne.2) then
     write(*,'(X,4(A,I0),A)') 'etaxi must be (',n,' x ',2,&
                              ') but is (',SIZE(etaxi,1),' x ',SIZE(etaxi,2),')'
     stop ' make_ct_gpu(): wrong size: etaxi'
  endif
  if (SIZE(ge,1).ne.n .or. SIZE(ge,2).ne.ne) then
     write(*,'(X,4(A,I0),A)') 'ge must be (',n,' x ',ne,&
                              ') but is (',SIZE(ge,1),' x ',SIZE(ge,2),')'
     stop ' make_ct_gpu(): wrong size: ge'
  endif
  if (SIZE(vo,1).ne.n .or. SIZE(vo,2).ne.nspv) then
     write(*,'(X,4(A,I0),A)') 'vo must be (',n,' x ',nspv,&
                              ') but is (',SIZE(vo,1),' x ',SIZE(vo,2),')'
     stop ' make_ct_gpu(): wrong size: vo'
  endif
  if (SIZE(cveta).ne.ne) then
     write(*,'(X,2(A,I0),A)') 'cveta must be (',ne,') but is (',SIZE(cveta),')'
     stop ' make_ct_gpu(): wrong size: cveta'
  endif
  if (SIZE(cvxi).ne.ne) then
     write(*,'(X,2(A,I0),A)') 'cvxi must be (',ne,') but is (',SIZE(cvxi),')'
     stop ' make_ct_gpu(): wrong size: cvxi'
  endif
  if (SIZE(ct,1).ne.(2*nt+1) .or. SIZE(ct,2).ne.ne) then
     write(*,'(X,4(A,I0),A)') 'ct must be (',2*nt+1,' x ',ne,&
                              ') but is (',SIZE(ct,1),' x ',SIZE(ct,2),')'
     stop ' make_ct_gpu(): wrong size: ct'
  endif

! nj-dependent code
  if (nj.gt.0) then
     if (SIZE(gej,1).ne.n .or. SIZE(gej,2).ne.ne .or. SIZE(gej,3).ne.nj) then
        write(*,'(X,6(A,I0),A)') 'gej must be (',n,' x ',ne,' x ',nj,&
                                 ') but is (',SIZE(gej,1),' x ',SIZE(gej,2),&
                                 ' x ',SIZE(gej,3),')'
        stop ' make_ct_gpu(): wrong size: gej'
     endif
     if (SIZE(cvetaj,1).ne.ne .or. SIZE(cvetaj,2).ne.nj) then
        write(*,'(X,4(A,I0),A)') 'cvetaj must be (',ne,' x ',nj,&
                              ') but is (',SIZE(cvetaj,1),' x ',SIZE(cvetaj,2),')'
        stop ' make_ct_gpu(): wrong size: cvetaj'
     endif
     if (SIZE(cvxij,1).ne.ne .or. SIZE(cvxij,2).ne.nj) then
        write(*,'(X,4(A,I0),A)') 'cvxij must be (',ne,' x ',nj,&
                              ') but is (',SIZE(cvxij,1),' x ',SIZE(cvxij,2),')'
        stop ' make_ct_gpu(): wrong size: cvxij'
     endif
     if (SIZE(ctj,1).ne.(2*nt+1) .or. SIZE(ctj,2).ne.ne .or. SIZE(ctj,3).ne.nj) then
        write(*,'(X,6(A,I0),A)') 'ctj must be (',2*nt+1,' x ',ne,' x ',nj,&
                                 ') but is (',SIZE(ctj,1),' x ',SIZE(ctj,2),&
                                 ' x ',SIZE(ctj,3),')'
        stop ' make_ct_gpu(): wrong size: ctj'
     endif
  endif

! Temporary arrays
  ncolcrci=max(nj,1)
  allocate(crveta(ne,ncolcrci), stat=st); &
           if(st/=0) stop ' make_ct_gpu(): trouble allocating crveta '
  allocate(civeta(ne,ncolcrci), stat=st); &
           if(st/=0) stop ' make_ct_gpu(): trouble allocating civeta '
  allocate(crvxi(ne,ncolcrci), stat=st); &
           if(st/=0) stop ' make_ct_gpu(): trouble allocating crvxi '
  allocate(civxi(ne,ncolcrci), stat=st); &
           if(st/=0) stop ' make_ct_gpu(): trouble allocating civxi '
  allocate(qri(n,2), stat=st); &
           if(st/=0) stop ' make_ct_gpu(): trouble allocating qri '

! Overwrite pt0 with etaxi before propagation
  pt0(:,1,1:2)=etaxi(:,1:2)
  !$acc update device(pt0(:,1,1:2)) !!! add 'ft' here if previously computed on host

  !$acc data copyin(ge,gej) &
  !$acc create(vo,cveta,cvxi,cvetaj,cvxij,crveta,civeta,crvxi,civxi,qri) &
  !$acc copyout(ct,ctj)

! vks never changes, so compute expvks before propagation
  call make_expvks

  !$acc data present(ct)
  !$acc parallel loop gang vector collapse(2) async
  do ie=1,ne
     do it=-nt,nt
        ct(it,ie)=0.d0
     enddo
  enddo
  !$acc end data

  if (nj.gt.0) then
     !$acc data present(ctj)
     !$acc parallel loop gang vector collapse(3) async
     do ij=1,nj
        do ie=1,ne
           do it=-nt,nt
              ctj(it,ie,ij)=0.d0
           enddo
        enddo
      enddo
      !$acc end data
  endif

! Create a cuda graph to minimize kernel latency. Note that vo1_make()
! and ct_make() have it-dependent pointers, so they cannot be captured
  if (makecudagraph) then
     st=cudaStreamBeginCapture(stream_g,0)
     call cv_prep
     if (nj.gt.0) call cvj_prep
     call propdt_pt(n,1,.TRUE.)
     st=cudaStreamEndCapture(stream_g,graph)
     st=cudaGraphInstantiate(graph_exec,graph,error_node,buffer,buffer_len)
  endif

! Propagate
  do it=0,nt
     if (block_gam_alg) then
        call vo1_make_block
     else
        call vo1_make_random
     endif
     if (makecudagraph) then
        st=cudaGraphLaunch(graph_exec,stream_g)
        st=cudaStreamSynchronize(stream_g)
     else
        call cv_prep
        if (nj.gt.0) call cvj_prep
        call propdt_pt(n,1,.TRUE.)
     endif
     call ct_make
     if (nj.gt.0) call ctj_make
  enddo
  !$acc wait
  !$acc end data

  deallocate(crveta,civeta,crvxi,civxi, stat=st); &
   if(st/=0) stop ' make_ct_gpu(): trouble deallocating crveta,civeta,crvxi,civxi '
  deallocate(qri, stat=st); &
   if(st/=0) stop ' make_ct_gpu(): trouble deallocating qri '

contains

  subroutine make_expvks

    implicit none
    integer :: i,j
    complex*16, parameter :: ci = (0d0,1d0)
    complex*16 :: fac

    fac=-ci*dt/2d0

!   Only one column of expvks needed for the selected spin direction
    !$acc data present(expvks,vks)
    !$acc parallel loop gang vector async 
    do i=1,n
       expvks(i,1,1)=exp(fac*vks(i,sp0))
    enddo
    !$acc end data

  end subroutine make_expvks

  subroutine vo1_make_block

    implicit none
    integer :: j,k,ssj
    integer :: i,js,jf,nseggrp,nsamples
    real*8  :: sf,ftr,fti,qqr,qqi
    complex*16, parameter :: ci = (0d0,1d0)

    nseggrp=(n+seg-1)/seg
    nsamples=ngam/nseggrp
    sf=seg_fctr/dble(ngam_nzero_blk)

!   This version assumes that segsh and segw are blocked, not randomly generated

    !$acc data present(vo,gam,ft,qri,seg_sh,seg_w)
    !$acc parallel loop gang vector_length(256) private(js,jf,k,qqr,qqi) async
    do i=1,n
       js=(i-1)/seg * nsamples + 1
       jf=js+nsamples-1
       k=MOD(i-1,seg)+1
       qqr=0.d0 ; qqi=0.d0
       !$acc loop seq
       do j=js,jf
          ftr=dble(ft(it,j)) ; fti=aimag(ft(it,j))
          qqr=qqr+gam(k,j)*ftr ; qqi=qqi+gam(k,j)*fti
       enddo
       vo(i,1) = (qqr+ci*qqi) * sf
    enddo
    !$acc end data

  end subroutine vo1_make_block

  subroutine vo1_make_random

    implicit none
    integer :: j,k,ssj
    real*8  :: sf,ftr,fti,qqr,qqi
    complex*16, parameter :: ci = (0d0,1d0)

    sf=seg_fctr/dble(ngam_nzero_blk)

    !$acc data present(vo,gam,ft,qri,seg_sh,seg_w)
    !$acc parallel loop gang vector collapse(2) async
    do j=1,2
       do k=1,n
          qri(k,j)=0.d0
       enddo
    enddo

!   Caution: segments can overlap, so qri accesses must be atomic
!   No complex*16 atomic support yet in CUDA, so do two real reductions
    !$acc parallel loop gang vector_length(256) private(ftr,fti,ssj) async
    do j=1,ngam
       ftr=dble(ft(it,j))
       fti=aimag(ft(it,j))
       ssj=seg_sh(j)
       do k=1,seg_w(j)
          !$acc atomic
          qri(k+ssj,1) = qri(k+ssj,1) + gam(k,j)*ftr          
          !$acc atomic
          qri(k+ssj,2) = qri(k+ssj,2) + gam(k,j)*fti
       enddo
    enddo

    !$acc parallel loop gang vector async
    do k=1,n
       vo(k,1) = (qri(k,1)+ci*qri(k,2)) * sf
    enddo
    !$acc end data

  end subroutine vo1_make_random

  subroutine cv_prep

    implicit none
    integer :: i,iout,ie
    complex*16 :: tmpeta,tmpxi,fac
    complex*16, parameter   :: ci=(0d0,1d0)

    ! Array etaxi(:,1) is in pt0(:,1,1); etaxi(:,2) in pt0(:,1,2)
    !$acc data present(cvxi,cveta,pt0,vo,ge,crveta,civeta,crvxi,civxi)
    !$acc parallel loop gang vector async
    do ie=1,ne
       crveta(ie,1) = 0.d0; civeta(ie,1)=0.d0 ; crvxi(ie,1)=0.d0 ; civxi(ie,1)=0.d0
    enddo

    !$acc parallel loop private(tmpeta,tmpxi,fac) vector_length(256) collapse(2) async
    do ie=1,ne
        do iout=1,n,256
           tmpeta=(0.d0,0.d0)
           tmpxi=(0.d0,0.d0)
           !$acc loop vector reduction(+:tmpeta) reduction(+:tmpxi) shortloop
           do i=iout,min(iout+255,n)
              fac=vo(i,1)*ge(i,ie)
              tmpeta=tmpeta+fac*conjg(pt0(i,1,1))
              tmpxi =tmpxi +fac*      pt0(i,1,2)
           enddo
           !$acc atomic
           crveta(ie,1)=crveta(ie,1)+dble(tmpeta)
           !$acc atomic
           civeta(ie,1)=civeta(ie,1)+aimag(tmpeta)
           !$acc atomic
           crvxi(ie,1)=crvxi(ie,1)+dble(tmpxi)
           !$acc atomic
           civxi(ie,1)=civxi(ie,1)+aimag(tmpxi)
        enddo
    enddo

    !$acc parallel loop gang vector async
    do ie=1,ne
       cveta(ie)=(crveta(ie,1)+ci*civeta(ie,1))*dv ; &
       cvxi(ie)=(crvxi(ie,1)+ci*civxi(ie,1))*dv
    enddo

    !$acc end data

  end subroutine cv_prep

  subroutine cvj_prep

    implicit none
    integer :: i,iout,ie,ij
    complex*16 :: tmpeta,tmpxi,fac
    complex*16, parameter   :: ci=(0d0,1d0)

    ! Array etaxi(:,1) is in pt0(:,1,1); etaxi(:,2) in pt0(:,1,2)
    !$acc data present(cvxij,cvetaj,pt0,vo,ge,crveta,civeta,crvxi,civxi)
    !$acc parallel loop gang vector collapse(2) async
    do ij=1,nj
       do ie=1,ne
          crveta(ie,ij) = 0.d0; civeta(ie,ij)=0.d0 ; crvxi(ie,ij)=0.d0 ; civxi(ie,ij)=0.d0
       enddo
    enddo

    !$acc parallel loop private(tmpeta,tmpxi,fac) vector_length(256) collapse(3) async
    do ij=1,nj
       do ie=1,ne
          do iout=1,n,256
             tmpeta=(0.d0,0.d0)
             tmpxi=(0.d0,0.d0)
             !$acc loop vector reduction(+:tmpeta) reduction(+:tmpxi) shortloop
             do i=iout,min(iout+255,n)
                fac=vo(i,1)*gej(i,ie,ij)
                tmpeta=tmpeta+fac*conjg(pt0(i,1,1))
                tmpxi =tmpxi +fac*      pt0(i,1,2)
             enddo
             !$acc atomic
             crveta(ie,ij)=crveta(ie,ij)+dble(tmpeta)
             !$acc atomic
             civeta(ie,ij)=civeta(ie,ij)+aimag(tmpeta)
             !$acc atomic
             crvxi(ie,ij)=crvxi(ie,ij)+dble(tmpxi)
             !$acc atomic
             civxi(ie,ij)=civxi(ie,ij)+aimag(tmpxi)
          enddo
       enddo
    enddo

    !$acc parallel loop gang vector collapse(2) async
    do ij=1,nj
       do ie=1,ne
          cvetaj(ie,ij)=(crveta(ie,ij)+ci*civeta(ie,ij))*dv ; &
          cvxij(ie,ij)=(crvxi(ie,ij)+ci*civxi(ie,ij))*dv
       enddo
    enddo

    !$acc end data

  end subroutine cvj_prep

  subroutine ct_make

    implicit none
    integer :: ie
    real*8  :: wg

    wg = wgf(it); if (it==0) wg=wg/2d0
    !$acc data present(cvxi,cveta,ct)
    !$acc parallel loop gang vector async
    do ie=1,ne
       ct( it,ie) = ct( it,ie) + cvxi(ie)*normxi*wg ; ct(-it,ie) = ct(-it,ie) - cveta(ie)*normeta*wg
    enddo
    !$acc end data

  end subroutine ct_make

  subroutine ctj_make

    implicit none
    integer :: ie,ij
    real*8  :: wg

    wg = wgf(it); if (it==0) wg=wg/2d0
    !$acc data present(cvxij,cvetaj,ctj)
    !$acc parallel loop gang vector collapse(2) async
    do ij=1,nj
       do ie=1,ne
          ctj( it,ie,ij) = ctj( it,ie,ij) + cvxij(ie,ij)*normxi*wg ; &
          ctj(-it,ie,ij) = ctj(-it,ie,ij) - cvetaj(ie,ij)*normeta*wg
       enddo
    enddo
    !$acc end data

  end subroutine ctj_make

end subroutine make_ct_gpu
