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

subroutine vo_make_gpu
  use gwm
  use simple_mpi, only : rank
  use device_mem_module

  implicit none
  integer :: it, st
  complex*16, parameter :: ci=(0.d0,1.d0)
  real*8, allocatable :: norm(:),spnfacs(:), vht0(:,:), rho_pp(:,:), gamvrd(:,:,:)
  complex*16, allocatable :: crho(:,:)
  complex*8, allocatable :: fw(:,:) ! *8 or *16

! Error handling
  if (gamflg/=2.or.flgdyn) then
     write(*,*) 'For GPU calculation, must have flgdyn=.FALSE., gamflg=2'
     stop ' vo_make_gpu(): bad flgdyn, gamflg '
  endif
  if (ngam<1) stop ' vo_make_gpu(): ngam < 1 '
  if (.not.device_setup) stop ' vo_make_gpu(): device not setup'
  if (.not.ffts_setup) stop ' vo_make_gpu(): ffts not setup'
  if (.not.propdtpt_setup) call init_propdtpt_device 
  if (.not.vomake_setup) call init_vomake_device
  if (.not.allocated(pt)) stop ' vo_make_gpu(): pt not allocated'
  if (.not.allocated(del)) stop ' vo_make_gpu(): del not allocated'
  if (.not.allocated(vk)) stop ' vo_make_gpu(): vk not allocated'
  if (SIZE(pt,1).ne.n .or. SIZE(pt,2).ne.ns) then
     write(*,'(X,4(A,I0),A)') 'pt must be (',n,' x ',ns,&
                              ') but is (',SIZE(pt,1),' x ',SIZE(pt,2),')'
     stop ' vo_make_gpu(): wrong size: pt'
  endif
  if (SIZE(del).ne.n) then
     write(*,'(X,2(A,I0),A)') 'del must be (',n,') but is (',SIZE(del),')'
     stop ' vo_make_gpu(): wrong size: del'
  endif
  if (SIZE(vk).ne.(n*(scale_vh)**3)) then
     write(*,'(X,2(A,I0),A)') 'vk must be (',n*(scale_vh)**3,') but is (',SIZE(vk),')'
     stop ' vo_make_gpu(): wrong size: vk'
  endif
  allocate(norm(2),spnfacs(ns), stat=st); &
           if(st/=0) stop ' vo_make_gpu(): trouble allocating norm,spinfacs '
  allocate(vht0(n,2),rho_pp(n,2),crho(n*(scale_vh)**3,2),gamvrd(ngam,nsp,2), stat=st); &
           if(st/=0) stop ' vo_make_gpu(): trouble allocating vht0,rho_pp,crho,gamvrd '
  allocate(fw(0:nw-1,cnk), stat=st); &
           if(st/=0) stop ' vo_make_gpu(): trouble allocating fw '
  call makespinfacs
! Put global 'pt' array into local-to-module 'pt0'; the latter has the added 
! dimension for running the perturbed and unperturbed propagations in parallel
  pt0(:,:,1)=pt(:,:)
  !$acc update device(pt0(:,:,1))
  !$acc data copyin(spnfacs,del,vk)  &
  !$acc create(norm,vht0,rho_pp,crho,gamvrd,fw)

  call vr_make_gpu      ! propagation
  call vr_to_vo_gam_gpu ! time-ordering

  !$acc end data
  ! acc update host(ft) !!! Add '$' before 'acc' to copy 'ft' to host after run

  deallocate(norm, stat=st);
             if(st/=0) stop ' vo_make_gpu(): trouble deallocating norm '
  deallocate(spnfacs, stat=st);
             if(st/=0) stop ' vo_make_gpu(): trouble deallocating spnfacs '
  deallocate(vht0,rho_pp,crho,gamvrd, stat=st); &
             if(st/=0) stop ' vo_make_gpu(): trouble deallocating vht0,rho_pp,crho,gamvrd '
  deallocate(fw, stat=st);
             if(st/=0) stop ' vo_make_gpu(): trouble deallocating fw '

contains

  subroutine makespinfacs

    implicit none
    integer :: is

!   Compute the spin factors
    do is=1,ns
       if (det_tddft) then
          spnfacs(is) = occ_det(is) ! no 1/ns, no 2.
       else
          spnfacs(is) = 2d0/nsp*normpt(is)**2d0/ns_blk ! spin: 2.  
       end if
    enddo

  end subroutine makespinfacs

  subroutine vr_make_gpu

    implicit none
    type (cudaGraph) :: graph
    type (cudaGraphExec) :: graph_exec
    type (cudaGraphNode) :: error_node
    character(c_char) :: buffer
    integer(c_size_t) :: buffer_len
    integer :: i,j,k
    complex*16, parameter :: ci=(0.d0,1.d0)
    complex*16 :: fac

    fac=-ci*sm

!   Perturb the eta fxns (eq. 14)
    !$acc data present(pt0,del,qr,qi)
    !$acc parallel loop gang vector collapse(2) async
    do j=1,ns
       do i=1,n
          pt0(i,j,2)=exp(fac*del(i) )*pt0(i,j,1)
       enddo
    enddo

!   Set qr,qi to zero before propagating to ensure that the correct 
!   subset of the grid is processed by propnl() 
    !$acc parallel loop gang vector collapse(3) async
    do k=1,2
       do j=1,ns
          do i=1,n
             qr(i,j,k)=0.d0 ; qi(i,j,k)=0.d0
          enddo
       enddo
    enddo
    !$acc end data

!   Create a cuda graph. Note that 'it' is set to 1 to avoid running the 
!   kernel to fill vht0 in makevt(); also, the pointers in fr1_gama()
!   depend on 'it' so this kernel cannot be captured.
    if (makecudagraph) then
       it=1
       st=cudaStreamBeginCapture(stream_g,0)
       call makerho
       call makevt
       call propdt_pt(n,ns,.FALSE.)
       call gamvrmake
       st=cudaStreamEndCapture(stream_g,graph)
       st=cudaGraphInstantiate(graph_exec,graph,error_node,buffer,buffer_len)
    endif

!   Propagate both unperturbed, perturbed etas (post-eq. 14)
    do it=0,nt
      if (it.eq.0 .or. (.not.makecudagraph)) then
         call makerho
         call makevt
         call propdt_pt(n,ns,.FALSE.)
         call gamvrmake
      else
         st=cudaGraphLaunch(graph_exec,stream_g)
         st=cudaStreamSynchronize(stream_g)
      endif
      call fr1_gama
    enddo
    !$acc wait

  end subroutine vr_make_gpu

  subroutine makerho
    implicit none
    integer :: i,j,is,iout
    real*8  :: ntmp,val

    val=rnel/dv

!   Sum the spin-weighted densities for each spatial point in eta
    !$acc data present(rho_pp,norm)
    !$acc data present(pt0,spnfacs)
    !$acc parallel loop gang collapse(2) async
    do j=1,2
       do i=1,n
          rho_pp(i,j)=0.d0
          !$acc loop seq
          do is=1,ns
             rho_pp(i,j) = &
             rho_pp(i,j)+spnfacs(is)*((dble(pt0(i,is,j)))**2+(aimag(pt0(i,is,j)))**2 )
          end do
       enddo
    enddo
    !$acc end data

!   Normalize density for each spatial point at this time step (eq. 15)
    !$acc parallel loop gang vector async
    do j=1,2
       norm(j)=0.d0
    enddo

    !$acc parallel loop gang private(ntmp) collapse(2) vector_length(256) async
    do j=1,2
       do iout=1,n,256
          ntmp=0.d0
          !$acc loop vector reduction(+:ntmp) shortloop
          do i=iout,min(iout+255,n)
             ntmp=ntmp+rho_pp(i,j)
          enddo
          !$acc atomic
          norm(j)=norm(j)+ntmp
       enddo
    enddo

    !$acc parallel loop gang vector collapse(2) async
    do j=1,2
       do i=1,n
          rho_pp(i,j)=rho_pp(i,j) * (val/norm(j))
       enddo
    enddo
    !$acc end data

  end subroutine makerho

  subroutine makevt
    use gwm
    implicit none
    integer :: i,j,sp
    complex*16, parameter :: ci=(0.d0,1.d0)
    complex*16 :: fac

    fac=-ci*dt/2.d0

!   Apply 1/r operator to density (pre-eq. 1)
    call rho_p_to_vht_deviceonly

    !$acc data present(vht0,rho_pp)
    if (it.eq.0) then
       !$acc parallel loop gang vector collapse(2) async
       do j=1,2
          do i=1,n
             vht0(i,j)=rho_pp(i,j)
          enddo
       enddo
    endif

!   Difference used here and in gamvrmake()
!   'rho_pp' contains v_H^(lam)(r,t) after this kernel
    !$acc parallel loop gang vector collapse(2) async
    do j=1,2
       do i=1,n
          rho_pp(i,j)=rho_pp(i,j)-vht0(i,j)
       enddo
    enddo

    !$acc data present(expvks,vks)
    !$acc parallel loop gang vector collapse(3) async
    do j=1,2
       do sp=1,nsp
          do i=1,n
             expvks(i,sp,j)=exp(fac*(vks(i,sp)+rho_pp(i,j)))
          enddo
       enddo
    enddo
    !$acc end data
    !$acc end data

  end subroutine makevt
  
  subroutine rho_p_to_vht_deviceonly
    use gwm, only : nx, ny, nz, n, scale_vh
    use openacc
    use cudafor
    use cufft

    implicit none
    integer    :: i,j,rho_os,crho_os,nbig
    integer    :: nxy,nxyb,nxb,nyb,nzb,mx,my,mz,ix,iy,iz
    integer(4) :: cerr

    nxb = scale_vh * nx
    mx = nxb/2 - nx/2
    nyb = scale_vh * ny
    my = nyb/2 - ny/2
    nzb = scale_vh * nz
    mz = nzb/2 - nz/2

    nxy=nx*ny
    nxyb=nxb*nyb
    nbig=n*(scale_vh)**3

    !$acc data present(crho,rho_pp)
    if (scale_vh.eq.1) then
       !$acc parallel loop gang vector collapse(2) async
       do j=1,2
          do i=1,n
             crho(i,j) = rho_pp(i,j)
          enddo
       enddo
    else
       !$acc parallel loop gang vector collapse(2) async
       do j=1,2
          do i=1,nbig
             crho(i,j) = 0.d0
          enddo
       enddo

       !$acc parallel loop gang private(rho_os,crho_os) collapse(3) async
       do j=1,2
          do iz=1,nz
             do iy=1,ny
                rho_os = (iz-1)*nxy + (iy-1)*nx
                crho_os = (iz-1+mz)*nxyb + (iy-1+my)*nxb + mx
                !$acc loop vector
                do ix=1,nx
                   crho(crho_os+ix,j)=rho_pp(rho_os+ix,j)
                enddo
             enddo
          enddo
       enddo
    endif

    !$acc host_data use_device(crho)
    cerr = cufftExecZ2Z(cufft_plan_rho, crho, crho, CUFFT_FORWARD)
    !$acc end host_data
    if (cerr/=0) stop 'rho_p_to_vht_deviceonly(): problem cutfft_forward'
  
    !$acc data present(vk)
    !$acc parallel loop gang vector collapse(2) async
    do j=1,2
       do i=1,nbig
          crho(i,j) = crho(i,j)*vk(i)
       enddo
    enddo
    !$acc end data

    !$acc host_data use_device(crho)
    cerr = cufftExecZ2Z(cufft_plan_rho, crho, crho, CUFFT_INVERSE)
    !$acc end host_data
    if (cerr/=0) stop 'rho_p_to_vht_deviceonly(): problem cutfft_inverse'

    if (scale_vh.eq.1) then
       !$acc parallel loop gang vector collapse(2) async
       do j=1,2
          do i=1,n
             rho_pp(i,j) = dble(crho(i,j))
          enddo
       enddo
    else
       !$acc parallel loop gang private(rho_os,crho_os) collapse(3) async
       do j=1,2
          do iz=1,nz
             do iy=1,ny
                rho_os = (iz-1)*nxy + (iy-1)*nx
                crho_os = (iz-1+mz)*nxyb + (iy-1+my)*nxb + mx
                !$acc loop vector
                do ix=1,nx
                   rho_pp(rho_os+ix,j)=dble(crho(crho_os+ix,j))
                enddo
             enddo
          enddo
       enddo
    endif
    !$acc end data

  end subroutine rho_p_to_vht_deviceonly

  subroutine gamvrmake
    use gwm
    implicit none
    integer :: i,j,k,iout
    real*8  :: dot

!   Compute overlaps <xi_j | u_xi(t)>, i.e. compute
!   v_xi^(lam)(r,t) in 'gamvr' from v_H^(lam)(r,t) in 'rho_pp'
    !$acc data present(rho_pp,seg_sh,seg_w,gam,gamvrd)

!    ! acc parallel loop gang private(dot) vector_length(32) collapse(2) async
!    do i=1,2
!       do j=1,ngam
!          dot=0d0
!          ! acc loop reduction(+:dot)
!          ! acc loop seq
!          do k=1,seg_w(j)
!            dot=dot+gam(k,j)*rho_pp(k+seg_sh(j),i)
!          enddo
!          gamvrd(j,1,i)=dot*dv
!          if (nsp>1) gamvrd(j,nsp,i) = gamvrd(j,1,i)
!       enddo
!    enddo
!    ! acc end parallel


    !$acc parallel loop gang vector collapse(2) async
    do i=1,2
       do j=1,ngam
          gamvrd(j,1,i)=0.d0
       enddo
    enddo

    !$acc parallel loop gang private(dot) vector_length(32) collapse(3) async
    do i=1,2
       do j=1,ngam
          do iout=1,seg,32
             dot=0d0
             !$acc loop seq
             do k=iout,min(iout+31,seg_w(j))
                dot=dot+gam(k,j)*rho_pp(k+seg_sh(j),i)
             enddo
             !$acc atomic
             gamvrd(j,1,i)=gamvrd(j,1,i)+dot*dv
          enddo
       enddo
    enddo

    if (nsp>1) then
       !$acc parallel loop gang vector collapse(2) async
       do i=1,2
          do j=1,ngam
             gamvrd(j,nsp,i) = gamvrd(j,1,i)
          enddo
       enddo
    endif

    !$acc end data

  end subroutine gamvrmake
  
  subroutine fr1_gama
    implicit none
    integer :: j
    real*8  :: wg

!   NOTE: the CPU version differs from below by storing the unperturbed
!   (real*8) gamvr(:,sp0) in (complex*8) ft(it,:), which loses half of
!   precision in the real part. To reproduce the CPU calculation, use:
!    ft(it,j) = wg*(dble(gamvrd(j,sp0,2)) - dble(sngl(gamvrd(j,sp0,1))))

!   Accumulate u^R(r,t)
    wg = wgf(it)/sm
    !$acc data present(gamvrd,ft)
    !$acc parallel loop gang vector async
    do j=1,ngam
!       ft(it,j) = wg*(dble(gamvrd(j,sp0,2)) - dble(gamvrd(j,sp0,1)))
       ft(it,j) = wg*(dble(gamvrd(j,sp0,2)) - dble(sngl(gamvrd(j,sp0,1))))
    enddo
    !$acc end data

  end subroutine fr1_gama

  subroutine vr_to_vo_gam_gpu
    use openacc
    use cudafor
    use cufft

    implicit none
    integer :: i,j,k,ih,clim
    integer(4) :: cerr

    !$acc data present(ft,fw)
    do i=1,ngam,cnk

       ih=min(i-1+cnk,ngam)
       clim=min(cnk,ngam-i+1)

!      Skip dividing fw(0,:) elements by 2 since they are guaranteed to be zero 
       !$acc parallel loop gang vector collapse(2) async
       do k=1,cnk
          do j=0,nw-1
             fw(j,k)=0.d0
          enddo
       enddo
       !$acc parallel loop gang vector collapse(2) async
       do k=1,clim
          do j=0,nt
             fw(j,k)=ft(j,i+k-1)/dble(nw)
          enddo
       enddo

       !$acc host_data use_device(fw)
       cerr = cufftExecC2C(cufft_plan_ft, fw, fw, CUFFT_INVERSE) ! C2C or Z2Z
       !$acc end host_data
       if (cerr/=0) stop 'vr_to_vo_gam_gpu(): problem cutfft_forward'

       !$acc parallel loop gang vector collapse(2) async
       do k=1,cnk
          do j=nw/2,nw-1
             fw(j,k)=conjg(fw(j,k))
          enddo
       enddo

       !$acc host_data use_device(fw)
       cerr = cufftExecC2C(cufft_plan_ft, fw, fw, CUFFT_FORWARD) ! C2C or Z2Z
       !$acc end host_data
       if (cerr/=0) stop 'vr_to_vo_gam_gpu(): problem cutfft_inverse'

       !$acc parallel loop gang vector collapse(2) async
       do k=1,clim
          do j=0,nt
             ft(j,i+k-1)=fw(j,k)
          enddo
       enddo

    enddo
    !$acc wait
    !$acc end data

  end subroutine vr_to_vo_gam_gpu

end subroutine vo_make_gpu
