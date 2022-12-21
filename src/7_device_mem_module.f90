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
module device_mem_module
  use simple_mpi, only : rank
  use iso_c_binding
  use openacc
  use cudafor
  use cufft

  implicit none
  save
  integer(4) :: cufft_plan_rho, cufft_plan(2), cufft_plan_ft, cufft_plan_pteta
  integer(kind=cuda_stream_kind) :: stream_g
  integer :: nalj
  integer, parameter :: cnk = 200
  integer, allocatable    :: ialj(:), jss(:), nbtt(:), map_sp_pteta(:)
  real*8,  allocatable    :: qr(:,:,:), qi(:,:,:), nrmarray(:,:)
  real*8,  allocatable    :: qcr(:,:), qci(:,:), shscvks(:,:), scpsi(:)
  real*8  :: two_ov_dh
  complex*16, allocatable :: pe(:,:), po(:,:), pa(:,:), ptmp(:,:)
  complex*16, allocatable :: pt0(:,:,:), expvks(:,:,:), ovdev(:)
  logical :: device_setup=.FALSE.
  logical :: ffts_setup=.FALSE.
  logical :: filter_setup=.FALSE.
  logical :: propdtpt_setup=.FALSE.
  logical :: vomake_setup=.FALSE.
  logical :: makecudagraph=.TRUE.

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine precompute_ovdev(chebfacs)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Allocates arrays for device propdt_pt module

      use gwm,    only : na, dv, dt, dh
      use kb_mod, only : lpptop, lpploc, mapai, iub, ngs, vp
      use kb_mod, only : vp, pvp
      logical, intent(in)     :: chebfacs
      complex*16, allocatable :: ovdevtmp(:)
      integer, allocatable    :: iatmp(:), jsstmp(:), nbttmp(:)
      integer :: ma, igg, k, l, j, ib, it, ia, nbt, jb, js
      real*8, parameter  :: toll_o = 1d-8
      real*8  :: dot
      complex*16, parameter   :: ci=(0d0,1d0)
      complex*16              :: fac
      
      fac=-ci*dt/2d0
      two_ov_dh=2d0/dh

      ALLOCATE(ovdevtmp(16*na),iatmp(16*na),jsstmp(16*na),nbttmp(16*na))

      k=1
      do ia=1,na
         ma=mapai(ia)
         ! determine starting point "u" for different l
         jb=iub(ia)
         do l=0,lpptop(ma)
            if (l/=lpploc(ma))then
               ib=l**2     ! 0,1,4,9
               it=ib+2*l   ! 0,3,8,15
               nbt=it-ib+1 ! 1,3,5,7

               ! ov=<vphi|vphi> in 3d.
               do j=ib,it
                  js=jb+(j-ib)

                  if (chebfacs) then
                     ovdevtmp(k) = dv * two_ov_dh / pvp(l,ma)
                  else
                     dot=0d0
                     do igg=1,ngs(ia)
                        dot= dot + vp(js+(igg-1)*nbt)**2
                     enddo
                     ovdevtmp(k) = max(dot*dv,toll_o)
                     ovdevtmp(k) = dv / ovdevtmp(k) * &
                     (exp(fac*ovdevtmp(k) / pvp(l,ma))-1d0)
                  endif

                  iatmp(k)=ia
                  jsstmp(k)=js
                  nbttmp(k)=nbt
                  k=k+1
              enddo
              jb=jb+nbt*ngs(ia)
            endif
         enddo
      enddo

!     Resize atom-l-j arrays
      nalj=k-1
      ALLOCATE(ialj(nalj),jss(nalj),nbtt(nalj))
      ialj(1:nalj)=iatmp(1:nalj)
      jss(1:nalj)=jsstmp(1:nalj)
      nbtt(1:nalj)=nbttmp(1:nalj)
      if (chebfacs) then
         allocate(scpsi(nalj))
         scpsi(1:nalj)=REAL(ovdevtmp(1:nalj))
      else
         allocate(ovdev(nalj))
         ovdev(1:nalj)=ovdevtmp(1:nalj)
      endif

      DEALLOCATE(ovdevtmp,iatmp,jsstmp,nbttmp)

      end subroutine precompute_ovdev

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine init_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets up device (GPU)

      implicit none
      integer :: thegpu,ndev

      ndev=acc_get_num_devices(acc_device_nvidia)
      thegpu=MOD(rank,ndev)
      call acc_set_device_num(thegpu,acc_device_nvidia)
!      write(*,'(X,3(A,I0),A)') 'MPI rank ',rank,' on device (',thegpu+1,'/',ndev,')'

!     cuFFT plans
      call init_ffts_device

      device_setup = .TRUE.

      end subroutine init_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine flush_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets up device (GPU)

      implicit none

!     cuFFT plans
      call flush_ffts_device

      device_setup = .FALSE.

      end subroutine flush_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine init_filter_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Allocates arrays for device Chebyshev filter module

      use gwm,    only : n, ns, nsp, sp0, dh, havg
      use gwm,    only : vks, ekn, map_sp_pt
      use kb_mod, only : mapkbg, ngs, vp

      implicit none
      integer :: st

!     Error handling
      if (.not.allocated(vks)) stop ' init_filter_device(): vks not allocated'
      if (.not.allocated(ekn)) stop ' init_filter_device(): ekn not allocated'

      if (SIZE(vks,1).ne.n .or. SIZE(vks,2).ne.nsp) then
         write(*,'(X,4(A,I0),A)') 'vks must be (',n,' x ',nsp,&
                                  ') but is (',SIZE(vks,1),' x ',SIZE(vks,2),')'
         stop ' init_filter_device(): wrong size: vks'
      endif
      if (SIZE(ekn).ne.n) then
         write(*,'(X,2(A,I0),A)') 'ekn must be (',n,') but is (',SIZE(ekn),')'
         stop ' init_filter_device(): wrong size: ekn'
      endif
      if (SIZE(map_sp_pt).ne.ns) then
         write(*,'(X,2(A,I0),A)') 'map_sp_pt must be (',ns,') but is (',SIZE(map_sp_pt),')'
         stop ' init_filter_device(): wrong size: map_sp_pt'
      endif

      allocate(shscvks(n,nsp), stat=st); &
               if(st/=0) stop ' init_filter_device(): trouble allocating shscvks '
      allocate(pe(n,ns+1),po(n,ns+1),ptmp(n,ns+1), stat=st); &
               if(st/=0) stop ' init_filter_device(): trouble allocating pa '
      allocate(qcr(n,ns+1),qci(n,ns+1), stat=st); if(st/=0) stop ' allocate  qr, qi '
      allocate(map_sp_pteta(ns+1), stat=st); &
               if(st/=0) stop ' init_filter_device(): trouble allocating map_sp_pteta '

      shscvks(:,:)=2d0/dh*(vks(:,:)-havg)
      map_sp_pteta(1:ns)=map_sp_pt(1:ns) ! spin mapping (pt)
      map_sp_pteta(ns+1)=sp0             ! spin mapping (eta)

      call precompute_ovdev(.TRUE.)

      !$acc enter data create(qcr,qci,pe,po,ptmp) &
      !$acc copyin(shscvks,ekn,map_sp_pteta) &
      !$acc copyin(mapkbg,ngs,vp,scpsi,ialj,jss,nbtt)

      filter_setup=.TRUE.

      end subroutine init_filter_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine flush_filter_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Dellocates arrays for device propdt_pt module

      use gwm,    only : ekn
      use kb_mod, only : mapkbg, ngs, vp

      implicit none
      integer :: st

      !$acc exit data delete(pe,po,ptmp,qcr,qci) & 
      !$acc delete(shscvks,ekn,map_sp_pteta) &
      !$acc delete(mapkbg,ngs,vp,scpsi,ialj,jss,nbtt)

      deallocate(shscvks, stat=st); if (st/=0) stop ' deallocate shscvks '
      deallocate(pe,po,ptmp, stat=st); if (st/=0) stop ' deallocate pe,po,ptmp '
      deallocate(qcr,qci, stat=st); if (st/=0) stop ' deallocate qr, qi '
      deallocate(map_sp_pteta, stat=st); if(st/=0) stop ' deallocate map_sp_pteta '
      deallocate(scpsi, stat=st); if (st/=0) stop ' deallocate scpsi '
      deallocate(ialj, stat=st); if (st/=0) stop ' deallocate ialj '
      deallocate(jss, stat=st); if (st/=0) stop ' deallocate jss '
      deallocate(nbtt, stat=st); if (st/=0) stop ' deallocate nbtt '

      filter_setup=.FALSE.
      
      end subroutine flush_filter_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine init_propdtpt_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Allocates arrays for device propdt_pt module

      use gwm,    only : n, ns, nsp, na
      use gwm,    only : expnk, map_sp_pt
      use kb_mod, only : mapkbg, ngs, vp

      implicit none
      integer :: st
!      integer :: thegpu,ndev

!      ndev=acc_get_num_devices(acc_device_nvidia)
!      thegpu=MOD(rank,ndev)
!      call acc_set_device_num(thegpu,acc_device_nvidia)
!      write(*,'(X,3(A,I0),A)') 'MPI rank ',rank,' on device (',thegpu+1,'/',ndev,')'

      allocate(pt0(n,ns,2), stat=st); if(st/=0) stop ' allocate pt0 '
      allocate(expvks(n,nsp,2),stat=st); if(st/=0) stop ' allocate expvks '
      allocate(qr(n,ns,2),qi(n,ns,2), stat=st); if(st/=0) stop ' allocate  qr, qi '
      allocate(nrmarray(ns,2), stat=st); if(st/=0) stop ' allocate nrmarray '

      call precompute_ovdev(.FALSE.)

!     Check for dimension mismatches
      if (SIZE(expnk).ne.n) then
         write(*,'(X,2(A,I0),A)') 'expnk must be (',n,') but is (',SIZE(expnk),')'
         stop ' init_propdtpt_device(): wrong size: expnk'
      endif
      if (SIZE(map_sp_pt).ne.ns) then
         write(*,'(X,2(A,I0),A)') 'map_sp_pt must be (',ns,') but is (',SIZE(map_sp_pt),')'
         stop ' init_propdtpt_device(): wrong size: map_sp_pt'
      endif

      !$acc enter data create(pt0,expvks,qr,qi,nrmarray) &
      !$acc copyin(expnk,map_sp_pt) &
      !$acc copyin(mapkbg,ngs,vp,ovdev,ialj,jss,nbtt)

      propdtpt_setup=.TRUE.

      end subroutine init_propdtpt_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine flush_propdtpt_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Dellocates arrays for device propdt_pt module

      use gwm,    only : expnk, map_sp_pt
      use kb_mod, only : mapkbg, ngs, vp

      implicit none
      integer :: st

      !$acc exit data delete(pt0,expvks,qr,qi,nrmarray) & 
      !$acc delete(expnk,map_sp_pt) &
      !$acc delete(mapkbg,ngs,vp,ovdev,ialj,jss,nbtt)

      deallocate(pt0, stat=st); if (st/=0) stop ' deallocate pt0 '
      deallocate(expvks, stat=st); if (st/=0) stop ' deallocate expvks '
      deallocate(qr,qi, stat=st); if (st/=0) stop ' deallocate qr, qi '
      deallocate(nrmarray, stat=st); if (st/=0) stop ' deallocate nrmarray '
      deallocate(ovdev, stat=st); if (st/=0) stop ' deallocate ovdev '
      deallocate(ialj, stat=st); if (st/=0) stop ' deallocate ialj '
      deallocate(jss, stat=st); if (st/=0) stop ' deallocate jss '
      deallocate(nbtt, stat=st); if (st/=0) stop ' deallocate nbtt '

      propdtpt_setup=.FALSE.
      
      end subroutine flush_propdtpt_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine init_vomake_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Allocates arrays for device vo_make module

      use gwm, only : vks, gam, seg_sh, seg_w, ft
      use gwm, only : n, nsp, seg, ngam, nt

      implicit none

!     Check global arrays for allocation status; allocate local arrays
      if (.not.allocated(vks)) stop ' vks not allocated'
      if (.not.allocated(gam)) stop ' gam not allocated'
      if (.not.allocated(seg_sh)) stop ' seg_sh not allocated'
      if (.not.allocated(seg_w)) stop ' seg_w not allocated'
      if (.not.allocated(ft)) stop ' ft not allocated'

      if (SIZE(vks,1).ne.n .or. SIZE(vks,2).ne.nsp) then
         write(*,'(X,4(A,I0),A)') 'vks must be (',n,' x ',nsp,&
                                  ') but is (',SIZE(vks,1),' x ',SIZE(vks,2),')'
         stop ' init_vomake_device(): wrong size: vks'
      endif
      if (SIZE(gam,1).ne.seg .or. SIZE(gam,2).ne.ngam) then
         write(*,'(X,4(A,I0),A)') 'gam must be (',seg,' x ',ngam,&
                                  ') but is (',SIZE(gam,1),' x ',SIZE(gam,2),')'
         stop ' init_vomake_device(): wrong size: gam'
      endif
      if (SIZE(seg_sh).ne.ngam) then
         write(*,'(X,2(A,I0),A)') 'seg_sh must be (',ngam,') but is (',SIZE(seg_sh),')'
         stop ' init_vomake_device(): wrong size: seg_sh'
      endif
      if (SIZE(seg_w).ne.ngam) then
         write(*,'(X,2(A,I0),A)') 'seg_w must be (',ngam,') but is (',SIZE(seg_w),')'
         stop ' init_vomake_device(): wrong size: seg_w'
      endif
      if (SIZE(ft,1).ne.(nt+1) .or. SIZE(ft,2).ne.ngam) then
         write(*,'(X,4(A,I0),A)') 'ft must be (',nt+1,' x ',gam,&
                                  ') but is (',SIZE(ft,1),' x ',SIZE(ft,2),')'
         stop ' init_vomake_device(): wrong size: ft'
      endif

!     Used by both 8_vr.f90 and 8_ct.f90:
      !$acc enter data copyin(vks,gam,seg_sh,seg_w) &
      !$acc create(ft)

      vomake_setup=.TRUE.

      end subroutine init_vomake_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine flush_vomake_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Deallocates arrays for device vomake module

      use gwm, only : vks, gam, seg_sh, seg_w, ft

      implicit none

!     Used by both 8_vr.f90 and 8_ct.f90:
      !$acc exit data delete(vks,gam,seg_sh,seg_w,ft)

      vomake_setup=.FALSE.

      end subroutine flush_vomake_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine init_ffts_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initialize cuFFT plans

      use gwm, only: n,ns,nx,ny,nz,nt,nw,ngam,scale_vh

      implicit none
      integer(4) :: rnk,rnk1,cerr
      integer    :: nb,istat
      integer    :: fftsize(3),fftsizebig(3)
      integer, pointer :: dummy=Null()

      rnk = 3
      rnk1 = 1
      fftsize = (/nz,ny,nx/)
      fftsizebig = (/nz*scale_vh,ny*scale_vh,nx*scale_vh/)
      nb = n*(scale_vh)**3

      istat = cudaStreamCreateWithFlags(stream_g, cudaStreamNonBlocking)
      call acc_set_cuda_stream(0, stream_g)

!     FFT plan (many vsn) for pteta ('ns' cols for pt +  1 col for eta)
      cerr = cufftPlanMany(cufft_plan_pteta,rnk,fftsize,dummy,1,n,dummy,1,n,CUFFT_Z2Z,ns+1)
      if (cerr/=0) stop 'Error creating cuFFTplanMany (pteta)'
      cerr = cufftSetStream(cufft_plan_pteta,stream_g)
      if (cerr/=0) stop ' problem in cufft cufftSetStream for cufft_plan_pteta '

!     FFT plan (many vsn) for rho
      cerr = cufftPlanMany(cufft_plan_rho,rnk,fftsizebig,dummy,1,nb,dummy,1,nb,CUFFT_Z2Z,2)
      if (cerr/=0) stop 'Error creating cuFFTplanMany (rho)'
      cerr = cufftSetStream(cufft_plan_rho,stream_g)
      if (cerr/=0) stop ' problem in cufft cufftSetStream for cufft_plan_rho '

!     FFT plan (many vsn) for propagation ('2*ns' or '2' for parallel propagation)
      cerr = cufftPlanMany(cufft_plan(1),rnk,fftsize,dummy,1,n,dummy,1,n,CUFFT_Z2Z,2*ns)
      if (cerr/=0) stop 'Error creating cuFFTplanMany (propagation-1)'
      cerr = cufftSetStream(cufft_plan(1),stream_g)
      if (cerr/=0) stop ' problem in cufft cufftSetStream for cufft_plan-1 '
      cerr = cufftPlanMany(cufft_plan(2),rnk,fftsize,fftsize,1,n*ns,fftsize,1,n*ns,CUFFT_Z2Z,2)
      if (cerr/=0) stop 'Error creating cuFFTplanMany (propagation-2)'
      cerr = cufftSetStream(cufft_plan(2),stream_g)
      if (cerr/=0) stop ' problem in cufft cufftSetStream for cufft_plan-2 '

!     FFT plan (many vsn) for time-ordering ('cnk' is batch size)
      cerr = cufftPlanMany(cufft_plan_ft,rnk1,nw,dummy,1,nw,dummy,1,nw,CUFFT_C2C,cnk) ! C2C or Z2Z
      if (cerr/=0) stop 'Error creating cuFFTplanMany (ft)'
      cerr = cufftSetStream(cufft_plan_ft,stream_g)
      if (cerr/=0) stop ' problem in cufft cufftSetStream for cufft_plan_ft '

      ffts_setup=.TRUE.

      end subroutine init_ffts_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine flush_ffts_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Destroy cuFFT plans

      implicit none
      integer(4) :: cerr

!     FFT plan (many vsn) for filtering
      cerr = cufftDestroy(cufft_plan_pteta)
      if (cerr/=0) stop 'Error flushing cuFFTplanMany (pteta)'

!     FFT plan (many vsn) for rho
      cerr = cufftDestroy(cufft_plan_rho)
      if (cerr/=0) stop 'Error flushing cuFFTplanMany (rho)'

!     FFT plan (many vsn) for propagation
      cerr = cufftDestroy(cufft_plan(1))
      if (cerr/=0) stop 'Error flushing cuFFTplanMany (propagation-1)'
      cerr = cufftDestroy(cufft_plan(2))
      if (cerr/=0) stop 'Error flushing cuFFTplanMany (propagation-2)'

!     FFT plan (many vsn) for time-ordering
      cerr = cufftDestroy(cufft_plan_ft)
      if (cerr/=0) stop 'Error flushing cuFFTplanMany (ft)'

      ffts_setup=.FALSE.

      end subroutine flush_ffts_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module device_mem_module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
