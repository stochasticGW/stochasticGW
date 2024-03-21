!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE FFTINTERPOLATE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Module for interpolating from one grid to another via FFT

      USE PROCOPTIONS
      USE UTILS
      USE, intrinsic :: iso_c_binding

      implicit none
      include 'fftw3.f03_include'

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine fft_interpolate3(v,w,nv,nw)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Interpolate from v-grid to w-grid. v,w must be allocated before calling

      implicit none
      real*8, intent(in)    :: v(:)
      real*8, intent(inout) :: w(:)
      integer, intent(in)   :: nv(3),nw(3)
      real, allocatable     :: w3(:,:,:),waux3(:,:,:)
      complex*16, allocatable :: vaux(:), waux(:)
      integer :: nvh(3),nwb(3),fac(3)
      integer :: i,j,k,idx,jdx,kdx
      real*8  :: norm
      logical :: debug=.false.

!     All dimensions equal: copy w <- v and exit
      if (ALL(nv.eq.nw)) then
         w(:)=v(:)
         return
      endif

!     Find grid-length scaling factors and scaled grid size
      do i=1,3
         nvh(i)=nv(i)/2
         fac(i)=1
         do
           if (nw(i)*fac(i) .ge. nv(i)) exit
           fac(i)=2*fac(i)
         enddo
         nwb(i)=nw(i)*fac(i)
      enddo

!     Forward transform 'vaux' spatial -> frequency
      allocate(vaux(PRODUCT(nv)))
      vaux(:)=v(:)
      call fft3d_forward(vaux,nv)

      if (debug) then
         write(*,*) 'Unpadded Fourier coefs'
         do i=1,size(vaux)
            write(*,*) i,vaux(i)
         enddo
      endif

!     Copy Fourier coefs 'vaux' -> 'waux'
      allocate(waux(PRODUCT(nwb)))
      waux(:)=(0.d0,0.d0)

      call zero_pad_dcube(vaux,waux,nv,nwb)

      if (debug) then
         write(*,*) 'Padded Fourier coefs'
         do i=1,size(waux)
            write(*,*) i,waux(i)
         enddo
      endif

!     Reverse transform 'waux' frequency -> spatial
      call fft3d_reverse(waux,nwb)

      if (debug) then
         write(*,*) 'Extended-grid values'
         do i=1,size(waux)
            write(*,*) i,waux(i)
         enddo
      endif

!     Copy (possibly-downsampled) values 'waux' -> w
      allocate(w3(nw(1),nw(2),nw(3)),waux3(nwb(1),nwb(2),nwb(3)))
      waux3=reshape(REAL(waux),(/nwb(1),nwb(2),nwb(3)/))
      do k=1,nw(3)
         kdx=(k-1)*fac(3)+1
         do j=1,nw(2)
            jdx=(j-1)*fac(2)+1
            do i=1,nw(1)
               idx=(i-1)*fac(1)+1
               w3(i,j,k)=waux3(idx,jdx,kdx)
            enddo
         enddo
      enddo
      w=reshape(w3,(/nw(1)*nw(2)*nw(3)/))

      if (debug) then
         write(*,*) 'Downsampled values (unnormalized)'
         do i=1,size(w)
            write(*,*) i,w(i)
         enddo
      endif

      norm=1.d0/REAL(PRODUCT(nv))
      w(:)=norm*w(:)

      deallocate(vaux,waux,w3,waux3)

      end subroutine fft_interpolate3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine fft_interpolate1(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Interpolate from v-grid to w-grid. v,w must be allocated before calling

      implicit none
      real*8, intent(in)    :: v(:)
      real*8, intent(inout) :: w(:)
      complex*16, allocatable :: vaux(:), waux(:)
      integer :: i,nv,nw,nvh,nvhh,fac
      logical :: debug=.false.

      nv=SIZE(v)
      nw=SIZE(w)
      nvh=nv/2
      nvhh=(nv+1)/2

!     Same size: just copy input data
      if (nv.eq.nw) then
         w(1:nw)=v(1:nv)
         return
      endif

!     Find the grid-length scaling factor
      fac=1
      do
        if (nw*fac .ge. nv) exit
        fac=2*fac
      enddo

!     Forward transform 'vaux' spatial -> frequency
      allocate(vaux(nv))
      vaux(1:nv)=v(1:nv)
      call fft1d_forward(vaux)

      if (debug) then
         write(*,*) 'Unpadded Fourier coefs:'
         do i=1,nv
            write(*,*) i,vaux(i)
         enddo
      endif

      allocate(waux(nw*fac))
      waux(:)=(0.d0,0.d0)
      do i=1,nvh+1
         waux(i)=vaux(i)
      enddo
      do i=1,nvh
         waux(nw*fac+1-i)=vaux(nv+1-i)
      enddo

!     For even-length vaux, symmetrize Nyquist frequency in waux
      if (mod(nv,2).eq.0) then
         waux(nvh+1)=0.5d0*waux(nvh+1)
         waux(nw*fac+1-nvh)=0.5d0*waux(nw*fac+1-nvh)
      endif

      if (debug) then
         write(*,*) 'Padded Fourier coefs (extended):'
         do i=1,nw*fac
            write(*,*) i,waux(i)
         enddo
      endif

!     Reverse transform 'waux' frequency -> spatial
      call fft1d_reverse(waux)

      if (debug) then
         write(*,*) 'Extended values:'
         do i=1,nw*fac
            write(*,*) i,waux(i)
         enddo
      endif

!     Copy (possibly-downsampled) values 'waux' -> w
      do i=1,nw
         w(i)=waux((i-1)*fac+1)
      enddo

      if (debug) then
              write(*,*) 'Downsampled values (unnormalized):'
         do i=1,nw
            write(*,*) i,w(i)
         enddo
      endif

      w(:)=1/REAL(nv)*w(:)

      if (debug) then
              write(*,*) 'Downsampled values (normalized):'
         do i=1,nw
            write(*,*) i,w(i)
         enddo
      endif

      deallocate(vaux,waux)

      end subroutine fft_interpolate1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine fft3d_forward(v,vdim)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for FFTW 3D forward transform

      implicit none
      complex*16, intent(inout) :: v(:)
      integer, intent(in) :: vdim(:)
      integer*8 :: plan
      integer   :: d

      d=SIZE(vdim)
      if (d.ne.3) then
         write(*,*) 'fft3d_forward(): d must be 3 but is ',d
         stop
      endif
      
      call dfftw_plan_dft_3d(plan,vdim(1),vdim(2),vdim(3),&
                             v,v,FFTW_FORWARD,FFTW_ESTIMATE)
      call dfftw_execute_dft(plan,v,v)
      call dfftw_destroy_plan(plan)

      end subroutine fft3d_forward

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine fft3d_reverse(v,vdim)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for FFTW 3D reverse transform

      implicit none
      complex*16, intent(inout) :: v(:)
      integer, intent(in) :: vdim(:)
      integer*8 :: plan
      integer   :: d

      d=SIZE(vdim)
      if (d.ne.3) then
         write(*,*) 'fft3d_reverse(): d must be 3 but is ',d
         stop
      endif

      call dfftw_plan_dft_3d(plan,vdim(1),vdim(2),vdim(3),&
                             v,v,FFTW_BACKWARD,FFTW_ESTIMATE)
      call dfftw_execute_dft(plan,v,v)
      call dfftw_destroy_plan(plan)

      end subroutine fft3d_reverse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine fft1d_forward(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for FFTW 1D forward transform

      implicit none
      complex*16, intent(inout) :: v(:)
      integer*8  :: plan
      integer    :: n

      n=SIZE(v)
      call dfftw_plan_dft_1d(plan,n,v,v,FFTW_FORWARD,FFTW_ESTIMATE)
      call dfftw_execute_dft(plan,v,v)
      call dfftw_destroy_plan(plan)

      end subroutine fft1d_forward

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine fft1d_reverse(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for FFTW 1D reverse transform

      implicit none
      complex*16, intent(inout) :: v(:)
      integer*8  :: plan
      integer    :: n

      n=SIZE(v)
      call dfftw_plan_dft_1d(plan,n,v,v,FFTW_BACKWARD,FFTW_ESTIMATE)
      call dfftw_execute_dft(plan,v,v)
      call dfftw_destroy_plan(plan)

      end subroutine fft1d_reverse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine zero_pad_dcube(v,w,vdim,wdim)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copy Fourier components from multi-indexed array v to expanded multi-
! indexed array w, overwriting w. v,w must be allocated before calling.

      implicit none
      complex*16, intent(in)    :: v(:)
      complex*16, intent(inout) :: w(:)
      integer, intent(in)  :: vdim(:),wdim(:)
      integer, allocatable :: c(:),csign(:),vh(:),vvals(:),vindex(:)
      integer, allocatable :: stridev(:),stridew(:)
      integer, allocatable :: vrange(:,:),wrange(:,:)
      integer, allocatable :: v1range(:,:),w1range(:,:)
      real*8, allocatable  :: nymul(:)
      integer :: h,i,j,ndim,ncnk,imod,nmod,v1d,w1d
      real*8  :: nyfac
      logical, parameter :: debug=.false.

!     Number of hypercubes to copy, which includes beginning and end of
!     of range for each degree of freedom (DOF)
      ndim=SIZE(vdim)
      ncnk=2**ndim

!     Error checking
      if (SIZE(wdim).ne.ndim) then
         write(*,*) 'copy_dcube(): vdim, wdim have ',ndim,SIZE(wdim),&
         ' dimensions, respectively; must match'
         stop 'v,w must have same dimensions'
      endif
      do j=1,ndim
         if (wdim(j).lt.vdim(j)) then
            write(*,*) 'copy_dcube(): mode ',j,': wdim(j) = ',wdim(j),&
            ' < vdim(j) = ',vdim(j)
            stop
         endif
      enddo
      if (PRODUCT(vdim).ne.SIZE(v)) then
         write(*,*) &
         'copy_dcube(): dimensions of v are inconsistent with vdim'
         stop
      endif
      if (PRODUCT(wdim).ne.SIZE(w)) then
         write(*,*) &
         'copy_dcube(): dimensions of w are inconsistent with wdim'
         stop
      endif

      allocate(c(ndim),csign(ndim),vh(ndim),vvals(ndim),vindex(ndim))
      allocate(stridev(ndim),stridew(ndim),nymul(ndim))
      allocate(vrange(ndim,2),wrange(ndim,2))
      allocate(v1range(ndim,2),w1range(ndim,2))
      vh(:)=vdim(:)/2

!     Strides for 1D arrays
      stridev(1)=1
      stridew(1)=1
      do j=2,ndim
         stridev(j)=stridev(j-1)*vdim(j-1)
         stridew(j)=stridew(j-1)*wdim(j-1)
      enddo

!     Multiplier for symmetrizing Nyquist frequency
!     (1.0 for odd dimensions, 0.5 for even)
      do j=1,ndim
         nymul(j)=1.d0/REAL(-MOD(vdim(j),2)+2)
      enddo

      if (debug) then
         write(*,*) 'ichunk, which mode segment'
         write(*,*) 'vdim :',(vdim(j),j=1,ndim)
         write(*,*) 'wdim :',(wdim(j),j=1,ndim)
         write(*,*) 'nymul:',(nymul(j),j=1,ndim)
      endif

      w(:)=(0.d0,0.d0)

      do i=1,ncnk

!       Determine ranges of multi-dimensional hypercube to copy
        imod=i-1
        nmod=ncnk
        do j=ndim,1,-1
           nmod=nmod/2
           c(j)=imod/nmod
           imod=MOD(imod,nmod)

           ! For even sizes, the ranges are selected so that the Nyquist
           ! frequency is divided by 2 and copied into 2 slots in w
           vvals(j)=vh(j)+1-c(j)
           csign(j)=(-1)**c(j)
           vrange(j,1)=c(j)*(vdim(j)-1)
           vrange(j,2)=vrange(j,1)+csign(j)*(vvals(j)-1)
           wrange(j,1)=c(j)*(wdim(j)-1)
           wrange(j,2)=wrange(j,1)+csign(j)*(vvals(j)-1)

           do h=1,2
              v1range(j,h)=vrange(j,h)*stridev(j)
              w1range(j,h)=wrange(j,h)*stridew(j)
           enddo
        enddo

        if (ANY(vvals.eq.0)) cycle

        if (debug) then
           write(*,'(I4,X,3(X,I1))') i, (c(j),j=1,ndim)
           write(*,'(A,I2,5(A,3(X,I2)),A)')  &
                'Chunk ',i,': copy v(', &
                (vrange(j,1),j=1,ndim),'):v(',(vrange(j,2),j=1,ndim),&
                ') -> w(',(wrange(j,1),j=1,ndim),'):w(',&
                (wrange(j,2),j=1,ndim),'): ',(vvals(j),j=1,ndim),' dims'
        endif

!       Initialize index counters
        nyfac=1.d0
        vindex(:)=1
        vindex(1)=0
        v1d=v1range(1,1)-csign(1)*stridev(1)+1
        w1d=w1range(1,1)-csign(1)*stridew(1)+1
        do j=2,ndim
           v1d=v1d+v1range(j,1)
           w1d=w1d+w1range(j,1)
           if (vindex(j).eq.vvals(j)) nyfac=nyfac*nymul(j)
        enddo

!       Multidimensional hypercube copy
        do
!          Find the DOF to increment
           do j=1,ndim
              h=j
              if (vindex(j).lt.vvals(j)) exit
           enddo

!          Exit when last DOF has reached max index value
           if (h.eq.ndim .and. vindex(h).eq.vvals(h)) exit

!          Reset earlier DOFs; increment first DOF not at limit
           do j=1,h-1
              vindex(j)=1
              nyfac=nyfac/nymul(j)
              v1d=v1d-v1range(j,2)+v1range(j,1)
              w1d=w1d-w1range(j,2)+w1range(j,1)
           enddo
           vindex(h)=vindex(h)+1
           if (vindex(h).eq.vvals(h)) nyfac=nyfac*nymul(h)
           v1d=v1d+csign(j)*stridev(h)
           w1d=w1d+csign(j)*stridew(h)
           w(w1d)=w(w1d)+v(v1d)*nyfac !!! Power of 0.5 for Nyquist-even

           if (debug) then
              write(*,'(A,2(3(X,I0),A),2(I0,A))') &
              'Copying v(',&
              (vrange(j,1)+csign(j)*(vindex(j)-1),j=1,ndim),&
              ') -> w(',(wrange(j,1)+csign(j)*(vindex(j)-1),j=1,ndim),&
              '), i.e. copy v(',v1d,') -> w(',w1d,')'
              write(*,*) 'nyfac = ',nyfac
           endif

        enddo ! Multidimensional hypercube copy

      enddo ! segments

      deallocate(c,csign,vrange,wrange,vh,vvals,vindex)
      deallocate(stridev,stridew,v1range,w1range,nymul)

      end subroutine zero_pad_dcube

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine test_1d_interpolation()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests 1D interpolation routine

      implicit none
      integer :: i,ni,nf
      real*8  :: L,st,pi,x,ae,avgerr,maxerr
      real*8, allocatable :: v(:), w(:), wref(:)

      ni=67
      nf=32
      allocate(v(ni),w(nf),wref(nf))

      pi=dacos(-1.d0)
      L=4*pi
      st=2.d0

      write(*,*) '---------------------'
      write(*,*) '1D Interpolation Test'
      write(*,*) '---------------------'
      write(*,*)

!     Values on initial grid
      write(*,*) 'Initial values:'
      do i=1,ni
         x=st+(i-1)*L/REAL(ni)
         v(i)=test_f(x,L)
         write(*,*) i,x,v(i)
      enddo

!     Reference values on final grid
      do i=1,nf
         x=st+(i-1)*L/REAL(nf)
         wref(i)=test_f(x,L)
      enddo

!     Interpolate initial grid -> final grid
      call fft_interpolate1(v,w)

      write(*,*) 'Interpolated values (new grid):'
      avgerr=0.d0
      maxerr=0.d0
      do i=1,nf
         x=st+(i-1)*L/REAL(nf)
         ae=abs(w(i)-wref(i))
         write(*,*) i,x,w(i),wref(i),ae
         avgerr=avgerr+ae
         if (ae.gt.maxerr) maxerr=ae
      enddo
      avgerr=avgerr/REAL(nf)

      write(*,*) 'Avg Error = ',avgerr,'; Max Error = ',maxerr
      write(*,*)

      deallocate(v,w,wref)

      end subroutine test_1d_interpolation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine test_1dN_interpolation()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests 1D interpolation routine using band-limited test function

      implicit none
      integer :: i,ni,nf,nt
      real*8  :: L,st,pi,x,ae,avgerr,maxerr
      real*8, allocatable :: v(:), w(:), wref(:)

      ni=16
      nf=9
      nt=ni
      allocate(v(ni),w(nf),wref(nf))

      pi=dacos(-1.d0)
      L=4*pi
      st= 0.d0

      write(*,*) '---------------------'
      write(*,*) '1D Interpolation Test'
      write(*,*) '(from Fourier coefs)'
      write(*,*) '---------------------'
      write(*,*)

!     Values on initial grid
      write(*,*) 'Initial values:'
      do i=1,ni
         x=st+(i-1)*L/REAL(ni)
         v(i)=test_fN(x,L,nt)
         write(*,*) i,x,v(i)
      enddo

!     Reference values on final grid
      do i=1,nf
         x=st+(i-1)*L/REAL(nf)
         wref(i)=test_fN(x,L,nt)
      enddo

!     Interpolate initial grid -> final grid
      call fft_interpolate1(v,w)

      write(*,*) 'Interpolated values (new grid):'
      avgerr=0.d0
      maxerr=0.d0
      do i=1,nf
         x=st+(i-1)*L/REAL(nf)
         ae=abs(w(i)-wref(i))
         write(*,*) i,x,w(i),wref(i),ae
         avgerr=avgerr+ae
         if (ae.gt.maxerr) maxerr=ae
      enddo
      avgerr=avgerr/REAL(nf)

      write(*,*) 'Avg Error = ',avgerr,'; Max Error = ',maxerr
      write(*,*)

      deallocate(v,w,wref)

      end subroutine test_1dN_interpolation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine test_3d_interpolation()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests 3D interpolation routine using band-limited test function

      implicit none
      integer :: i,j,k,m,ni(3),nf(3)
      real*8  :: L,st,pi,x,y,z,tx,ty,tz,ae,avgerr,maxerr
      real*8, allocatable :: v(:), w(:), wref(:)

      ni=(/16,15,12/)
      nf=(/6,10,6/)
      allocate(v(PRODUCT(ni)),w(PRODUCT(nf)),wref(PRODUCT(nf)))

      pi=dacos(-1.d0)
      L=4*pi
      st=2.d0

      write(*,*) '---------------------'
      write(*,*) '3D Interpolation Test'
      write(*,*) '(from Fourier coefs)'
      write(*,*) '---------------------'
      write(*,*)

!     Values on initial grid
!      write(*,*) 'Initial values:'
      m=1
      do k=1,ni(3)
         z=(k-1)*L/REAL(ni(3))
         tz=test_fN(z,L,ni(3))
         do j=1,ni(2)
            y=(j-1)*L/REAL(ni(2))
            ty=test_fN(y,L,ni(2))
            do i=1,ni(1)
               x=(i-1)*L/REAL(ni(1))
               tx=test_fN(x,L,ni(1))
               v(m)=tx*ty*tz !tx+ty+tz !+ tx*ty*tz
!               write(*,*) m,z,y,x,v(m)
               m=m+1
            enddo
         enddo
      enddo

!     Reference values on final grid
      m=1
      do k=1,nf(3)
         z=(k-1)*L/REAL(nf(3))
         tz=test_fN(z,L,ni(3))
         do j=1,nf(2)
            y=(j-1)*L/REAL(nf(2))
            ty=test_fN(y,L,ni(2))
            do i=1,nf(1)
               x=(i-1)*L/REAL(nf(1))
               tx=test_fN(x,L,ni(1))
               wref(m)=tx*ty*tz
               m=m+1
            enddo
         enddo
      enddo

!     Interpolate initial grid -> final grid
      call fft_interpolate3(v,w,ni,nf)

!      write(*,*) 'Interpolated values (new grid):'
      avgerr=0.d0
      maxerr=0.d0
      m=1
      do k=1,nf(3)
         do j=1,nf(2)
            do i=1,nf(1)
               ae=abs(w(m)-wref(m))
!               write(*,*) k,j,i,w(m),wref(m),ae
               avgerr=avgerr+ae
               if (ae.gt.maxerr) maxerr=ae
               m=m+1
            enddo
         enddo
      enddo

      avgerr=avgerr/REAL(PRODUCT(nf))

      write(*,*) 'Avg Error = ',avgerr,'; Max Error = ',maxerr
      write(*,*)

      deallocate(v,w,wref)

      end subroutine test_3d_interpolation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine test_dcube()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Test routine for zero_pad_dcube() copy routine

      implicit none
      complex*16, allocatable :: v(:), w(:)
      real*8, allocatable     :: v3(:,:,:), w3(:,:,:)
      integer, allocatable    :: vdim(:), wdim(:)
      integer :: i,j,k,l
      real*8  :: val
      character(len=64) :: frmt

      allocate(vdim(3),wdim(3))
      vdim=(/5,1,2/)
      wdim=(/8,1,2/)

      allocate(v(PRODUCT(vdim)),w(PRODUCT(wdim)))
      v(:)=(1.d0,0.d0)

      l=1
      do k=1,vdim(3)
         do j=1,vdim(2)
            do i=1,vdim(1)
               val=8.d0 !100*(k-1)+10*(j-1)+(i-1)
               v(l)=val
               l=l+1
            enddo
         enddo
      enddo

      call zero_pad_dcube(v,w,vdim,wdim)

      allocate(v3(vdim(1),vdim(2),vdim(3)),w3(wdim(1),wdim(2),wdim(3)))
      v3=reshape(REAL(v),vdim(1:3))
      w3=reshape(REAL(w),wdim(1:3))

      write(frmt,'(A,I0,A)') '(',vdim(1),'(F4.0,X))'
      write(*,*) 'v, by slices'
      write(*,*)
      do k=1,vdim(3)
         write(*,*) 'v, slice ',k
         do j=1,vdim(2)
            write(*,frmt) (v3(i,j,k),i=1,vdim(1))
         enddo
         write(*,*)
      enddo

      write(frmt,'(A,I0,A)') '(',wdim(1),'(F4.0,X))'
      write(*,*) 'w, by slices'
      write(*,*)
      do k=1,wdim(3)
         write(*,*) 'w, slice ',k
         do j=1,wdim(2)
            write(*,frmt) (w3(i,j,k),i=1,wdim(1))
         enddo
         write(*,*)
      enddo

      deallocate(v,w,vdim,wdim)

      end subroutine test_dcube

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine show_zcut(v,nxyz,dxyz,lxyz,inds)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Print cut in z for 3D cube

      implicit none
      real*8, intent(in)  :: v(:), dxyz(:), lxyz(:)
      integer, intent(in) :: nxyz(:), inds(:)
      integer :: j,k,idx,stride,offset
      real*8, allocatable :: vtest(:),wtest(:)
      real*8  :: xyz(3)

      allocate(vtest(nxyz(3)),wtest(nxyz(3)))

!     x,y position
      do j=1,2
         xyz(j)=(inds(j)-1)*dxyz(j)-lxyz(j)/2.d0
      enddo

!     For z-position in big array
      stride=nxyz(1)*nxyz(2)
      offset=inds(1)-1+(inds(2)-1)*nxyz(1)

      do k=1,nxyz(3)
         idx=offset+(k-1)*stride+1
         vtest(k)=v(idx)
         xyz(3)=(k-1)*dxyz(3)-lxyz(3)/2.d0
         write(*,*) (xyz(j),j=1,3), v(idx)
      enddo

      call fft_interpolate1(vtest,wtest)

      deallocate(vtest,wtest)

      end subroutine show_zcut

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      real*8 function test_f(x,L) result(f)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Test function for interpolation

      implicit none
      real*8, intent(in) :: x,L
      real*8  :: pi

      pi=dacos(-1.d0)

      f = sin(22*pi*x/L) &
          - 0.3d0*(cos(5*pi*x/L)**2) &
          + 3.d0*exp(-2*(x-5.d0)**2) &
          - 2.d0*exp(-4*(x-8.d0)**2) &
          + 3.d0

      end function test_f

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      real*8 function test_fN(x,L,N) result(g)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Test function for interpolation, band-limited

      implicit none
      real*8, intent(in)  :: x,L
      integer, intent(in) :: N
      integer    :: i
      real*8     :: pi,fac
      complex*16 :: t,f

      pi=dacos(-1.d0)
      f=(1.d0,0.d0)
      do i=1,N/2
         t=2*i*pi*x/L
         fac=REAL(N-i+1)/REAL(N)
         f = f + fac*((1.d0,0.d0)*cos(t) + (0.d0,1.d0)*sin(t))
         f = f + fac*((1.d0,0.d0)*cos(-t) + (0.d0,1.d0)*sin(-t))
      enddo
      if (mod(N,2).eq.0) then
         t=2*(N/2)*pi*x/L
         fac=REAL(N-N/2+1)/REAL(N)
         f = f - fac*((1.d0,0.d0)*cos(t) + (0.d0,1.d0)*sin(t))
      endif
      g=REAL(f)

      end function test_fN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE FFTINTERPOLATE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
