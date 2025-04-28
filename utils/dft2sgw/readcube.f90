!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE READCUBE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Module for reading cube files

      USE PROCOPTIONS
      USE UTILS
      USE FFTINTERPOLATE

      implicit none
      character*2 :: an(118)

      TYPE CUBE
        LOGICAL :: assigned = .false.
        INTEGER :: na,np,spin,indx
        INTEGER :: nxyz(3)
        REAL*8  :: dxyz(3)
        REAL*8  :: nrg,dv,norm
        REAL*8, ALLOCATABLE :: coords(:,:)
        REAL*8, ALLOCATABLE :: phi(:) 
        CHARACTER*2, ALLOCATABLE :: sym(:)
        CHARACTER(len=9)   :: label
        character(len=256) :: filename
      END TYPE CUBE

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine read_cubefile(filename,thecube)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads cube file, extracts and returns data

      implicit none
      TYPE (CUBE) :: thecube
      character(len=256), intent(in) :: filename
      logical :: file_flg
      integer :: i,u,za,nx,ny,nz
      real*8  :: dummy,dx,dy,dz

      call prep_symbol()

      thecube%filename=filename
      u=openfile(filename)

!     Read grid info
      read(u,*)
      read(u,*)
      read(u,*) thecube%na
      read(u,*) nx, dx
      read(u,*) ny, dummy, dy
      read(u,*) nz, dummy, dummy, dz

!     Read geometry
      allocate(thecube%sym(thecube%na),thecube%coords(thecube%na,3))
      do i=1,thecube%na
         read(u,*) za, dummy, thecube%coords(i,:)
         thecube%sym(i)=an(za)
      enddo
      thecube%np=nx*ny*nz
      thecube%dv=dx*dy*dz
      thecube%nxyz=(/nx,ny,nz/)
      thecube%dxyz=(/dx,dy,dz/)

!     Read orbital
      allocate(thecube%phi(thecube%np))
      read(u,*) thecube%phi

      close(u)

      end subroutine read_cubefile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine assign_cube(thecube,label,indx,nrg,spin)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Destroys CUBE type

      implicit none
      TYPE (CUBE) :: thecube
      integer, intent(in) :: indx,spin
      real*8, intent(in)  :: nrg
      character(len=9), intent(in) :: label

      thecube%label=label
      thecube%indx=indx
      thecube%nrg=nrg
      thecube%spin=spin
      thecube%assigned=.true.

      end subroutine assign_cube

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine process_cubefile(thecube,theoptions)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Processes CUBE type

      implicit none
      TYPE (CUBE) :: thecube
      TYPE (PROCOPS), intent(in) :: theoptions
      real*8, allocatable :: phi3(:,:,:)
      real*8, allocatable :: phitmp(:)
      real*8  :: lxyz(3),norm
      integer :: i,j

      IF (theoptions%shift) THEN
!        Shift atomic coordinates into cell if outside due to low input data precision
         lxyz(:)=dble(thecube%nxyz(:))*thecube%dxyz(:)
         do i=1,thecube%na
            do j=1,3
               thecube%coords(i,j)=thecube%coords(i,j)-lxyz(j)/2
               if (thecube%coords(i,j) < -lxyz(j)/2) &
                  thecube%coords(i,j)=thecube%coords(i,j)+lxyz(j)
               if (thecube%coords(i,j) >  lxyz(j)/2) &
                  thecube%coords(i,j)=thecube%coords(i,j)-lxyz(j)
            enddo
         enddo
      ENDIF

      IF (theoptions%resort) THEN
!        Reshape orbital from (slowest) x-y-z (fastest)
!                          to (slowest) z-y-x (fastest) index ordering
         allocate(phi3(thecube%nxyz(1),thecube%nxyz(2),thecube%nxyz(3)))
         phi3=reshape(thecube%phi,thecube%nxyz,order=(/3,2,1/))
         thecube%phi=reshape(phi3,(/thecube%np/))
         deallocate(phi3)
         write(*,*) 'Orbital data resorted: {x-y-z} -> {z-y-x}'
      ENDIF

      IF (theoptions%squareroot.and.(thecube%indx.gt.0)) THEN
         do j=1,thecube%np
            if (thecube%phi(j)<0) then
               thecube%phi(j)=-sqrt(-thecube%phi(j))
            else
               thecube%phi(j)=sqrt(thecube%phi(j))
            endif
         enddo
      ENDIF

!     Interpolate grid if input file dims differ from cube file dims
      IF (ANY(theoptions%nxyz.ne.thecube%nxyz)) THEN
         write(*,'(X,A,7(I0,A))') &
         'Regridding state ',thecube%indx,': (',thecube%nxyz(1),',',&
         thecube%nxyz(2),',',thecube%nxyz(3),') -> (',&
         theoptions%nxyz(1),',',theoptions%nxyz(2),',',&
         theoptions%nxyz(3),')'

         allocate(phitmp(PRODUCT(theoptions%nxyz)))
         call fft_interpolate3(thecube%phi,phitmp,&
                               thecube%nxyz,theoptions%nxyz)
         deallocate(thecube%phi)
         allocate(thecube%phi(PRODUCT(theoptions%nxyz)))
         thecube%phi(:)=phitmp(:)
         deallocate(phitmp)
 
         thecube%np=1
         thecube%dv=1
         do j=1,3
            thecube%nxyz(j)=theoptions%nxyz(j)
            thecube%dxyz(j)=lxyz(j)/dble(theoptions%nxyz(j))
            thecube%np=thecube%np*theoptions%nxyz(j)
            thecube%dv=thecube%dv*thecube%dxyz(j)
         enddo
      ENDIF

      IF (theoptions%normalize.and.(thecube%indx.gt.0)) THEN
         norm=sum(thecube%phi**2)*thecube%dv
         if (thecube%assigned) then
            write(*,*) 'norm of state ',thecube%indx,' = ',norm
         else
            write(*,*) 'norm ',thecube%indx,' = ',norm
         endif
         norm=1.d0/sqrt(norm)
         thecube%phi=thecube%phi*norm
      ENDIF

!     Calculate the norm
      IF (thecube%indx.gt.0) THEN
         thecube%norm=sum(thecube%phi**2)*thecube%dv
      ELSE ! Density: normalizes to number of electrons
         thecube%norm=sum(thecube%phi)*thecube%dv
      ENDIF

      end subroutine process_cubefile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine show_cubefile(thecube)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Show stats of CUBE file

      implicit none
      TYPE (CUBE) :: thecube
      integer :: i

      write(*,*) '***********************'
      write(*,*) '***   Cube stats    ***'
      write(*,*) '***********************'
      write(*,*) 'file read: ',TRIM(ADJUSTL(thecube%filename))
      write(*,'(X,A,7X,A,10X,I0)') 'atoms',':',thecube%na
      write(*,'(X,A,3X,3(X,I8))') 'orbital dims: ',(thecube%nxyz(i),i=1,3)
      write(*,'(X,A,3X,3(X,f8.6))') 'orbital grid: ',(thecube%dxyz(i),i=1,3)
      write(*,'(10X,A,4X,I0,A,f12.10)') 'np = ',thecube%np,'; dv = ',thecube%dv
      write(*,'(X,A,f18.12)') 'orbital norm: ',thecube%norm
      IF (thecube%assigned) then
         write(*,'(X,A,3X,A,3X,I6,A,f16.8,A,I0)') &
               thecube%label,':',thecube%indx,'; E  = ',&
               thecube%nrg,'; spin = ',thecube%spin
      endif
      write(*,*)

      end subroutine show_cubefile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine flush_cubefile(thecube)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Destroys CUBE type

      implicit none
      TYPE (CUBE) :: thecube

      thecube%na=0
      thecube%np=0
      thecube%dv=0.d0
      thecube%norm=0.d0
      thecube%nxyz(:)=0
      thecube%dxyz(:)=0.d0
      thecube%label=''
      thecube%filename=''
      thecube%indx=-1
      thecube%nrg=-1.d99
      thecube%spin=0
      thecube%assigned=.false.
      if (allocated(thecube%coords)) deallocate(thecube%coords)
      if (allocated(thecube%phi)) deallocate(thecube%phi)
      if (allocated(thecube%sym)) deallocate(thecube%sym)

      end subroutine flush_cubefile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function make_cube_name(prefix,indx,filecase)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets array of symbols of the chemical elements

      implicit none

      integer, intent(in) :: indx, filecase
      character(len=150), intent(in) :: prefix
      character(len=256) :: make_cube_name

      select case (filecase)
      case(0) ! CP2K

         if (indx.gt.0) then
            write(make_cube_name,'(A,A,I5.5,A)') &
                  TRIM(ADJUSTL(prefix)),'-WFN_',&
                  indx,'_1-1_0.cube'
         else
            write(make_cube_name,'(A,A)') TRIM(ADJUSTL(prefix)),&
                  '-ELECTRON_DENSITY-1_0.cube'
         endif

      case(1) ! Quantum Espresso

         if (indx.gt.0) then
            write(make_cube_name,'(A,I0)') 'tmp.pp_orb',indx
         else
            write(make_cube_name,'(A)') 'tmp.pp_dens'
         endif

      case(2) ! RMGDFT

         if (indx.gt.0) then
            write(make_cube_name,'(A,I0,A)') 'kpt0_mo',indx-1,'.cube'
         else
            write(make_cube_name,'(A)') 'density.cube'
         endif

      case default
         write(*,'(A,I0,A)') 'Wrapper for files of type ',filecase,&
         ' are not yet implemented.'
         stop
      end select

      end function make_cube_name

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine prep_symbol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets array of symbols of the chemical elements

      implicit none
      an( 1: 2)=&
        (/'H ',                                                                                'He'/)
      an( 3:10)=&
        (/'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne'/)
      an(11:18)=&
        (/'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar'/)
      an(19:36)=&
        (/'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr'/)
      an(37:54)=&
        (/'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe'/)
      an(55:56)=&
        (/'Cs','Ba'/)
      an(72:86)=&
                   (/'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn'/)
      an(87:88)=&
        (/'Fr','Ra'/)
      an(104:118)=     (/'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/)

      an(57:71)=       (/'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'/)
      an(89:103)=      (/'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'/)

      end subroutine prep_symbol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE READCUBE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
