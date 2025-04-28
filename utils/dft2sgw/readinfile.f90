!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE READINFILE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Module for reading dft2sgw input files

      USE PROCOPTIONS
      USE UTILS
      USE iso_fortran_env

      implicit none

      TYPE DFT2SGWIN
        character*50 :: outputfile
        character*16 :: filetype
        integer :: norbs, filecase=-99
        integer, allocatable :: orbs(:)
      END TYPE DFT2SGWIN

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine read_dft2sgw(filename,theinput,theoptions)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads dft2sgwin type

      implicit none
      TYPE (DFT2SGWIN) :: theinput
      TYPE (PROCOPS)   :: theoptions
      character(len=256), intent(in) :: filename
      integer, parameter :: norbmax=1024
      integer, allocatable :: orbs(:)
      character*9 :: ch
      logical :: file_flg
      integer :: i,u,norbs,stat
      real*8  :: dummy

      allocate(orbs(norbmax))
      orbs(:)=0

      u=openfile(filename)
      read(u,*) ch, theinput%outputfile
      read(u,*) ch, theinput%filetype
      read(u,*,iostat=stat) ch,theoptions%nxyz(1),&
                            theoptions%nxyz(2),theoptions%nxyz(3)
      read(u,*,iostat=stat) ch,orbs(:)
      if (stat == iostat_end) continue

!     Filetype,filecase: for reading from different elec. str. packages
!     Set processing options here also
      if (compareStringsForEQ(theinput%filetype,'cp2k')) then
         theinput%filecase=0
         theoptions%shift=.TRUE.
         theoptions%resort=.TRUE.
         theoptions%squareroot=.FALSE.
         theoptions%normalize=.TRUE.
      elseif (compareStringsForEQ(theinput%filetype,'qe')) then
         theinput%filecase=1
         theoptions%shift=.TRUE.
         theoptions%resort=.TRUE.
         theoptions%squareroot=.TRUE.
         theoptions%normalize=.TRUE.
      elseif (compareStringsForEQ(theinput%filetype,'rmg') .or. &
              compareStringsForEQ(theinput%filetype,'rmgdft') ) then
         theinput%filecase=2
         theoptions%shift=.TRUE.
         theoptions%resort=.TRUE.
         theoptions%squareroot=.FALSE.
         theoptions%normalize=.TRUE.
      else
         write(*,*) 'Wrapper for files of type ',&
                    TRIM(ADJUSTL(theinput%filetype)),&
                   ' are not yet implemented.'
         stop
      endif

!     Count number orbitals read and copy list to type
      norbs=0
      do i=1,norbmax
         if (orbs(i).lt.1) exit
         norbs=norbs+1
      enddo
      
      theinput%norbs=norbs
      allocate(theinput%orbs(norbs))
      theinput%orbs(1:norbs)=orbs(1:norbs)
      deallocate(orbs)

      close(u)

      end subroutine read_dft2sgw

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine flush_dft2sgw(theinput)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Destroys DFT2SGWIN type

      implicit none
      TYPE (DFT2SGWIN) :: theinput

      theinput%outputfile=''
      theinput%filetype=''
      theinput%filecase=-99
      theinput%norbs=0
      if(allocated(theinput%orbs)) deallocate(theinput%orbs)

      end subroutine flush_dft2sgw

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE READINFILE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
