!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE READDFTOUT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Module for reading cp2k output files

      USE UTILS
      USE iso_fortran_env

      implicit none

      TYPE DFTOUT
        integer :: norb, nocc
        real*8  :: homo, lumo
        integer, allocatable :: occ(:)
        real*8, allocatable  :: nrg(:)
        character(len=150)   :: prefix
      END TYPE DFTOUT

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine read_dftout(filename,theoutput,filecase)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (DFTOUT) :: theoutput
      character(len=256), intent(in) :: filename
      integer, intent(in) :: filecase

      select case(filecase)

      case(0) ! CP2K
         call read_cp2kout(filename,theoutput)
      case(1) ! QE
         call read_qeout(filename,theoutput)
      case default
         write(*,'(A,I0,A)') 'Wrapper for files of type ',filecase,&
         ' are not yet implemented.'
         stop
      end select

      end subroutine read_dftout

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine read_cp2kout(filename,theoutput)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads cp2k output and fills DFTOUT type

      implicit none
      TYPE (DFTOUT) :: theoutput       
      character(len=256), intent(in) :: filename
      character(len=150) :: line
      character*50       :: s1
      character*3        :: dummy_c
      logical :: file_flg
      integer :: i,u,dummy_i,nocc,norb
      real*8  :: dummy_r, occ_r, energy

      u=openfile(filename)

!     Get prefix
      s1="GLOBAL| Project name"
      call searchstring(u,s1,line)
      theoutput%prefix = get_c_field(line,s1)

      write(*,*) 'prefix = ',TRIM(ADJUSTL(theoutput%prefix))

!     Get number of occupied and total orbitals
      s1="Number of occupied orbitals:"
      call searchstring(u,s1,line)
      nocc = get_i_field(line,s1)
      s1="Number of molecular orbitals:"
      call searchstring(u,s1,line)
      norb = get_i_field(line,s1)

      allocate(theoutput%nrg(norb),theoutput%occ(norb))
      theoutput%norb=norb
      theoutput%nocc=nocc
      theoutput%nrg(:)=-1.d-99
      theoutput%occ(:)=-1

!     Read number of states
      s1='*** SCF run converged in'
      call searchstring(u,s1,line)
      s1='MO| EIGENVALUES AND OCCUPATION NUMBERS'
      call searchstring(u,s1,line)
      read(u,*)
      read(u,*)
      do i=1,norb
         read(u,*) dummy_c , dummy_i, energy, dummy_r, occ_r
         theoutput%nrg(i)=energy
         theoutput%occ(i)=NINT(occ_r)
      enddo

      theoutput%homo=theoutput%nrg(nocc)
      theoutput%lumo=theoutput%nrg(nocc+1)

      close(u)

!     Need at least one empty orbital for HOMO-LUMO-gap
      if (norb.le.nocc) then
         write(*,'(X,2(A,I0),A)') 'norb (',norb,&
           ') must be greater than nocc (',nocc,')'
         write(*,*) 'Make sure your CP2K input contains the following:'
         write(*,*) '&GLOBAL -> PRINT_LEVEL  MEDIUM'
         write(*,*) '&FORCE_EVAL -> &DFT -> &PRINT -> &MO -> EIGENVALUES'
         write(*,*) '&FORCE_EVAL -> &DFT -> &PRINT -> &MO -> OCCUPATION_NUMBERS'
         write(*,*) '&FORCE_EVAL -> &DFT -> &SCF -> ADDED_MOS 1 1'
         stop
      endif
      theoutput%homo=theoutput%nrg(nocc)
      theoutput%lumo=theoutput%nrg(nocc+1)

      end subroutine read_cp2kout

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine read_qeout(filename,theoutput)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads Quantum Espresso output and fills DFTOUT type

      implicit none
      TYPE (DFTOUT) :: theoutput
      character(len=256), intent(in) :: filename
      real*8, parameter  :: eV2au = 0.03674932217565499
      character(len=150) :: line,subline
      character*50       :: s1
      integer :: i, u, norb, nocc, whereis

      u=openfile(filename)

!     Read number of states   
      s1 = 'number of Kohn-Sham states='
      call searchstring(u,s1,line)
      norb = get_i_field(line,s1)

      allocate(theoutput%nrg(norb),theoutput%occ(norb))
      theoutput%norb=norb
      theoutput%nrg(:)=-1.d-99
      theoutput%occ(:)=-1

!     Read energies from QE output
      s1 = 'End of self-consistent'
      call searchstring(u,s1,line)
      do i=1,3
          read(u,*)
      enddo
      read(u,*) theoutput%nrg

!     Read HOMO/LUMO energies
      s1 = 'highest occupied, lowest unoccupied level (ev):'
      call searchstring(u,s1,line)
      write(*,*) line
      whereis=index(line,TRIM(ADJUSTL(s1)))+len(TRIM(ADJUSTL(s1)))
      read(line(1+whereis:),*) theoutput%homo,theoutput%lumo

!     Prefix not used in QE vsn:
      s1 = " Writing all to output data dir"
      call searchstring(u,s1,line)
      subline = get_c_field(line,s1)
      theoutput%prefix =''

      close(u)

!     Find occupations by searching for HOMO in energy array
      nocc=0
      do i=1,norb
         if (theoutput%nrg(i).eq.theoutput%homo) then
            nocc=i
            exit
         endif
      enddo

!     Make sure enough orbitals have been computed
      if (nocc.eq.0) then
         write(*,*) 'ERROR: energy for HOMO not found'
         stop
      elseif (nocc.eq.norb) then
         write(*,*) 'ERROR: need at least nocc+1 orbitals in QE output'
         write(*,'(X,A,I0,A)') 'Make sure nbnd = ',nocc+1,&
         ' or higher in your QE input'
         stop
      endif

      theoutput%nocc=nocc
      theoutput%occ(:nocc)=2
      theoutput%occ(nocc+1:)=0
      theoutput%nrg=theoutput%nrg*eV2au !convert from eV to Hartree
      theoutput%homo=theoutput%homo*eV2au
      theoutput%lumo=theoutput%lumo*eV2au

      end subroutine read_qeout

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine flush_dftout(theoutput)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Destroys DFTOUT type

      implicit none
      TYPE (DFTOUT) :: theoutput

      theoutput%nocc=0
      theoutput%homo=-1.d99
      theoutput%lumo=-1.d99
      theoutput%prefix=''
      if(allocated(theoutput%nrg)) deallocate(theoutput%nrg)
      if(allocated(theoutput%occ)) deallocate(theoutput%occ)      

      end subroutine flush_dftout

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      integer function get_i_field(line,tag) result(val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds integer field value in file

        implicit none
        character(len=*), intent(in) :: line,tag
        integer :: whereis

        whereis=index(line,TRIM(ADJUSTL(tag)))+len(TRIM(ADJUSTL(tag)))
        read(line(1+whereis:),*) val

      end function get_i_field

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function get_c_field(line,tag)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds integer field value in file

        implicit none
        character(len=*), intent(in) :: line, tag
        character(len=150) :: get_c_field
        integer :: whereis

        whereis=index(line,TRIM(ADJUSTL(tag)))+len(TRIM(ADJUSTL(tag)))
        read(line(1+whereis:),*) get_c_field

      end function get_c_field

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE READDFTOUT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
