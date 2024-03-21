!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      program DFT2SGW

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Converts DFT wavefunctions/density to StochasticGW format

      use READINFILE
      use READDFTOUT
      use READCUBE
      use PROCOPTIONS
      use UTILS
      use WRITESGW
      use FFTINTERPOLATE

      implicit none
      TYPE (DFT2SGWIN)   :: theinput
      TYPE (DFTOUT)      :: theoutput
      TYPE (CUBE)        :: thecube
      TYPE (PROCOPS)     :: theoptions
      character(len=256) :: filename
      integer :: i,u

!!! TESTING
!      call test_1d_interpolation()
!      call test_1dN_interpolation()
!      call test_dcube()
!      call test_3d_interpolation()
!!! END TESTING

!     Read input file for dft2sgw
      filename='dft2sgw.in'
      write(*,*) 'reading ',TRIM(ADJUSTL(filename)),' ...'
      call read_dft2sgw(filename,theinput,theoptions)

!     Read output for the electronic structure code of choice
!     Read HOMO/LUMO energies and occupancy from DFT output
      filename=TRIM(ADJUSTL(theinput%outputfile))
      write(*,*) 'reading ',TRIM(ADJUSTL(filename)),' ...'
      call read_dftout(filename,theoutput,theinput%filecase)

!     Read orbitals from cube file; write data to sgwinp.txt
      u=LookForFreeUnit()
      open(u,file='sgwinp.txt')

      do i=1,theinput%norbs
         filename = make_cube_name(theoutput%prefix,theinput%orbs(i),&
                                   theinput%filecase)
         write(*,*) 'reading ',TRIM(ADJUSTL(filename)),' ...'
         call read_cubefile(filename,thecube)
         call assign_cube(thecube,'ORB      ',theinput%orbs(i),&
                          theoutput%nrg(theinput%orbs(i)),1)
         call process_cubefile(thecube,theoptions)
         call show_cubefile(thecube)
         call write_orbital(u,theoutput,thecube)
         call flush_cubefile(thecube)
      enddo

!     Read density and write to sgwinp.txt
      filename = make_cube_name(theoutput%prefix,0,theinput%filecase)
      write(*,*) 'reading ',TRIM(ADJUSTL(filename)),' ...'
      call read_cubefile(filename,thecube)
      call assign_cube(thecube,'DENS     ',0,0.d0,0)
      call process_cubefile(thecube,theoptions)
      call show_cubefile(thecube)
      call write_orbital(u,theoutput,thecube)
      call flush_cubefile(thecube)

!     Clean up
      call flush_dftout(theoutput)
      call flush_dft2sgw(theinput)

      close(u)

      end program DFT2SGW

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
