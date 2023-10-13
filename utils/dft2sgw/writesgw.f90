!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE WRITESGW

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Module for writing 'sgwinp.txt' input for stochasticGW

      use READCUBE
      use READDFTOUT
      use UTILS

      implicit none
      logical :: header_written = .FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine write_header(u,theoutput,thecube)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes header portion of 'sgwinp.txt'

      implicit none
      TYPE (DFTOUT), intent(in) :: theoutput
      TYPE (CUBE), intent(in)   :: thecube
      integer, intent(in)       :: u
      integer :: i,j

      write(u,*) '$GEOMETRY'
      do i=1,thecube%na
         write(u,*) thecube%sym(i),(thecube%coords(i,j),j=1,3)
      enddo

      write(u,*) '$GRID'
      write(u,*) 'nx       ',thecube%nxyz(1)
      write(u,*) 'ny       ',thecube%nxyz(2)
      write(u,*) 'nz       ',thecube%nxyz(3)
      write(u,*) 'dx       ',thecube%dxyz(1)
      write(u,*) 'dy       ',thecube%dxyz(2)
      write(u,*) 'dz       ',thecube%dxyz(3)
      write(u,*) 'nsp      ',1
      write(u,*) 'homo     ',theoutput%homo
      write(u,*) 'lumo     ',theoutput%lumo
      write(u,*) '$ORBITALS'

      header_written = .TRUE.

      end subroutine write_header

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine write_orbital(u,theoutput,thecube)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes orbital/density portion of 'sgwinp.txt'

      implicit none
      TYPE (DFTOUT), intent(in) :: theoutput
      TYPE (CUBE), intent(in)   :: thecube
      integer, intent(in)       :: u

      if (.not.thecube%assigned) then
         write(*,*) 'ERROR: cube file not assigned.'
         stop
      endif

      if (.not.header_written) call write_header(u,theoutput,thecube)

      write(u,*) thecube%label,thecube%indx,thecube%nrg,thecube%spin
      write(u,*) thecube%phi

      end subroutine write_orbital

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE WRITESGW

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

