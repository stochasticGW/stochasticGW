!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE UTILS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Utils module

      implicit none

      CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      integer function openfile(filename)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Open file called 'filename' with checks; returns file unit

      implicit none
      character(len=256), intent(in) :: filename
      logical :: file_flg

      openfile=LookForFreeUnit()
      inquire(file=TRIM(ADJUSTL(filename)),exist=file_flg)
      if (.not. file_flg) then
         write(*,*) 'file: "',TRIM(ADJUSTL(filename)),'" not found'
         stop
      endif
      open(unit=openfile,file=TRIM(ADJUSTL(filename)))

      end function openfile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      integer function LookForFreeUnit()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Set the smallest integer equal or greater than 20 that is 
! available as unit i/o. 
! returns the integer number of the smallest free unit.


      implicit none
      logical            :: File_Opened
      integer, parameter :: UnitMax = 300

      DO LookForFreeUnit = 20, UnitMax ! Look for non-opened I/O unit
         INQUIRE( UNIT=LookForFreeUnit, OPENED=File_Opened )
         IF (.NOT. File_Opened) EXIT
      END DO

!     If an available unit has been found, use as function results
      if (LookForFreeUnit.eq.UnitMax) then
         write(*,*) "No free I/O unit available"
         stop
      endif

      end function LookForFreeUnit

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine searchstring(u, string, line)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Find a string in file u

      implicit none
      integer, intent(in) :: u
      character(len=150), intent(out) :: line
      character(len=25), intent(in)   :: string
      integer :: stat, whereis

      gdo : do while (.true.)

          read(u,'(A)',iostat=stat) line
          whereis=index(line,trim(string))
          if(whereis.ne.0) then
              exit gdo
          elseif(stat<0) then
              write(*,*) 'ERROR: string not found: "',string,&
                         '" in file unit,',u
              stop
          endif

        enddo gdo

      end subroutine searchstring

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function compareStringsForEQ( s1, s2 ) RESULT( equal_result )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares strings and does not differentiate between lower and upper
! case. Also, leading and/or trailing spaces are ignored.

      implicit none
      character(LEN=*), intent(in) :: s1, s2
      logical                      :: equal_result
      integer(4)                   :: is, offset1, offset2, L, distance
      integer(4)                   :: length1, length2
      character(LEN=1)             :: char1, char2
      logical                      :: equal

      equal = .TRUE.

!     Determine the length of each string, NOT counting empty spaces.
      length1 = LEN_TRIM(ADJUSTL(s1))
      length2 = LEN_TRIM(ADJUSTL(s2))

      IF ( length1 /= length2 ) then
         equal = .FALSE.
      ELSE

!        Left adjust each string.
         offset1 = VERIFY( s1, " " ) - 1
         offset2 = VERIFY( s2, " " ) - 1

!        Lexographic distance between 'a' and 'A' in ASCII list.
         L = ABS( ICHAR('a') - ICHAR('A') )

         do is = 1 , LEN_TRIM(s1)
!           Next characters in strings.
            char1 = s1(offset1+is:offset1+is)
            char2 = s2(offset2+is:offset2+is)

            distance = ABS( ICHAR(char1) - ICHAR(char2) )
            if ( distance /= 0 .AND. distance /= L ) then
               equal = .FALSE.
               EXIT
            endif
         enddo

      ENDIF

      equal_result = equal

      end function compareStringsForEQ

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE UTILS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

