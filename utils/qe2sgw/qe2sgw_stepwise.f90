program qe2sgw
   use :: iso_fortran_env
   implicit none
   character*9   :: ch
   character(len=9) :: output
   integer       :: nxyz(3),nx,ny,nz,nstates,na,nocc
   integer       :: whereis,stat,i,j,one=1,xstat,cstat
   real*8        :: dxyz(3),dx,dy,dz
   real*8        :: norm,factor 
   real*8, allocatable :: energies(:),orb(:),geom(:,:)
   integer,parameter    :: njmax=1024
   integer              :: nj
   integer, allocatable :: orbj(:),za(:)
   logical       :: flg_bin, file_flg, first_hit, cubespresent
   logical       :: write_uni, write_wf, write_orbj, write_target, write_cnt
   character*50  :: outputfile,ppinput,prefix,qe_loc,s1
   character*3   :: ibnd
   character*2   :: an(118)
   character(len=150) :: line

   call prep_symbol

   allocate(orbj(njmax))
   orbj(:)=0

!  Setting variables      
   inquire(file="qe2sgw.in",exist=file_flg)
   if (.not. file_flg) stop 'Input file qe2sgw.in not found!'
   open(unit=130,file="qe2sgw.in")
   read(130,*) ch, outputfile
   read(130,*) ch, prefix
   read(130,'(A,A)') ch, qe_loc
   read(130,*) ch, flg_bin
   read(130,*) ch, output
   read(130,*) ch, nocc
   read(130,*,iostat=stat) ch,orbj(:)
   if (stat == iostat_end) continue 
!  Count number orbitals to read from input
   nj=0
   do i=1,njmax
      if (orbj(i).lt.1) exit
      nj=nj+1
   enddo

   write(6,*) 'outputfile   ', outputfile
   write(6,*) 'prefix       ', prefix
   write(6,*) 'qe bin exe   ', qe_loc
   write(6,*) 'flg_bin      ', flg_bin
   write(6,*) 'nocc         ', nocc
   write(6,*) 'output       ', TRIM(ADJUSTL(output))
   write(6,*) 'orb          ', (orbj(i),i=1,nj)

!  Output style
   write_orbj=.FALSE.
   if (TRIM(ADJUSTL(output)) == 'unified') then
      write_cnt=.FALSE.
      write_uni=.TRUE.
      write_wf=.FALSE.
      write_target=.FALSE.
      if (nj.lt.1) stop 'ERROR: orb list required for output=unified'
   elseif (TRIM(ADJUSTL(output)) == 'targeted') then
      write_cnt=.TRUE.
      write_uni=.FALSE.
      write_wf=.FALSE.
      write_target=.TRUE.
      if (nj.gt.0) write_orbj=.TRUE.
   elseif (TRIM(ADJUSTL(output)) == 'full') then
      write_cnt=.TRUE.
      write_uni=.FALSE.
      write_wf=.TRUE.
      write_target=.FALSE.
      if (nj.gt.0) stop 'ERROR: orb list not allowed for output=full'
   else
      write(6,*) 'ERROR: unrecognized option "',TRIM(ADJUSTL(output)),&
                 '" for "output" flag'
      write(6,*) 'Valid arguments are: "unified", "targeted", and "full"'
      stop 'wrong "output" value'
   endif

   if (.not.write_wf .and. flg_bin) then
      write(6,*) 'ERROR: set flg_bin=T only for output=full'
      stop 'wrong flg_bin value'
   endif

!  Read number of states, energies, from QE output 
   s1='number of Kohn-Sham states'
   open(unit=1,file=outputfile)
   call searchstring(s1)
   call get_nstates

   factor=0.03674932217565499
   allocate(energies(nstates))
   s1='End of self-consistent'
   call searchstring(s1)
   do i=1,3
       read(1,*)
   enddo
   read(1,*) energies
   energies=energies*factor !convert from eV to Hartree
   close(1)

!  Check for cube files of queried states
!  If any are missing, generate any missing inputs for pp.x
   cubespresent=.TRUE.
   do i=1,nstates
      if (write_wf .or. ANY(orbj.eq.i).or. &
          (write_target.and.(i.eq.nocc .or. i.eq.nocc+1))) then

         ! Cube file check
         write(outputfile,'(A,I0)') 'tmp.pp_orb',i
         inquire(file=trim(adjustl(outputfile)),exist=file_flg)
         if (.not. file_flg) then
            cubespresent=.FALSE.
            write(*,*) '* Cube file ',trim(adjustl(outputfile)),' not found!'
            ! Input for pp.x check
            write(ppinput ,'(A,I0,A)') 'pp_orb_',i,'.in'
            inquire(file=trim(adjustl(ppinput)),exist=file_flg)
            if (.not. file_flg) then
               call make_input_for_pp_orb(i)
               write(*,*) '-> generating ',trim(adjustl(ppinput)),' for pp.x'
            endif
         endif

      endif
   enddo

!  Check for density cube file; write input if needed
   if (write_uni .or. write_target) then

      write(outputfile,'(A)') 'tmp.pp_dens'
      inquire(file=trim(adjustl(outputfile)),exist=file_flg)
      if (.not. file_flg) then

         cubespresent=.FALSE.
         write(*,*) '* Cube file ',trim(adjustl(outputfile)),' not found!'
         ! Input for pp.x check
         write(ppinput ,'(A)') 'pp_dens.in'
         inquire(file=trim(adjustl(ppinput)),exist=file_flg)
         if (.not. file_flg) then 
            call make_input_for_pp_orb(0)
            write(*,*) '-> generating ',trim(adjustl(ppinput)),' for pp.x'
         endif
      endif

   endif
   if (.not.cubespresent) stop ' cube files missing!'

!  Read cube files and write outputs
   first_hit = .TRUE.
   do i=1,nstates
      if (write_wf .or. ANY(orbj.eq.i).or. & 
          (write_target.and.(i.eq.nocc .or. i.eq.nocc+1))) then

         write(outputfile,'(A,I0)') 'tmp.pp_orb',i
         call read_cubefile(outputfile)

!        First pass: allocate arrays, set variables and write headers
         if (first_hit) then
            nx=nxyz(1)
            ny=nxyz(2)
            nz=nxyz(3)
            dx=dxyz(1)
            dy=dxyz(2)
            dz=dxyz(3)
            call write_headers
            first_hit = .FALSE.
         endif

!        Process orbital
         do j=1,nx*ny*nz
            if (orb(j)<0) then
               orb(j)=-sqrt(-orb(j))
            else
               orb(j)=sqrt(orb(j))
            endif
         enddo
         norm=sum(orb**2)*dx*dy*dz
         write(*,*) 'norm of state',i,'is',norm
         norm=1.d0/sqrt(norm)
         orb=orb*norm

         if (write_uni) then
            if (ANY(orbj.eq.i)) then
               write(222,*) 'ORB      ',i,energies(i),1
               write(222,*) orb
            endif
         endif

         if (write_wf) then
            if (flg_bin) then
               write(1) i,one
               write(1) orb
            else
               write(1,*) i,one
               write(1,*) orb
            endif
         endif

         if (write_target) then
            if (i==nocc) write(122,*) orb
            if (i==nocc+1) write(123,*) orb
         endif

         if (write_orbj) then
            if(ANY(orbj.eq.i)) write(125,*) orb
         endif
       
         deallocate(orb,geom,za)
       endif ! orbital in set
   enddo

! Read and write density if needed

   if (write_uni .or. write_target) then
      outputfile=trim('tmp.pp_dens')
      call read_cubefile(outputfile)!,na,za,nxyz,dxyz,geom,orb)

      if (write_uni) then
         write(222,*) 'DENS     ',0,0.d0,0
         write(222,*) orb
      endif

      if (write_target) then
         write(124,*) orb
      endif
      deallocate(orb,geom,za)
   endif

   deallocate(energies)

contains

    subroutine read_cubefile(fnm)

        implicit none
        character*50, intent(in) :: fnm
        real*8, allocatable :: orb3d(:,:,:)
        integer :: i,k
        real*8  :: junk,lxyz(3)

        open(unit=4,file=fnm)
        rewind(4)
        read(4,*)
        read(4,*)
        read(4,*) na
        read(4,*) nxyz(1),dxyz(1)
        read(4,*) nxyz(2),junk,dxyz(2)
        read(4,*) nxyz(3),junk,junk,dxyz(3)

        lxyz(:)=dble(nxyz(:))*dxyz(:)
        allocate(geom(na,3),za(na))
!       Read and shift atoms into centered cell, making sure atoms do not fall
!       outside of cell due to low input data precision
        do i=1,na
           read(4,*) za(i),junk,geom(i,:)
           geom(i,:)=geom(i,:)-lxyz(:)/2
           do k=1,3
              if (geom(i,k)<-lxyz(k)/2) geom(i,k)=geom(i,k)+lxyz(k)
              if (geom(i,k)> lxyz(k)/2) geom(i,k)=geom(i,k)-lxyz(k)
           enddo
        enddo
        allocate(orb(nxyz(1)*nxyz(2)*nxyz(3)))
        allocate(orb3d(nxyz(1),nxyz(2),nxyz(3)))
        read(4,*) orb
        close(4)

!       Reshape orbital to SGW index ordering
        orb3d=reshape(orb,(/nxyz(1),nxyz(2),nxyz(3)/),order=(/3,2,1/))
        orb=reshape(orb3d,(/nxyz(1)*nxyz(2)*nxyz(3)/))
        deallocate(orb3d)

    end subroutine read_cubefile

    subroutine make_input_for_pp_orb(i)

      implicit none
      integer, intent(in) :: i
      integer :: funit
      character(len=150)  :: infile

      funit=i+300

      if (i.eq.0) then
         write(infile,'(A)')  'pp_dens.in'
      else
         write(infile ,'(A,I0,A)') 'pp_orb_',i,'.in'
      endif

!     Generate input file for pp.x
      open(unit=funit,file=trim(adjustl(infile)))
      write(funit,*) "&inputpp"
      write(funit,*) "prefix='"//trim(adjustl(prefix))//"'"
      if (i.gt.0) then
         write(funit,*) "plot_num=7"
         write(funit,*) "kband(1)=",i
         write(funit,*) "kband(2)=",i
         write(funit,*) "kpoint=1"
         write(funit,*) "lsign=.true."
      else
         write(funit,*) "plot_num=0"
      endif
      write(funit,*) "/"
      write(funit,*) "&plot"
      write(funit,*) "iflag=3"
      write(funit,*) "output_format=6"
      if (i.eq.0) then
         write(funit,*) "fileout='tmp.pp_dens'"
      else
         write(funit,'(X,A,I0,A)') "fileout='tmp.pp_orb",i,"'"
      endif

      write(funit,*) "/"
      close(funit)

    end subroutine make_input_for_pp_orb

    subroutine write_headers
 
      implicit none
      integer :: i

      if (write_uni) then
         open(unit=222,file='sgwinp.txt')
         write(222,*) '$GEOMETRY'
      endif
      if (write_cnt) then
         open(unit=2,file='cnt.ini')
      endif
      do i=1,na
         if (write_uni) write(222,*) an(za(i)),geom(i,:)
         if (write_cnt) write(2,*) an(za(i)),geom(i,:)
      enddo
      if (write_cnt) close(2)

!     Write header of orbital(s)-containing file
      if (write_uni) then
         write(222,*) '$GRID'
         write(222,*) 'nx       ',nx
         write(222,*) 'ny       ',ny
         write(222,*) 'nz       ',nz
         write(222,*) 'dx       ',dx
         write(222,*) 'dy       ',dy
         write(222,*) 'dz       ',dz
         write(222,*) 'nsp      ',1
         write(222,*) 'homo     ',energies(nocc)
         write(222,*) 'lumo     ',energies(nocc+1)
         write(222,*) '$ORBITALS'
      endif

      if (write_wf) then
         if (flg_bin) then
            open(unit=1,file='wf.bin',form='unformatted')
            write(1) 'nx       ',nx
            write(1) 'ny       ',ny
            write(1) 'nz       ',nz
            write(1) 'dx       ',dx
            write(1) 'dy       ',dy
            write(1) 'dz       ',dz
            write(1) 'nsp      ',one
            write(1) 'nstates  ',nstates
            write(1) 'evls     '
            write(1) energies
            write(1) 'orbitals '
         else
            open(unit=1,file='wf.txt')
            write(1,*) 'nx       ',nx
            write(1,*) 'ny       ',ny
            write(1,*) 'nz       ',nz
            write(1,*) 'dx       ',dx
            write(1,*) 'dy       ',dy
            write(1,*) 'dz       ',dz
            write(1,*) 'nsp      ',one
            write(1,*) 'nstates  ',nstates
            write(1,*) 'evls     '
            write(1,*) energies
            write(1,*) 'orbitals '
         endif
      endif

      if (write_target) then
         open(unit=122,file='homo.txt')
         write(122,*) 'nx       ',nx
         write(122,*) 'ny       ',ny
         write(122,*) 'nz       ',nz
         write(122,*) 'dx       ',dx
         write(122,*) 'dy       ',dy
         write(122,*) 'dz       ',dz
         write(122,*) 'nsp       ',1
         write(122,*) 'orb       ',energies(nocc), 1

         open(unit=123,file='lumo.txt')
         write(123,*) 'nx       ',nx
         write(123,*) 'ny       ',ny
         write(123,*) 'nz       ',nz
         write(123,*) 'dx       ',dx
         write(123,*) 'dy       ',dy
         write(123,*) 'dz       ',dz
         write(123,*) 'nsp       ',1
         write(123,*) 'orb       ',energies(nocc+1), 1

         open(unit=124,file='dens.txt')
         write(124,*) 'nx       ',nx
         write(124,*) 'ny       ',ny
         write(124,*) 'nz       ',nz
         write(124,*) 'dx       ',dx
         write(124,*) 'dy       ',dy
         write(124,*) 'dz       ',dz
         write(124,*) 'nsp       ',1
         write(124,*) 'dens       '
      endif

      if (write_orbj) then
         open(unit=125,file='orbj.txt')
         write(125,*) 'nx       ',nx
         write(125,*) 'ny       ',ny
         write(125,*) 'nz       ',nz
         write(125,*) 'dx       ',dx
         write(125,*) 'dy       ',dy
         write(125,*) 'dz       ',dz
         write(125,*) 'isp       ',1
         write(125,*) 'orbj      '
      endif

    end subroutine write_headers

    subroutine get_nstates
        implicit none

        whereis=index(line,'=')+len('=')
        read(line(1+whereis:),*) nstates
    end subroutine get_nstates

    subroutine searchstring(string)
        implicit none
        character(len=25) :: string
        
        gdo : do while (.true.)
            
            read(1,'(A)',iostat=stat) line
            whereis=index(line,trim(string))
            if(whereis.ne.0) then
                exit gdo
            elseif(stat<0) then
                write(*,*) 'error : string not found: ', string
            endif
        enddo gdo

    end subroutine searchstring

    subroutine prep_symbol
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
end program qe2sgw      
