program qe2sgw
   use :: iso_fortran_env
   implicit none
   character*9   :: ch
   character(len=9) :: output
   integer       :: nx,ny,nz,nstates,nlines,na,nocc
   integer       :: whereis,stat,i,j,Za,one=1
   real*8        :: dx,dy,dz,junk,lx,ly,lz,geom(3)
   real*8        :: norm,factor 
   real*8, allocatable :: energies(:),orb(:),orb3d(:,:,:),dens(:)
   integer,parameter    :: njmax=1024
   integer              :: nj
   integer, allocatable :: orbj(:)
   logical       :: flg_bin, file_flg
   logical       :: write_uni, write_wf, write_orbj, write_target, write_cnt
   character*50  :: outputfile,prefix,qe_loc,s1
   character*3   :: ibnd
   character*2   :: an(118)
   character(len=150) :: line

   allocate(orbj(njmax))
   orbj(:)=0

!  Setting variables      
   inquire(file="qe2sgw.in",exist=file_flg)
   if (.not. file_flg) stop 'please provide a qe2sgw.in file'
   open(unit=130,file="qe2sgw.in")
   read(130,*) ch, outputfile
   read(130,*) ch, prefix
   read(130,'(A,A)') ch, qe_loc
   read(130,*) ch, flg_bin
   read(130,*) ch, output
   read(130,*) ch, nocc
   read(130,*,iostat=stat) ch,orbj(:)
   if (stat == iostat_end) continue 
!  Count number orbitals to read
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

!  Read number of states   
   s1='number of Kohn-Sham states'
   open(unit=1,file=outputfile)
   call searchstring(s1)
   call get_nstates

!  Generate input file for pp.x
   open(unit=2,file='pp_orb.in')
   write(2,*) "&inputpp"
   write(2,*) "prefix='"//trim(prefix)//"'"
   write(2,*) "plot_num=7"
   write(2,*) "kband(1)=1"
   write(2,*) "kband(2)=",nstates
   write(2,*) "kpoint=1"
   write(2,*) "lsign=.true."
   write(2,*) "/"
   write(2,*) "&plot"
   write(2,*) "iflag=3"
   write(2,*) "output_format=6"
   write(2,*) "/"
   close(2)

!  Execute pp.x
   write(line,*) trim(qe_loc)//"pp.x -i pp_orb.in"
   call execute_command_line(line)

!  Read energies from QE output
   factor=0.03674932217565499
   nlines=(nstates-1)/8+1
   allocate(energies(nstates))
  
   s1='End of self-consistent'
   call searchstring(s1)

   do i=1,3
       read(1,*)
   enddo

   read(1,*) energies
   energies=energies*factor !convert from eV to Hartree
   close(1)

!  Read grid information   
   open(unit=1,file='tmp.pp_K001_B001')
   read(1,*)
   read(1,*)
   read(1,*) na
   read(1,*) nx,dx
   read(1,*) ny,junk,dy
   read(1,*) nz,junk,junk,dz

!  Read and write cnt.ini file   
   lx=dble(nx)*dx
   ly=dble(ny)*dy
   lz=dble(nz)*dz

   call prep_symbol

   if (write_uni) then
      open(unit=222,file='sgwinp.txt')
      write(222,*) '$GEOMETRY'
   endif
   if (write_cnt) then
      open(unit=2,file='cnt.ini')
   endif
   do i=1,na
      read(1,*) za,junk,geom
      geom(1)=geom(1)-lx/2
      geom(2)=geom(2)-ly/2
      geom(3)=geom(3)-lz/2
!     Map atomic coordinates into cell if outside due to low input data precision
      if (geom(1)<-lx/2) geom(1)=geom(1)+lx
      if (geom(2)<-ly/2) geom(2)=geom(2)+ly
      if (geom(3)<-lz/2) geom(3)=geom(3)+lz
      if (geom(1)>lx/2) geom(1)=geom(1)-lx
      if (geom(2)>ly/2) geom(2)=geom(2)-ly
      if (geom(3)>lz/2) geom(3)=geom(3)-lz
      if (write_uni) write(222,*) an(za),geom
      if (write_cnt) write(2,*) an(za),geom
   enddo
   if (write_cnt) close(2)
   close(1)

!  Write header of orbital(s)-containing file
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

!  Read and write orbitals and dens
   allocate(orb(nx*ny*nz),orb3d(nx,ny,nz),dens(nx*ny*nz))
   dens =0d0
   do i=1,nstates
       if(i<10) then
           write(ibnd,'(A2,I1)') '00',i
       else if(i<100) then
           write(ibnd,'(A1,I2)') '0',i
       else if(i<1000) then
           write(ibnd,'(I3)') i
       else 
           stop 'too many bands!'
       endif
     
       outputfile=trim('tmp.pp_K001_B'//ibnd)
       write(*,*) outputfile

       open(unit=2,file=outputfile)
       rewind(2)
       do j=1,na+6
          read(2,*)
       enddo

       read(2,*) orb
       close(2)

       orb3d=reshape(orb,(/nx,ny,nz/),order=(/3,2,1/))
       orb=reshape(orb3d,(/nx*ny*nz/))

       do j=1,nx*ny*nz
          if (orb(j)<0) then
             orb(j)=-sqrt(-orb(j))
          else
             orb(j)=sqrt(orb(j))
          endif
       enddo
!!!    PT modified: normalize orbital
       norm=sum(orb**2)*dx*dy*dz
       write(*,*) 'norm of state',i,'is',norm
       norm=1.d0/sqrt(norm)
       orb=orb*norm
!!!    END PT modified

       if (i.le.nocc) dens = dens + 2d0*abs(orb)**2d0

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
   enddo

   if (write_uni) then
      write(222,*) 'DENS     ',0,0.d0,0
      write(222,*) dens
   endif

   if (write_target) then
      write(124,*) dens
   endif

   deallocate(energies)

   write(line,*) "rm -f tmp.pp_K001*"
   call execute_command_line(line)

contains   
    subroutine get_nstates
        implicit none

        whereis=index(line,'=')+len('=')
        read(line(1+whereis:),*) nstates
    end subroutine

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
