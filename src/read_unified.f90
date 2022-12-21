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
subroutine read_unified
! This subroutine must be called AFTER read_input_vars!
  use simple_mpi, only : rank, sync_mpi
  use simple_mpi, only : bcast_i, bcast_r8, bcast_scalar_i, bcast_scalar_r8, bcast_scalar_L
  use gwm, only : use_unified_input, read_homo_input, read_lumo_input, read_jidx_input
  use gwm, only : nj, det_dft, filter_cheby, gapped_flg, orb_indx

  implicit none
  integer :: ip,inpstat
  integer, allocatable :: inpsections(:,:)
  character(len=16), dimension(3) :: tagsections
  character(len=1024)  :: line

  ip=19
  if (rank.eq.0) then
!    Avoid using sgwinp.txt in these cases:
!    det_dft=T requires reading from wf.txt/wf.bin
!    orb_indx=-1 requires reading from orb.txt
!    nj>0 without j-index input reads from orbj.txt
     if (det_dft .or. orb_indx.eq.-1 .or. &
         (nj.gt.0 .and. .not.read_jidx_input)) then
        use_unified_input=.false.
     else
        open(ip,file='sgwinp.txt',status='old',iostat=inpstat)
        if (inpstat /=0) then
           write(6,*) ' Input file "sgwinp.txt" not found '
           use_unified_input=.false.
        else
           
           write(6,'(2X,A)') 'Input file "sgwinp.txt" detected: reading '
           write(6,'(2X,A)') 'grids, orbitals, density, and atomic coordinates...'
           use_unified_input=.true.
        endif
     endif
  endif
  call bcast_scalar_L(use_unified_input)

  if (use_unified_input) then

!    Find section headers
     call get_sections_in_unified

!    Read grid section
     call read_grid_from_unified

!    Read orbitals and density
     call read_orbs_from_unified

!    Read coordinates section
     call read_cnt_from_unified

     if (rank.eq.0) then
        deallocate(inpsections)
        close(19)
     endif

  else

!    Make sure homo and lumo were read if using gapped
     if (rank.eq.0) then
        if (filter_cheby .and. gapped_flg .and. .not.read_homo_input) then
           write(6,*) ' "mu" variable missing in INPUT so gapped filtering will be used'
           write(6,*) ' HOMO energy is required in INPUT or in "sgwinp.txt" for gapped filtering!'
           stop ' missing value for "homo"'
        endif
        if (filter_cheby .and. gapped_flg .and. .not.read_lumo_input) then
           write(6,*) ' "mu" variable missing in INPUT so gapped filtering will be used'
           write(6,*) ' LUMO energy is required in INPUT or in "sgwinp.txt" for gapped filtering!'
           stop ' missing value for "lumo"'
        endif
     endif
  endif

contains

  subroutine get_sections_in_unified

     implicit none
     integer :: i,nsections,iline,linelen,readstat,thissection
     integer, allocatable :: taglengths(:)
     logical :: found

     if (rank==0) then

        rewind(ip)

!       Section tags: to add new sections, put tags here and increase
!       size of tagsections array above.
        tagsections=(/&
                     '$GEOMETRY',&
                     '$GRID',&
                     '$ORBITALS'&
                     /)
        nsections=SIZE(tagsections)
        allocate(inpsections(nsections,2),taglengths(nsections))
        inpsections(:,:)=0
        do i=1,nsections
           taglengths(i)=LEN_TRIM(ADJUSTL(tagsections(i)))
        enddo

!       Read file and mark line range for each section
        iline=0
        thissection=0
        do
           read(ip,*,iostat=readstat) line
           if (readstat /=0) exit
           iline=iline+1
           line=ADJUSTL(line)
           linelen=LEN_TRIM(line)

           if (line(1:1).eq.'$') then
              found=.false.
              do i=1,nsections
 
                 if (line(1:taglengths(i)) == TRIM(ADJUSTL(tagsections(i)))) then
                    found=.true.

                    if (inpsections(i,1) /= 0) then
                       write(*,*) 'get_sections_in_unified(): ',&
                               TRIM(ADJUSTL(tagsections(i))),' section duplicated!'
                       stop 'read_unified(): problem reading intput file sgwinp.txt'
                    endif

                    inpsections(i,1)=iline
!                   Terminate previous section and mark current section being read
                    if (thissection.gt.0) inpsections(thissection,2)=iline-1
                    thissection=i
                    exit
                 endif

              enddo

              if (.not.found) then ! unrecognized section
                 write(*,*) 'get_sections_in_unified(): ',&
                            'unrecognized section: ',TRIM(line)
                 stop 'read_unified(): problem reading input file sgwinp.txt'
              endif
           endif
        enddo
        inpsections(thissection,2)=iline ! terminate final section

!       Make sure each section is found
        do i=1,nsections
           if (inpsections(i,1) == 0) then
              write(*,*) 'get_sections_in_unified(): no ',&
                         TRIM(ADJUSTL(tagsections(i))),' section found!'
              stop 'read_unified(): problem reading input file sgwinp.txt'
           endif
        enddo

        deallocate(taglengths)

     endif ! rnk0

  end subroutine get_sections_in_unified

  subroutine read_cnt_from_unified
     use atoms, only : an => atom_name, un=> atom_upper_name
     use atoms, only : cnt, ch_a
     use gwm,   only : na, fu=> flg_units_nuc, angstrom_in_au

     implicit none
     logical :: found
     integer :: readstat,iline,ia,j,st
     character*2 :: ccha
     character(len=8) :: utag

     if (rank==0) then

        rewind(ip)

        select case (fu)
        case (1)
           utag='bohr'
        case (2)
           utag='Angstrom'
        case default
           stop ' Error: flg_units_nuc wrong'
        end select

        write(6,*)
        write(6,*) ' Positions in ',TRIM(ADJUSTL(utag)),', read from: sgwinp.txt '
        write(17,*) ' Positions in ',TRIM(ADJUSTL(utag)),', read from: sgwinp.txt '

!       Make sure there are enough lines to have atoms in this section
        na = inpsections(1,2)-inpsections(1,1)
        if (na.lt.1) then
           write(*,*) 'read_cnt_from_unified(): no atoms in input!'
           stop 'read_unified(): problem reading input file sgwinp.txt'
        endif

        write(6,*) ' Number of atoms   : ',na
        write(17,*)' Number of atoms   : ',na

     endif ! rnk0 (1)

     call bcast_scalar_i(na)
     allocate(cnt(3,na),ch_a(na),stat=st); call check0(st,' cnt_allocate ')

     if (rank==0) then

        iline=0
        ia=0
        do
           iline=iline+1
           if (iline.le.inpsections(1,1)) then
              read(ip,*,iostat=readstat) ccha
           elseif (iline.gt.inpsections(1,2)) then
              exit
           else
              ia=ia+1
              read(ip,*,iostat=readstat) ccha,(cnt(j,ia),j=1,3)

              found = .false.
              atom_type : do j=1,size(an)
                 if (ccha==an(j).or.ccha==un(j)) then
                    if (found) stop ' read-atom-label fits more than 1 atom-types '
                    found=.true.
                    ch_a(ia)=j
                 endif
              end do atom_type
              if (.not.found) then
                 write(6,*)' read-atom-label = ',ccha,' does not fit any element '
                 stop
              end if
           endif
        enddo

        select case(fu)
        case(1)
           cnt = cnt
        case(2)
           cnt = cnt* angstrom_in_au
        case default
           stop ' Error: flg_units_nuc wrong'
        end select

        write(17,*)' Atomic charges in molecule: ',ch_a
        call flush(6)
        call flush(17)

     endif ! rnk0 (2)

     call bcast_i(ch_a, size(ch_a), 0)
     call bcast_r8(cnt, size(cnt), 0)

  end subroutine read_cnt_from_unified

  subroutine read_grid_from_unified
     use gwm, only : nx, ny, nz, dx, dy, dz, n, dv, nsp 
     use gwm, only : homo, lumo

     implicit none
     integer :: readstat, j
     real*8  :: homo_rd, lumo_rd
     character*9  :: ch
     character*16 :: a
     character(36) :: homotag, lumotag

     call sync_mpi

     if (rank==0) then

        a='from sgwinp.txt '

!       Find beginning of $GRID section
        rewind(ip)
        do j=1,inpsections(2,1)
           read(ip,*,iostat=readstat) ch
        enddo

!       Read values
        read(ip,*) ch, nx; if (TRIM(ADJUSTL(ch)) /='nx') then; write(6,*)' ERROR, expect nx, got ',ch,a; stop; endif
        read(ip,*) ch, ny; if (TRIM(ADJUSTL(ch)) /='ny') then; write(6,*)' ERROR, expect ny, got ',ch,a; stop; endif
        read(ip,*) ch, nz; if (TRIM(ADJUSTL(ch)) /='nz') then; write(6,*)' ERROR, expect nz, got ',ch,a; stop; endif
        read(ip,*) ch, dx; if (TRIM(ADJUSTL(ch)) /='dx') then; write(6,*)' ERROR, expect dx, got ',ch,a; stop; endif
        read(ip,*) ch, dy; if (TRIM(ADJUSTL(ch)) /='dy') then; write(6,*)' ERROR, expect dy, got ',ch,a; stop; endif
        read(ip,*) ch, dz; if (TRIM(ADJUSTL(ch)) /='dz') then; write(6,*)' ERROR, expect dz, got ',ch,a; stop; endif
        read(ip,*) ch, nsp; if (TRIM(ADJUSTL(ch)) /='nsp') then; write(6,*)' ERROR, expect nsp, got ',ch,a; stop; endif
        read(ip,*) ch, homo_rd; if (TRIM(ADJUSTL(ch)) /='homo') then; write(6,*)' ERROR, expect homo, got ',ch,a; stop; endif
        read(ip,*) ch, lumo_rd; if (TRIM(ADJUSTL(ch)) /='lumo') then; write(6,*)' ERROR, expect lumo, got ',ch,a; stop; endif

!       Process values
        homotag=''
        lumotag=''
        if (.not.gapped_flg) then
           homotag='(ignored since "mu" is set in INPUT)'
           lumotag='(ignored since "mu" is set in INPUT)'
        else
           if (read_homo_input) then
              homotag='(overwritten by value from INPUT)'
           endif
           if (read_lumo_input) then
              lumotag='(overwritten by value from INPUT)'
           endif
        endif

        n = nx*ny*nz
        dv = dx*dy*dz
        call check_r(dble(n),dble(nx)*dble(ny)*dble(nz),' n-nxyz ')
        write(6,*)
        write(6,*) " ########## Grid parameters as read ",a," ########## "
        write(6,*)
        write(6,'(2X,A,23X,I6)') "nx",nx ; if (nx<1) stop ' ERROR: nx should be > 0 '
        write(6,'(2X,A,23X,I6)') "ny",ny ; if (ny<1) stop ' ERROR: ny should be > 0 '
        write(6,'(2X,A,23X,I6)') "nz",nz ; if (nz<1) stop ' ERROR: nz should be > 0 '
        write(6,'(2X,A,18X,I12)') "n",n
        if (max(mod(nx,2),mod(ny,2),mod(nz,2))/=0) stop ' ERROR: nx, ny, or nz not even ! '
        write(6,'(2X,A,12X,F16.8)') "dx ",dx  ;  if(dx.le.0d0) stop ' ERROR: dx should be > 0 '
        write(6,'(2X,A,12X,F16.8)') "dy ",dy  ;  if(dy.le.0d0) stop ' ERROR: dy should be > 0 '
        write(6,'(2X,A,12X,F16.8)') "dz ",dz  ;  if(dz.le.0d0) stop ' ERROR: dz should be > 0 '
        write(6,'(2X,A,12X,F16.8)') "dv ",dv
        write(6,'(2X,A,22X,I6)')"nsp",nsp ; if (nsp<1 .or. nsp>2) stop ' ERROR: nsp must be 1 or 2 '
        write(6,'(2X,A,4X,F16.8,X,A)') "homo energy",homo_rd,TRIM(ADJUSTL(homotag))
        write(6,'(2X,A,4X,F16.8,X,A)') "lumo energy",lumo_rd,TRIM(ADJUSTL(lumotag))

        if (gapped_flg) then
           if (.not.read_homo_input) then
              homo=homo_rd
           endif
           if (.not.read_lumo_input) then
              lumo=lumo_rd
           endif
           if (homo.gt.lumo) stop ' ERROR: homo energy > lumo energy'
        endif

     endif ! rnk0

     call bcast_scalar_i(nx)
     call bcast_scalar_i(ny)
     call bcast_scalar_i(nz)
     call bcast_scalar_i(n)
     call bcast_scalar_i(nsp)
     call bcast_scalar_r8(dx)
     call bcast_scalar_r8(dy)
     call bcast_scalar_r8(dz)
     call bcast_scalar_r8(dv)
     call bcast_scalar_r8(homo)
     call bcast_scalar_r8(lumo)

  end subroutine read_grid_from_unified

  subroutine read_orbs_from_unified
     use gwm, only : n, nj, nsp, sp0, eorb, orb_kind, dv
     use gwm, only : orbi=>orb_rd, orbj=>rdorbj, dens0, orbj_indx

     implicit none
     logical :: foundiorb,founddens,jerr
     integer :: readstat,ij,j,iorb,isp,iline,st,iidx
     real*8  :: orbnrg,nrm
     integer, allocatable    :: jidx(:), jspn(:)
     real*8, allocatable     :: jnrg(:), ur(:)
     complex*16, allocatable :: uc(:)
     character*9 :: ch

     if (.not.allocated(orbi)) then; allocate(orbi(n),stat=st); if(st/=0) stop ' orbi allocate problem '
     else; call check(size(orbi),n,' size_orbi, n ')
     endif
     if (.not.allocated(orbj)) then; allocate(orbj(n,nj),stat=st); if(st/=0) stop ' orbj allocate problem '
     else; call check(size(orbj,1),n,' size_orbj_1, n '); call check(size(orbj,2),nj,' size_orbj_2, nj ')
     endif
     if (.not.allocated(dens0)) then; allocate(dens0(n,nsp),stat=st); if(st/=0) stop ' dens0 allocate problem '
     else; call check(size(dens0,1),n,' size_dens0_1, n '); call check(size(dens0,2),nsp,' size_dens0_2, nsp ')
     endif

     if (rank==0) then

        allocate(jidx(nj),jspn(nj),jnrg(nj))
        jidx(:)=0

        write(6,*)
        write(6,*) " ############ READING DFT ORBITALS FROM sgwinp.txt ############"
        write(6,*)
        select case(orb_kind)
           case(1); allocate(ur(n),stat=st); call check0(st,' ur ')
           case(2); allocate(uc(n),stat=st); call check0(st,' uc ')
        case default
           write(6,*)' ERROR orb_kind not 1(=real*8) nor 2(=complex*16) but ',orb_kind
           call flush(6)
           stop
        end select

!       Find beginning of $ORBITALS section
        rewind(ip)
        iline=0
        do j=1,inpsections(3,1)
           read(ip,*,iostat=readstat) ch
           iline=iline+1
        enddo

        foundiorb=.false.
        founddens=.false.
        do
!          Read orbital header, followed by data
           read(ip,*,iostat=readstat) ch, iorb, orbnrg, isp
           select case (orb_kind)
              case(1); read(ip,*) ur
              case(2); read(ip,*) uc
           end select
           iline=iline+1

           if (TRIM(ADJUSTL(ch))=='ORB') then ! orbital read

!             Target orbital (orb_indx from INPUT)
              if (iorb.eq.orb_indx) then ! target orbital

                 if (foundiorb) then
                    write(6,*) 'read_orbs_from_unified(): orbital ',iorb,' is duplicated!'
                    stop 'read_unified(): problem reading input file sgwinp.txt'
                 endif
                 foundiorb=.true.
                 eorb=orbnrg
                 sp0=isp
                 select case (orb_kind)
                    case(1); orbi = ur
                    case(2); orbi = dble(uc)
                 end select

              elseif (iorb.lt.1) then

                 write(6,*) 'read_orbs_from_unified(): orbital index (=',iorb,&
                         ') must be > 0!'
                 stop 'read_unified(): problem reading input file sgwinp.txt'

              endif

!             Check for j orbital (note that target orbital could also be a j orbital)
              do ij=1,nj

                 if (iorb.eq.orbj_indx(ij)) then
                    if (iorb.eq.jidx(ij)) then
                       write(6,*) 'read_orbs_from_unified(): orbital ',iorb,' is duplicated!'
                       stop 'read_unified(): problem reading input file sgwinp.txt'
                    endif
                    jidx(ij)=iorb
                    jnrg(ij)=orbnrg
                    jspn(ij)=isp
                    select case (orb_kind)
                       case(1); orbj(:,ij) = ur
                       case(2); orbj(:,ij) = dble(uc)
                    end select
                    exit
                 endif

              enddo
             
           elseif (TRIM(ADJUSTL(ch))=='DENS') then ! density read
              
              if (founddens) then
                 write(6,*) 'read_orbs_from_unified(): only one "DENS" is allowed!'
                 stop 'read_unified(): problem reading input file sgwinp.txt'
              endif
              founddens=.true.
              select case (orb_kind)
                 case(1); dens0(1:n,1) = ur
                 case(2); dens0(1:n,1) = dble(uc)
              end select

           else ! bad input
              write(6,*) 'read_orbs_from_unified(): expected "ORB" or "DENS" but read ',ch
              stop 'read_unified(): problem reading input file sgwinp.txt'
           endif

           if (isp.gt.nsp) then
              write(6,*) 'read_orbs_from_unified(): orbital spin value isp (=',isp,&
                         ') may not exceed nsp (=',nsp,')'
              stop 'read_unified(): problem reading input file sgwinp.txt'
           endif

           iline=iline+1
           if (iline.ge.inpsections(3,2)) exit
        enddo

!       Make sure that the 'i' orbital and all nj 'j' orbitals were read
        jerr=.false.
        if (.not.foundiorb) then
           write(6,*) 'read_orbs_from_unified(): target orbital ( orb_indx = ',orb_indx,') not found!'
           jerr=.true.
        endif
        do ij=1,nj
           if (jidx(ij).eq.0) then
              write(6,*) 'read_orbs_from_unified(): j orbital (',orbj_indx(ij),') not found'
              jerr=.true.
           endif
        enddo
        if (.not.founddens) then
           write(6,*) 'read_orbs_from_unified(): density not found!'
           jerr=.true.
        endif
        if (jerr) stop 'read_unified(): problem reading input file sgwinp.txt'

        select case(orb_kind)
           case(1); deallocate(ur)
           case(2); deallocate(uc)
        end select

        write(17,*)
        nrm=sum(orbi(:)**2)*dv
        write(17,*) "ORB (",orb_indx,") -> sum(psi^2)",nrm
        write(6,'(X,2(A,I0),2(A,F16.8))') &
              " ORBI(",orb_indx,") read: spin = ",sp0,"; energy = ",eorb,'; norm = ',nrm

        if (abs(nrm-1.d0)>1d-7) &
           write(6,'(X,A)') " -> WARNING: norm of ORBI differs from 1"
        do j=1,nj
           nrm=sum(orbj(:,j)**2)*dv
           write(17,*) "ORB (",jidx(j),") -> sum(psi^2)",nrm
           write(6,'(X,2(A,I0),2(A,F16.8))') &
                 " ORBJ(",jidx(j),") read: spin = ",jspn(j),"; energy = ",jnrg(j),'; norm = ',nrm
           if (abs(nrm-1.d0)>1d-7) &
              write(6,'(X,A,I0,A)') " -> WARNING: norm of ORBJ(",jidx(j),") differs from 1"
        enddo
        nrm=sum(dens0)*dv
        write(6,*) " DENS read: Nr. electrons == integral(dens*d3r)= ",nrm
        write(17,*) " Nr. electrons == integral(dens*d3r)= ",nrm
        call flush(6)
        call flush(17)

        deallocate(jidx,jspn,jnrg)

     endif ! rnk0

     call bcast_scalar_i(sp0)
     call bcast_scalar_r8(eorb)
     call bcast_r8(orbi,size(orbi),0)
     call bcast_r8(orbj,size(orbj),0)
     call bcast_r8(dens0,size(dens0),0)

  end subroutine read_orbs_from_unified

end subroutine read_unified
