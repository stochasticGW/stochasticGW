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
subroutine pp_convert(fname,znum,valnum)
  use simple_mpi, only : rank
  implicit none
  character(50) :: fname
  character*3   :: b
  integer       :: fnl,valnum,znum

  if(rank==0) then
     fnl = len(trim(fname))
     if(fnl<3) then
        write(6,*)' ERROR: file name too short:',fname
        call flush(6)
        stop
     endif
     b=fname(fnl-2:fnl)
     write(17,*) " PP File Ending: ",b
     call lower_the_case(b)
     select case(b)
     case('upf')
        call read_upf(fname,znum,valnum)
     case default
        stop "only upf format supported!"
     end select

     call flush(17)
  endif
end subroutine pp_convert

subroutine read_upf(fname,znum,valnum)
  use ppm, only : lmx, matop
  use ppm, only : nrpp, nrppmx, lpploc, lpptop
  use ppm, only : rrpp, vpploc, rho_core_a, core_correction
  use ppm, only : dij_diag, phipp, lpp, nproj

  implicit none
  integer           :: i,ii,i1,i2,in,np
  character(50)     :: fname
  character(2)      :: name 
  character(80)     :: str1,str2
  character(500)    :: full,f2,tmp
  integer           :: ir,znum,lmax,lloc,sz,lngth,st,licn,il,ipt,valnum
  integer           :: ma
  integer           :: iread
  logical           :: found_cc, core_c
  real*8,allocatable:: grd(:), rho_cr(:)
  real*8,allocatable:: vloc_rd(:)
  real*8,allocatable:: dij(:,:),phi(:,:),nrmp(:),dr(:)
  real*8            :: test, pi, r_valnum, norm2

  pi = dacos(-1d0)

  open(001,file='PP/'//trim(fname))
  read(001,*) str1,str2

  if(str2(scan(str2,'"')+1:scan(str2,'"',.true.)-1)/="2.0.1") &
                               stop "UPF Version 2.0.1 required!"
  
  headerhead : do
     read(001,*)str1
     if(index(str1,"<PP_HEADER")/=0) then
        exit headerhead
     endif
  end do headerhead

  header: do
     read(001,"(A)")full
     read(full,*)str1
     str2=str1(1:scan(str1,'=')-1)
     call lower_the_case(str2)
     select case(trim(str2))
     case("element")
        i=scan(full,'"')
        call check_le(i+2,len_trim(full),' i+2, len_trim(full) ')
        if(full(i+1:i+1)==' ') then ! " H"
           name=full(i+2:i+2)//' '
        else                        ! Arbitrary
           name=full(i+1:i+2)
        endif
        call elname(name,znum) 
        write(17,*) "element:",name,znum
     case("z_valence")
        write(f2,'(A)')str1(scan(str1,'"')+1:scan(str1,'"',.true.)-1)
        read(f2,*) r_valnum
        valnum = nint(r_valnum) ! integer
        if(abs(valnum-r_valnum)>1d-4)then;
           write(6,*)' error: valnum, r_valnum ',valnum,r_valnum
           stop
        endif
        write(17,*) "valence:",valnum
     case("l_max")
        if(len(str1(scan(str1,'"')+1:scan(str1,'"',.true.)-1))/=1) then
           stop "lmax <1!"
        end if
        call str2int(str1(scan(str1,'"')+1:scan(str1,'"',.true.)-1),1,lmax)
        write(17,*) "L max=",lmax
     case("is_ultrasoft")
        if(str1(scan(str1,'"')+1:scan(str1,'"',.true.)-1)=="T") then
           write(17,*) "Ultrasoft not supported!"
           stop
        end if
     case("is_paw")
        if(str1(scan(str1,'"')+1:scan(str1,'"',.true.)-1)=="T") then
           write(17,*) "PAW not supported!"
           stop
        end if
     case("functional")
        str2=str1(scan(str1,'"')+1:scan(str1,'"',.true.)-1) 
        if(str2/='PW'.and.str2/='SLA-PW') then
           write(6,*) "Error: Functional not supported!"
           write(6,*) "Functional is:",str2
           stop
        end if
     case("l_local")
        if(len(str1(scan(str1,'"')+1:scan(str1,'"',.true.)-1))/=1) then
           stop "lloc <10!"
        end if
        call str2int(str1(scan(str1,'"')+1:scan(str1,'"',.true.)-1),1,lloc)
        !if(lloc/=2) stop "L=2 must be local!"
        write(17,*) "L loc=",lloc
     case("core_correction")
        if(str1(scan(str1,'"')+1:scan(str1,'"',.true.)-1)=="T") then
           core_c=.true.
        else if(str1(scan(str1,'"')+1:scan(str1,'"',.true.)-1)=="F") then
           core_c=.false.
        else
           stop ' core_correction ill defined, stopping '
        end if
     case("number_of_proj")
        sz = len(str1(scan(str1,'"')+1:scan(str1,'"',.true.)-1))
        call str2int(str1(scan(str1,'"')+1:scan(str1,'"',.true.)-1),sz,np)
     end select
     if(index(full,"/>")/=0)  exit header
  end do header
  call flush(17)
  !
  ! verifying that atom is among list of nuclei
  ! 
  call check_le(1,znum,' 1, znum ')
  call get_ma(znum,ma)
  if(ma==0) then
     close(001)
     return
  endif
  search_pp_r : do     
     read(001,"(A)") full
     if(index(full,"PP_R")/=0) exit search_pp_r
  enddo search_pp_r
  
  read(full,*)str1,str2,str2
  if(index(str1,"PP_R")==0) stop "Problem reading PP_R!"

  sz = len(str2(scan(str2,'"')+1:scan(str2,'"',.true.)-1))
  call str2int(str2(scan(str2,'"')+1:scan(str2,'"',.true.)-1),sz,lngth)

  write(17,*) "Number of points:",lngth; call flush(17)

  allocate(      grd(lngth),       stat=st); call check0(st,' grd     ')
  allocate(  vloc_rd(lngth),       stat=st); call check0(st,' vloc_rd ')
  allocate(   rho_cr(lngth),       stat=st); call check0(st,' r_core  ')
  allocate(         dr(lngth),    stat=st); call check0(st,' dr        ')
  allocate(      phi(lngth,np), stat=st); call check0(st,' phi     ')
  allocate(      dij(np,np),    stat=st); call check0(st,' dij     ')
  allocate(       nrmp(np),       stat=st); call check0(st,' nrmp      ')

  read(001,*) grd
  ! make dr
  dr(1)      = grd(2)-grd(1)
  dr(lngth)     = grd(lngth)-grd(lngth-1)
  dr(2:lngth-1) = (grd(3:lngth)-grd(1:lngth-2))/2d0

  core_correction(ma) = core_c
  if(core_c) then
     found_cc = .false.
     search_ppnlcc : do 
        read(001,*)str1
        if(str1=="<PP_NLCC") then
           found_cc = .true.
           exit search_ppnlcc
        endif
     end do search_ppnlcc
     if(.not.found_cc) stop ' core_correction but did not find PP_NLCC '
     read(001,*) rho_cr(:)
     write(17,*)' rho_core_integral ',sum(4d0*pi*grd**2*rho_cr)
     if(any(rho_cr<0d0)) stop ' rho_cr_negative somewhere, problem '
  endif

  searchpplocal : do 
     read(001,*) str1
     if(str1=="<PP_LOCAL") exit searchpplocal
  end do searchpplocal
  
  read(001,*) vloc_rd(:)

  !
  !  Next: nonlocal part,
  !
  do
     read(001,*) str1
     if(str1=="<PP_NONLOCAL>") exit
  end do

  !
  !  now need to read PP_BETA.1 ,  PP_BETA.2 , ... THEN: PP_DIJ ....
  ! 

  channeloop : do ii=1,np
     ! read PP_BETA.  Header, angular momentum, and index
     in = -1 ! index
     il = -1 ! l
     PP_Beta_header_search : do
        read(001,'(A)') full
        if(index(full,"<PP_BETA").ne.0) exit PP_Beta_header_search
     end do PP_Beta_header_search

     write(17,*) "Reading ichannel",ii
     PP_BETA_next : do
        i1 = index(full,"angular_momentum=")
        if(i1/=0) then
           write(f2,'(A)') full(i1+18:)
           write(tmp,'(A)') f2(:index(f2,'"')-1)
           read( tmp,*) il
        endif

        i1 = index(full," index")
        if(i1/=0) then
           write(f2,'(A)') full(i1+8:)
           write(tmp,'(A)') f2(:index(f2,'"')-1)
           read(tmp,*) in
        end if

        if(index(full,">").ne.0) exit PP_BETA_next
        read(001,'(A)') full ! reading the next one
     end do PP_BETA_next

     ! now check that l and the index are correct. 
     call check_lele(0,il,9,' 0-il-9 ')
     call check(ii,in,' ii, in ')
     lpp(in,ma)=il

     ! read w.f.
     read(001,*) phi(:,in)
     !
     ! print w.f. properties
     write(17,*)' atom index=',ma,' index=',in
     norm2 = sum(phi(:,in)**2*dr) !!! need dr
     write(17,*)' <phi(:,in)|phi(:,in> ',norm2
     if(abs(norm2-1d0)>1d-4) then
        write(17,*)' potential problem:  in, l, <phi(:,in)|phi(:,in> ',in,il,norm2
     endif
     ! *grd**2 built within phi.

     ! normalize phi
     nrmp(in)=sqrt(norm2)
     phi(:,in) = phi(:,in)/nrmp(in)

     read(001,'(A)') full
     if(index(full,'</PP_BETA')==0) then; write(6,*)' error , full=',full; stop
     end if
  end do channeloop

  !
  !  now read PP_DIJ
  ! 

  read(001,'(A)') full
  if(index(full,'<PP_DIJ')==0) stop ' error, no PP_DIJ '

  end_dij : if(index(full,'>')==0) then
     loopend : do
        read(001,'(A)') full
        if(index(full,'>')/=0) exit loopend
     end do loopend
  end if end_dij

  read(001,*)dij
  write(17,*)' unnormalized dij(in Rydberg) ',dij

  !
  ! check diagonality, extract diagonal value.  
  !
  call check_le(np, size(dij_diag,1),' np, diag_diag ')

  do i1=1,np
     do i2=1,np
        dij(i1,i2) = nrmp(i1)*nrmp(i2)* dij(i1,i2)  ! normalize
        if(i1==i2) then;
           dij_diag(i1,ma) = dij(i1,i1) / 2d0  ! Rydberg to Hartree
        else
           if(abs(dij(i1,i2))>1d-6) stop ' dij not diagonal '
        end if
     end do
  end do
  write(17,*)' normalized dij_diag ',dij_diag(:np,ma)

  if(lmax>3) then
     write(6,*)' ERROR: the value of lmax read from pp is: ',lmax, ' higher than 3, max now '
     stop
  endif

  lpptop(ma)=min(lmax,lmx)
  lpploc(ma)=lloc
  nproj( ma) = np
  call check_le(lngth,nrppmx,' lngth nrppmx ')

  vpploc(    1:lngth,  ma) = vloc_rd(1:lngth)/2.d0               ! Rydberg to Hartree     
  nrpp(                ma) = lngth
  rrpp(      1:lngth,  ma) = grd(1:lngth)

  rho_core_a(:,        ma) = 0d0
  if(core_correction(ma))  rho_core_a(1:lngth,  ma) = rho_cr(1:lngth)

  phipp( 1:lngth,:,ma) = 0d0
  do in=1,np
     do ir=1,lngth
        if(grd(ir)>1d-5) phipp(ir,in,ma) = phi(ir,in)/grd(ir)
     enddo
  end do

  call flush(17)
  close(001)
  deallocate(grd)
  deallocate(vloc_rd)
  deallocate(dr)
  deallocate(phi)
  deallocate(dij)
  deallocate(nrmp)

end subroutine read_upf

subroutine elname(ccha_inp,ic)
  use ppm, only : an => atom_name, un=> atom_upper_name
  implicit none
  character(2) :: ccha, ccha_inp
  integer      :: ic
  integer      :: j
  logical      :: found
  
  ccha(1:2)  = ccha_inp(1:2)
  if( ccha=='  ') stop ' ERROR: element with empty name '

  !
  ! convert ' H'  to 'H ' -- not needed, but keep
  !
  if (ccha(1:1)==' ') then
     ccha(1:1) = ccha(2:2)
     ccha(2:2) = ' '
  end if

  found = .false.
  ic = -99
  atom_type : do j=1,size(an)
     if(ccha==an(j).or.ccha==un(j)) then
        if(found)stop ' read-atom is fit to more than 1 atom-types '
        found=.true.
        ic = j
     endif
  end do atom_type
  if(.not.found) then
     write(6,*)' Error: atom from name in PP, ',ccha,', does not fit any element '
     stop
  end if
  if(ic<1) stop ' Error: didnt find atom ' ! superflous, but always check.
end subroutine elname

subroutine str2int(str,n,num)
  implicit none
  integer, intent(in) :: n
  character(n)     :: str
  integer          :: num
  read(str,*) num
end subroutine str2int

