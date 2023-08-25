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
  use ppm, only : lpptop
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
     case('fhi')
        call read_fhi(fname,znum,valnum)
     case default
        stop "only fhi or upf format supported!"
     end select

     call flush(17)
  endif
end subroutine pp_convert

subroutine read_fhi(fname,znum,valnum)
  use ppm, only : lmx
  use ppm, only : nrpp, nrppmx, lpploc, lpptop
  use ppm, only : rrpp, vlpp, phipp, vpploc
  use ppm, only : rho_core_a, core_correction
  implicit none
  character(50)     :: fname
  integer           :: ir, znum, dum, lmax, npt, lloc
  integer           :: il, ipt,valnum
  integer           :: ma, st
  real*8,allocatable:: rho_cr(:)
  real*4            :: rzn,rva
  real*8            :: r8a, r8b, r8c, pi

  pi = dacos(-1d0)

  open(001,file='PP/'//trim(fname))
  read(001,*)
  read(001,*) rzn,rva
  znum   = nint(rzn)
  valnum = nint(rva)
  call get_ma(znum,ma)
  if(ma==0) then
     close(001)
     return
  endif

  read(001,*) dum,dum,lmax,lloc

  write(17,*)' charge, valence, lmax, lloc ',znum,valnum,lmax,lloc
  call flush(17)

  if(lmax>3) then
     write(6,*)' ERROR: value of lmax read from pp is: ',lmax,' higher than max allowed now, 3'
     stop
  endif

  if(lmax<lloc) then
     write(6,*)' ERROR: lmax, lloc read from pp are ',lmax,lloc
     write(6,*)' & in fhi format we assume lmax.ge.lloc so change fhi file to have lmax=lloc.'
     stop
  endif

  read(001,*) r8a, r8b
  call get_ma(znum,ma)
  !if(r8b==0.d0) then
  if(abs(r8b)<1d-5) then
    core_correction(ma) = .false.
  else
    core_correction(ma) = .true.
    write(6,*) " For atom charge", znum," and internal index ",ma," a core correction is used!"
    call flush(6)
  end if

  do ir = 1,14
     read(001,*)
  end do

  nrpp(ma)=0
  lloop : do il= 0,lmax  
     read(001,*) npt,r8a
     call check_le(npt,nrppmx,' npt nrppmx ')
     
     if(nrpp(ma)==0) then
        nrpp(ma)=npt
     else
        call check(nrpp(ma),npt,' nrpp_ma, npt ')
     end if
     
     do ipt = 1,npt
        read(001,*) dum, r8a, r8b, r8c ; call check(ipt,dum,' ipt, dum ')

        if(il==0) rrpp(ipt,ma)=r8a
        if(il >0) call check_r(rrpp(ipt,ma),r8a,' rr_ipt_ma, r8a ')

        if(il.le.lmx) then
            phipp(ipt,il,ma) = r8b/rrpp(ipt,ma)  ! note, dividing phi by r.
            vlpp( ipt,il,ma) = r8c               ! no Rydberg factor of half 
        end if

        if(il==lloc) then
           vpploc(ipt, ma)  = r8c               ! no Rydberg factor of half
        end if
     end do
  end do lloop
  
  lpptop(ma)=min(lmx,lmax)
  lpploc(ma)=lloc

  call get_ma(znum,ma)
  rho_core_a(:,ma) = 0d0

  if(core_correction(ma)) then
     allocate(rho_cr(npt),stat=st)
     if(st/=0) stop "Problem allocating rho_cr in FHI PP!"

     do ipt = 1,npt
        read(001,*) r8a, r8b
        rho_cr(ipt) = r8b/(4.d0*pi)     !note factor of 1/4pi 
     end do

     if(npt.gt.size(rho_core_a,1)) stop "Problem with rho_core_a in FHI PP!"
     if(rho_cr(npt).gt.1.D-10)     stop "Too hight core charge in FHI PP!"
     rho_core_a(1:npt,ma) = rho_cr(1:npt)
     deallocate(rho_cr)
  end if

  close(001)

  write(17,*)' ma=atom_type, vpploc_min(for atom_type=ma) ',ma,minval(vpploc(1:nrpp(ma),ma))
  call flush(17)

end subroutine read_fhi

subroutine read_upf(fname,znum,valnum)
  use ppm, only : lmx
  use ppm, only : nrpp, nrppmx, lpploc, lpptop
  use ppm, only : rrpp, vpploc, vlpp, phipp, vpploc, rho_core_a, core_correction
  implicit none
  integer           :: i
  character(50)     :: fname
  character(2)      :: name 
  character(80)     :: str1,str2
  character(500)    :: full,f2
  integer           :: ir,znum,lmax,lloc,sz,lngth,st,licn,il,ipt,valnum
  integer           :: ma
  integer           :: iread
  logical           :: found_cc, core_c
  real*8,allocatable:: grd(:), rho_cr(:)
  real*8,allocatable:: phi(:,:)
  real*8,allocatable:: vl_rd(:,:)
  real*8,allocatable:: vloc_rd(:)
  real*8            :: test, pi, r_valnum

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
     if(index(full,"/>")/=0)  exit header
     str2=str1(1:scan(str1,'=')-1)
     call lower_the_case(str2)
     select case(trim(str2))
     case("element")   ! Three possibilies, generically : " H", or "Si","H "
        i=scan(full,'"')
        call check_le(i+2,len_trim(full),' i+2, len_trim(full) ')
        if(full(i+1:i+1)==' ') then ! " H"
           name=full(i+2:i+2)//' '
        else                        ! "Si" or "H "
           name=full(i+1:i+2)
        endif
        call elname(name,znum) 
        write(17,*) "element:",name,znum
     case("z_valence")
!        call str2int(str1(scan(str1,'"')+1:scan(str1,'"',.true.)-1),1,valnum)
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
     end select
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
  allocate(      phi(lngth,0:lmax),stat=st); call check0(st,' phi     ')
  allocate(    vl_rd(lngth,0:lmax),stat=st); call check0(st,' vl_rd   ')
  allocate(  vloc_rd(lngth),       stat=st); call check0(st,' vloc_rd ')
  allocate(   rho_cr(lngth),       stat=st); call check0(st,' r_core  ')
  read(001,*) grd

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

  do 
     read(001,*) str1
     if(str1=="<PP_SEMILOCAL>") exit
  end do

  do il = 0,lmax
     if(il /= lloc) then
        read(001,*) str1,str1,str1,str1,str2 
        call str2int(str2(scan(str2,'"')+1:scan(str2,'"',.true.)-1),1,licn)
        write(17,*) "Reading icnnel",licn
        call check(il,licn,' il, licn ')
        read(001,*) vl_rd(:,licn)
        write(17,*)"First and Last potential point read for channel",licn," are: ",vl_rd(1,licn),&
             vl_rd(lngth,licn)
        read(001,*)
     end if
  end do

  do 
     read(001,*) str1
     if(index(str1,"<PP_PSWFC>")/=0) exit
  end do

  do il = 0,lmax
     search_closer : do
        read(001,'(A)') full
        if(index(full,'>')/=0) exit search_closer
     end do search_closer
     read(001,*) phi(:,il)
     read(001,'(A)') full
     if(index(full,'</PP_CHI')==0) then 
        write(6,*)' Error.  Should have read </PP_CHI but instead read: ',trim(full)
        stop
     end if
  end do
  
  if(lmax>3) then
     write(6,*)' ERROR: the value of lmax read from pp is: ',lmax, ' higher than 3, max now '
     stop
  endif

  lpptop(ma)=min(lmax,lmx)
  lpploc(ma)=lloc

  call check_le(lngth,nrppmx,' lngth nrppmx ')

  vpploc(    1:lngth,  ma) = vloc_rd(1:lngth)/2.d0               ! Rydberg to Hartree     
  nrpp(                ma) = lngth
  rrpp(      1:lngth,  ma) = grd(1:lngth)
  phipp(     1:lngth,:,ma) = 0d0
  vlpp(      1:lngth,:,ma) = 0d0

  rho_core_a(:,        ma) = 0d0
  if(core_correction(ma))  rho_core_a(1:lngth,  ma) = rho_cr(1:lngth)

  !write(6,*) "DEBUG: at =",ma,"sum rho_cr=",sum(rho_cr(1:lngth))

  lloop : do il = 0,min(lmx,lmax)
     if(il/=lloc) then
        phipp(1:lngth,il,ma) =  phi(1:lngth,il)/grd(1:lngth)     ! note dividing by r 
        vlpp( 1:lngth,il,ma) =  vl_rd(1:lngth,il)/2.d0           ! Rydberg to Hartree     
     end if
  end do lloop

  call flush(17)
  close(001)
  deallocate(grd)
  deallocate(phi)
  deallocate(vl_rd)
  deallocate(vloc_rd)

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

