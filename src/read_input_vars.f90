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
subroutine read_input_vars !(inputfname).   
  !
  ! note that units_nuc was read earlier.
  !
  use gwm
  use simple_mpi, only : rank, color_size, nodes, sync_mpi
  use :: iso_fortran_env
  implicit none
  character(lngth_char) :: varn
  character(lngth_char) :: varn_orig
  character(2048)       :: line
  character(50)         :: ppfname_dummy
  character(3)          :: box
  integer       :: i,j,njtmp,jidx(1024),buf(1024)
  integer       :: inpstat, npsi
  integer       :: np
  logical       :: proj           =.false.
  logical       :: rdread_eorb    =.false.
  logical       :: rdmultiplicity =.false.
  logical       :: rdnpsi         =.false.
  logical       :: rdnvr          =.false.
  logical       :: rdchrg_net     =.false.
  logical       :: rdnmctot       =.false.
  logical       :: rdcolor_sz     =.false.
  logical       :: rdns           =.false. 
  logical       :: rdne           =.false.
  logical       :: rdflgdyn       =.false.
  logical       :: rdxc_type      =.false.
  logical       :: rdekcut        =.false.
  logical       :: rddt           =.false.
  logical       :: rdgama         =.false.
  logical       :: rdsm           =.false.
  logical       :: rdtoll         =.false.
  logical       :: rdscale_vh     =.false.
  logical       :: rdoindx        =.false.
  logical       :: rdokind        =.false.
  logical       :: rdodet_kind    =.false.
  logical       :: rdspin         =.false.
  logical       :: rdbinwf        =.false.
  logical       :: rdtrace        =.false.
  logical       :: rddim_periodic =.false.
  logical       :: rdsp0          =.false.
  logical       :: rdgamflg       =.false.
  logical       :: rdseg_frctn    =.false.
  logical       :: rdmxmem_vo     =.false.
  logical       :: rdbox          =.false.
  logical       :: rdfunctional   =.false.
  logical       :: rdmethod_exch  =.false.
  logical       :: rdfuzzy_vk     =.false.
  logical       :: rdfilter       =.false.
  logical       :: rdtp           =.false.
  logical       :: rdnchb         =.false.
  logical       :: rdgapped       =.false.
  logical       :: rdhomo         =.false.
  logical       :: rdlumo         =.false.
  logical       :: rdnj           =.false.
  logical       :: rdjidx         =.false.
  logical       :: rdusegpu       =.false.
  logical       :: rdblock_gam_alg=.false.
  call sync_mpi
  rnk0 : if(rank==0) then
     open(101,file=trim(inputfname),iostat=inpstat)
     if(inpstat /= 0) stop "PROBLEM READING INPUT FILE"
     rewind(101)

     nsp          = 1
     ekcut        = 25d0
     ngam_blk     = 10000
     mxmem_vo     = 6e8   ! 600megabyte
     gamflg       = 2     ! stochastic fragments
     seg_frctn    = 0.003
     npsi         = 10000
     nvr          = 100
     nmctot       = 1000     
     color_size   = 1
     np           = 15
     ne           = 1
     flgdyn       = .false.
     xc_type      = 2 ! lda
     dt           = 0.05d0
     gama         = 0.06d0
     sm           = 1d-4   ! changed
     read_eorb    = .false.

     toll         = 1d-10
     scale_vh     = 2
     chrg_net     = 0d0 ! relevant only if det_tddft. otherwise determined by density (readed)
     orb_indx     = -1      ! -1: use rdorb
     orb_kind     = 1
     orb_det_kind = 1 

     binwf        = .false.
     dim_periodic = 0
     box          = 'pos'
     funct_x      = 0
     funct_c      = 0

     trace        = .false.
     eorb         = -999999d0
     de0          = -999999d0
     sp0          = 1

     method_exch='no_fuzzy_no_0' ! 'fuzzy',  'no_fuzzy_no_0' ,  'v00'
     fuzzy_vk   =.false.

     nj=0
     jidx(:)=0

     filter_cheby = .true.
     havg         =    0d0
     dh           =    1d0 ! will be read or determined
     dhscl        = 1.05d0
     Tp           = 0.01

     usegpu       = .false.
     block_gam_alg= .false.

     write(6,*)' '
     call flush(6)
     
     do 
        read(101,*,iostat=inpstat) varn_orig
        if(inpstat /= 0) exit
        varn = varn_orig
        if(scan(varn(1:2),'#')/=0) varn='#' ! if first or 2nd characrers are #, ignore.
        call lower_the_case(varn)
        select case(varn)
        case('#') ! dont do anything.
        case('multiplicity')
           if(rdmultiplicity)stop 'ERROR: ALREADY READ multiplicity '; rdmultiplicity =.true.
           read(101,*,iostat=inpstat) nsp
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: multiplicity"
        case('binary')
           if(rdbinwf)stop 'ERROR: ALREADY READ binary '; rdbinwf =.true.
           read(101,*,iostat=inpstat) binwf
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: binary"
        case('nxi')
           if(rdnpsi)stop 'ERROR: ALREADY READ nxi '; rdnpsi =.true.
           read(101,*,iostat=inpstat) npsi
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: nxi"
        case('nmctot')
           if(rdnmctot)stop 'ERROR: ALREADY READ nmctot '; rdnmctot =.true.
           read(101,*,iostat=inpstat) nmctot
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: nmctot"
        case('buffer_size')
           if(rdcolor_sz)stop 'ERROR: ALREADY READ buffer_size '; rdcolor_sz =.true.
           read(101,*,iostat=inpstat) color_size
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: buffer_size"
        case('ntddft')
           if(rdns)stop 'ERROR: ALREADY READ number_tddft_states '; rdns =.true.
           read(101,*,iostat=inpstat) np
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: ntddft"
        case('ne')
           read(101,*,iostat=inpstat) ne
           if(rdne)stop 'ERROR: ALREADY READ ne '; rdne =.true.
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: ne"
           if(ne > 1) stop 'ERROR: at present ne should be 1 '
           if(ne .le. 0) stop "ERROR: ne should be > 0"
        case('xc_type')
           read(101,*,iostat=inpstat)xc_type 
           if(rdxc_type)stop 'ERROR: ALREADY READ  '; rdxc_type =.true.
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: xc_type"
        case('exchange_value')
           read(101,*,iostat=inpstat)exchange_value 
           if(rdexch)stop 'ERROR: ALREADY READ  '; rdexch =.true.
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: exchange_value"
        case('flg_dyn')
           read(101,*,iostat=inpstat) flgdyn
           if(rdflgdyn)stop 'ERROR: ALREADY READ  '; rdflgdyn =.true.
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: flg_dyn"
        case('dt')
           if(rddt)stop 'ERROR: ALREADY READ dt '; rddt =.true.
           read(101,*,iostat=inpstat) dt
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: dt"
        case('ekcut')
           if(rdekcut)stop 'ERROR: ALREADY READ ekcut '; rdekcut =.true.
           read(101,*,iostat=inpstat) ekcut
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: ekcut"
        case('gamma')
           if(rdgama)stop 'ERROR: ALREADY READ gamma '; rdgama =.true.
           read(101,*,iostat=inpstat) gama
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: gamma"
        case('sm')
           if(rdsm)stop 'ERROR: ALREADY READ sm '; rdsm =.true.
           read(101,*,iostat=inpstat) sm
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: sm" 
        case('read_eorb')
           if(rdsm)stop 'ERROR: ALREADY READ read_orb '; rdread_eorb =.true.
           read(101,*,iostat=inpstat) read_eorb
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: read_eorb" 
        case('toll')
           if(rdtoll)stop 'ERROR: ALREADY READ toll '; rdtoll =.true.
           read(101,*,iostat=inpstat) toll
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: toll"
        case('scale_vh')
           if(rdscale_vh)stop 'ERROR: ALREADY READ scale_vh '; rdscale_vh =.true.
           read(101,*,iostat=inpstat) scale_vh
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: scale_vh"
        case('orb_indx')
           if(rdoindx)stop 'ERROR: ALREADY READ rdoindx '; rdoindx =.true.
           read(101,*,iostat=inpstat) orb_indx
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: orb_indx"
        case('orb_kind')
           if(rdokind)stop 'ERROR: ALREADY READ orb_kind '; rdokind =.true.
           read(101,*,iostat=inpstat) orb_kind
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: orb_kind"
        case('trace')
           if(rdtrace)stop 'ERROR: ALREADY READ trace '
           rdtrace =.true.
           trace=.true.
           read(101,*,iostat=inpstat) eorb, de0
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT for trace.  Need eorb and de0 "
        case('dim_periodic')
           if(rddim_periodic)stop 'ERROR: ALREADY READ dim_periodic '
           rddim_periodic =.true.
           read(101,*,iostat=inpstat) dim_periodic
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: dim_periodic"
        case('functional')
           if(rdfunctional)stop 'ERROR: ALREADY READ functional '
           rdfunctional =.true.
           read(101,*,iostat=inpstat) funct_x, funct_c
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: functional"
        case('box')
           if(rdbox)stop 'ERROR: ALREADY READ box '
           rdbox =.true.
           read(101,*,iostat=inpstat) box
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: box"
        case('method_exch')
           if(rdmethod_exch)stop 'ERROR: ALREADY READ method_exch '
           rdmethod_exch =.true.
           read(101,*,iostat=inpstat) method_exch
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: method_exch"
        case('fuzzy_vk')
           if(rdfuzzy_vk)stop 'ERROR: ALREADY READ fuzzy_vk '
           rdfuzzy_vk =.true.
           read(101,*,iostat=inpstat) fuzzy_vk
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: fuzzy_vk"
        case('spin_direction')
           if(rdsp0)stop 'ERROR: ALREADY READ spin_direction '
           rdsp0 =.true.
           read(101,*,iostat=inpstat) sp0
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: spin_direction"
        case('segment_fraction')
           if(rdseg_frctn) stop ' ERROR: already read segmenet_fraction '
           rdseg_frctn = .true.
           read(101,*,iostat=inpstat) seg_frctn
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: segment_fraction "
        case('nstripes')
           if(rdnvr) stop ' ERROR: already read nstripes '
           rdnvr = .true.
           read(101,*,iostat=inpstat) nvr
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: nstripes "
        case('max_mem_vot')
           if(rdmxmem_vo)stop 'ERROR: ALREADY READ spin_direction '
           rdmxmem_vo =.true.
           read(101,*,iostat=inpstat) mxmem_vo
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: spin_direction"
        case('units_nuc')
           read(101,*)
        case('scratch_path')
           read(101,*)
        case ('pps')
           atomloop : do 
              read(101,*,iostat=inpstat) ppfname_dummy
              if(scan(ppfname_dummy,'#')/=0) exit atomloop
           enddo atomloop
        case('projection')
           if(rdfilter)stop 'ERROR: ALREADY READ projection'; rdfilter =.true.
           read(101,*,iostat=inpstat) proj
           filter_cheby = .not.proj
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: projection "
        case('tp')
           if(rdtp)stop 'ERROR: ALREADY READ Tp '; rdTp =.true.
           read(101,*,iostat=inpstat) Tp
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: Tp"
        case('mu')
           if(rdmu)stop 'ERROR: ALREADY READ mu '; rdmu =.true.
           read(101,*,iostat=inpstat) mu
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: mu"
        case('nchb')
           if(rdnchb) stop ' ERROR: already read nchb'
           rdnchb= .true.
           read(101,*,iostat=inpstat) nchbmx
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: nchb "
        case('homo')
           if(rdhomo) stop 'ERROR: ALREADY READ homo '; rdhomo =.true.
           read(101,*,iostat=inpstat) homo
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: homo"
        case('lumo')
           if(rdlumo) stop 'ERROR: ALREADY READ lumo '; rdlumo =.true.
           read(101,*,iostat=inpstat) lumo
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: lumo"
        case('nj')
           if(rdnj) stop 'ERROR: ALREADY READ nj '; rdnj =.true.
           read(101,*,iostat=inpstat) nj
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: nj"
        case('orbj_indx')
           if(rdjidx) stop 'ERROR: ALREADY READ orbj_indx '; rdjidx =.true.
           read(101,'(A)',iostat=inpstat) line
           i=1
           do
             read(line,*,iostat=inpstat) buf(1:i)
             if (inpstat /= 0) exit
             i=i+1
           enddo
           i=i-1
           if (i.eq.0) stop "ERROR: VALUES ABSENT FOR VAR: orbj_indx"
           jidx(1:i)=buf(1:i)
        case('usegpu')
           if(rdusegpu)stop 'ERROR: ALREADY READ usegpu '; rdusegpu =.true.
           read(101,*,iostat=inpstat) usegpu
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: usegpu"
        case('block_gam_alg')
           if(rdblock_gam_alg)stop 'ERROR: ALREADY READ block_gam_alg '; rdblock_gam_alg =.true.
           read(101,*,iostat=inpstat) block_gam_alg
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: block_gam_alg"
        case default
           write(6,*)" ERROR: unrecognized varn in input file; varn=",varn_orig
           call flush(6)
           stop
        end select
     end do

     close(101)

     det_tddft = .false.
     if(np==-1) det_tddft =.true.

     det_dft = .false.
     if(np==-1.or.(.not.filter_cheby)) det_dft = .true.

     if (rdmu) gapped_flg=.false. ! read mu to disable gapped filtering
     if (rdhomo) read_homo_input=.true.
     if (rdlumo) read_lumo_input=.true.

     ns_blk = max(-1,np*nsp)
     if(npsi<-1.or.npsi>1e8)             stop "ERROR: npsi not in 1 - 1E8 range"
     if(nsp<1.or.nsp>2)                  stop "ERROR: multiplicity not 1 nor 2"
     if(nmctot<1)                        stop "ERROR: nmctot <1 "
     if(color_size<1)                    stop "ERROR: buffer_size<1 "
     if(dt<0)                            stop "ERROR: dt<0"
     if(gama<0)                          stop "ERROR: gama<0"
     if(sm<1d-7.or.sm>1d-2)              stop "ERROR: sm should be in range 1d-2 to 1d-7 "
     if(toll<0d0)                        stop "ERROR: toll should be positive "
     if(scale_vh.ne.1.and.scale_vh.ne.2) stop 'ERROR: scale_vh SHOULD BE 1 OR 2 '
     if(orb_kind<1.or.orb_kind>2)        stop "ERROR: orb_kind: use  1(real*8) OR 2(complex*16)" 
     if(orb_det_kind<1.or.orb_det_kind>2)stop "ERROR: orb_det_kind: use 1(real*8) or 2(cmplx*16)" 
     if(orb_indx.ne.-1.and.orb_indx<1)   stop "ERROR: orb_indx: use -1 (read from orb.txt) or >0 "
     if(xc_type<1.or.xc_type>2)          stop "ERROR: xc_type at present should be 1 or 2 "   
     if(flgdyn.and.xc_type==1)           stop "ERROR: dynamic xc but xc_type==1,reading xc "     
     if(np/=-1.and.np.le.0)              stop "ERROR: ntddft should be -1 or >0 "
     if(rdsp0.and.(.not.trace))          stop "ERROR: dont give spin_direction if reading orbital" 
     if(rdsp0.and.nsp==1)                stop "ERROR: spin_direction: give only if multiplicity=2"
     if(.not.rdexch.and.trace)           stop "ERROR: if trace needs to read exchange_value "
     if(trace.and.(rdoindx.or.rdokind))  stop "ERROR: trace used, dont use orb_indx,orb_kind"
     if(seg_frctn<1d-7)                  stop "ERROR: filg_seg < 1d-7 "
     if(seg_frctn>0.8d0) then
        write(6,*) " WARNING:make seg_frctn<0.8d0 or 1; rescaling from ",real(seg_frctn)," to ",1
        seg_frctn = 1d0
     end if
     select case(method_exch)
     case('fuzzy','no_fuzzy_no_0','v00');    continue
     case default; stop " ERROR: method_exch isnt 'fuzzy','no_fuzzy_no_0','v00' "
     end select

     ifnpsi: if(npsi.ge.1) then 
        gamflg=2
        ngam_blk = npsi
        nvr = 0
     else
        gamflg=3
        ngam_blk = 0d0
        nvr = nvr
        if(nvr<1)nvr=1
        if(nvr>nvrl) then
           write(6,*)' rescaling number of strips written to scratch disk from',nvrl,' to ',nvr
           nvr = nvrl
        end if
     end if ifnpsi

     if(rdchrg_net.and.(.not.det_dft)) write(6,*) " WARNING:.not.det_dft->dens rd,chrg_net ignrd"
     stoch_x=.true.; if(rdexch.or.det_dft) stoch_x=.false.

     if(color_size>max(nodes,1)) then
        color_size = max(nodes,1); write(6,*)" WARNING: readjusting buffer size to ",color_size
     end if

     if(mod(max(nodes,1),color_size)/=0) then
        write(6,*)" ERROR: buffer_size= ",color_size,",should divide cores number =",nodes; stop
     end if

     select case(dim_periodic)
     case(0)
     case(2,3)
     case default
        write(6,*)" ERROR: dim_periodic should be 0, 2, 3 , not ", dim_periodic; stop
     end select

     if(dim_periodic.gt.0.and.scale_vh==2) then
        scale_vh=1;if(rdscale_vh)write(6,*)"NOTE: periodic, no Martyna-Tuckerman ; scale_vh-> 1"
     end if

     if(mxmem_vo<1e7) stop  ' mxmem_vo < 10Mbyte, too small '
     if(mxmem_vo>5e12) stop ' mxmem_vo >  5Tbyte, probably too big.  '

     select case(box);
     case('pos','"po'); boxpos=.true.
     case('sym','"sy'); boxpos=.false.
     case default; stop ' ERROR: box should be "pos" or "sym" '
     end select

     if(dhscl<1d0)                       stop "ERROR: dhscl SHOULD BE >= 1.0 "
     if(Tp<0d0)                          stop "ERROR: Tp should be positive "
     if(rdmu.and..not.filter_cheby)      stop "ERROR: dont give mu if projection is on"
     if(rdnchb) toll = 1d5 !to allow nchb to be whatever value set
     if (rdhomo .and. rdlumo .and. (lumo .le. homo)) stop " LUMO must be greater than HOMO"

     if (rdjidx) then
        read_jidx_input=.TRUE.
        njtmp=0
        do i=1,SIZE(jidx)
           if (jidx(i).lt.1) exit
           do j=1,i-1
              if (jidx(j).eq.jidx(i)) stop ' ERROR: duplicated entry in orbj_indx'
           enddo
           njtmp=njtmp+1
        enddo
        if (rdnj) then
           if (njtmp.ne.nj) stop ' ERROR: orbj_indx must have nj entries'
        else
           nj=njtmp
        endif
     else
        do i=1,nj
           jidx(i)=i
        enddo
     endif

#if GPU_ENABLED
     if(usegpu) then
       if(color_size.ne.1) stop "ERROR: buffer_size must be 1 for GPU runs "
     endif
#else
     if(usegpu) stop ' build not compiled for GPU'
#endif

#if LIBXC_ENABLED
#else
     if (funct_x.ne.0 .or. funct_c.ne.0) then
        write(*,'(X,2(A,I0,X,I0),A)') 'Functional = (',funct_x,funct_c,&
                ') but LIBXC is not enabled, so only (',0,0,&
                ') (PW-LDA + Slater exch.) is allowed!'
        write(*,'(X,A)') 'Note: (0,0) here corresponds to LIBXC = (1,12)'
        stop ' wrong functional'
     endif
#endif

 end if rnk0
 call bcast_variables
 call bcast_jidx(jidx)
end subroutine read_input_vars

subroutine bcast_variables
  use gwm
  use simple_mpi, only : bcast_scalar_i,bcast_scalar_r8, rank, nodes, sync_mpi
  use simple_mpi, only : bcast_scalar_L,bcast_scalar_char,bcast_char30
  use simple_mpi, only : color_size
  implicit none
  integer i, ir

  !integers
  call bcast_scalar_i(nsp)
  call bcast_scalar_i(color_size)
  call bcast_scalar_i(ngam_blk)
  call bcast_scalar_i(nvr)
  call bcast_scalar_i(nmctot)
  call bcast_scalar_i(ns_blk)
  call bcast_scalar_i(ne)
  call bcast_scalar_i(xc_type)
  call bcast_scalar_i(dim_periodic)
  call bcast_scalar_i(scale_vh)
  call bcast_scalar_i(orb_kind)
  call bcast_scalar_i(orb_det_kind)
  call bcast_scalar_i(orb_indx)
  call bcast_scalar_i(sp0)
  call bcast_scalar_i(gamflg)
  call bcast_scalar_i(nj)
  call bcast_scalar_i(funct_x)
  call bcast_scalar_i(funct_c)
  call bcast_scalar_i(nchbmx)

  !real*8
  call bcast_scalar_r8(dt)
  call bcast_scalar_r8(ekcut)
  call bcast_scalar_r8(gama)
  call bcast_scalar_r8(sm)   
  call bcast_scalar_r8(toll) 
  call bcast_scalar_r8(chrg_net)
  call bcast_scalar_r8(eorb)
  call bcast_scalar_r8(de0)
  call bcast_scalar_r8(mxmem_vo)
  call bcast_scalar_r8(seg_frctn)
  call bcast_scalar_r8(havg) 
  call bcast_scalar_r8(dh)  
  call bcast_scalar_r8(dhscl)  
  call bcast_scalar_r8(tp)   
  call bcast_scalar_r8(homo)
  call bcast_scalar_r8(lumo)

  !logical
  call bcast_scalar_L(flgdyn) 
  call bcast_scalar_L(binwf)
  call bcast_scalar_L(det_tddft)
  call bcast_scalar_L(det_dft)
  call bcast_scalar_L(trace)
  call bcast_scalar_L(boxpos)
  call bcast_scalar_L(filter_cheby)
  call bcast_scalar_L(rdmu)
  call bcast_scalar_L(gapped_flg)
  if(rdmu) call bcast_scalar_r8(mu)
  call bcast_scalar_L(read_homo_input)
  call bcast_scalar_L(read_lumo_input)
  call bcast_scalar_L(read_jidx_input)
  call bcast_scalar_L(usegpu)
  call bcast_scalar_L(block_gam_alg)
  call bcast_scalar_L(fuzzy_vk)

  call bcast_method_exch
  call sync_mpi
  
  !  logical then real*8
  call bcast_scalar_L(rdexch)
  call bcast_scalar_L(stoch_x)
  
contains
  subroutine bcast_method_exch ! overly complicated? maybe call bcast_char30(method_exch,0) ?
    implicit none
    integer i
    character(30) char30

    if(rank==0) then
       i=len_trim(method_exch)
       if (i>30) then
          write(6,*)' method_exch names length=',i
          stop ' method_exch length too long (maxlen=30) '
       endif
       char30=method_exch(1:i)
    endif

    call bcast_char30(char30, 0)
    method_exch=char30
  end subroutine bcast_method_exch

end subroutine bcast_variables

subroutine bcast_jidx(jidx)
  use gwm
  use simple_mpi, only : rank, bcast_i
  implicit none
  integer, intent(in) :: jidx(1024)

  if (nj.gt.0) then
     allocate(orbj_indx(nj))
     if (rank==0) orbj_indx(1:nj)=jidx(1:nj)
     call bcast_i(orbj_indx, size(orbj_indx),0)
  endif

end subroutine bcast_jidx
