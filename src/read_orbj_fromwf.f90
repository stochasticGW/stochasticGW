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
subroutine read_orbj_fromwf
  !
  ! for density and for exchange
  !
  !
  use gwm
  use simple_mpi, only : rank, nodes, sr_mpi_r8, crnk=>color_rank
  use simple_mpi, only : send_receive_r8_in_mpi, color_scatter_r8, bcast_scalar_i
  use simple_mpi, only : bcast_r8
  implicit none
  integer isf,i,j,nds,ir,color_ir,tag,rank0,isp,ic,ip,st
  character*9             :: ch
  real*8                  :: occ
  real*8,     allocatable :: ur(:), urb(:)
  complex*16, allocatable :: uc(:)

  rank0 = 0
  nds = max(nodes, 1)

  ip = 440
  call open_wf_file(ip)
  call read_xyz_wf(ip)
  call alloc_ur_uc
  call read_ntop
  call read_ch
  call check_lele(1,orb_indx,ntop,' one, orb_indx, ntop ')

  if(.not.allocated(rdorbj)) then;
     allocate(rdorbj(n,nj),stat=st)
     if(st/=0) stop ' rdorbj problems '
  endif

  do isf=1,maxval(orbj_indx)
     !
     ! reading 
     !
     rnk0: if(rank==0) then
        select case(nsp)
        case(1)
           if(binwf) then; read(ip,  end=99  )i
           else;           read(ip,*,end=99  )i
           endif
           isp = 1
        case(2)
           if(binwf) then; read(ip, end=99  )i, isp
           else;           read(ip,*,end=99 )i, isp
           endif
           call check_lele(1,isp,2,' 1-isp-2 ')
        case default; stop ' nsp not 1 nor 2 '
        end select
        call check(isf,i,' istate vs. istate_read ')           
        
        call read_ur

        do j=1,nj
           if (isf==orbj_indx(j)) then
              call nrmlz(ur)
              rdorbj(:,j)=ur
              exit
           endif
        enddo

     end if rnk0
  enddo

  call bcast_r8(rdorbj,size(rdorbj),0)

  if(allocated(ur)) deallocate(ur)
  if(allocated(uc)) deallocate(uc)
  call close_wf_file(ip)
  return
99 write(6,*)' problem, isf=',isf,' got state ',i,' but no spin direction '
  stop
contains
  
  subroutine alloc_ur_uc
    implicit none
    allocate(ur(n),stat=st); call check0(st,' ur ')
    select case(orb_det_kind)  ! real or complex
    case(1)
    case(2)       
       if(rank==0) then; allocate(uc(n),stat=st);call check0(st,' uc ')
       endif
    case default;  stop ' orb_det_kind is not 1 neither 2 '
    end select
    ur = 0d0
  end subroutine alloc_ur_uc

  subroutine read_ntop
    use simple_mpi, only : bcast_scalar_i
    implicit none
    if(rank==0) then
       if(binwf) then;  read(ip  )ch,ntop
       else;            read(ip,*)ch,ntop
       endif
       if(ch/='nstates') then; write(6,*)" ERROR reading wf.txt. need nstates,got ",ch; stop
       endif
    endif
    call bcast_scalar_i(ntop)
  end subroutine read_ntop

  subroutine read_ch
    implicit none
    if(rank==0) then
       if(binwf) then;  read(ip);   read(ip);   read(ip  )ch
       else;            read(ip,*); read(ip,*); read(ip,*)ch
       endif
       if(ch/='orbitals') then; write(6,*)" ERROR reading wf.txt. need orbitals,got ",ch; stop
       end if
    endif
  end subroutine read_ch


  subroutine read_ur
    implicit none
    if(binwf) then
       select case(orb_det_kind)
       case(1);     read(ip  )ur(:)
       case(2);     read(ip  )uc(:); call check_real(uc,size(uc)); ur = uc
       end select
    else
       select case(orb_det_kind)
       case(1);  read(ip,*)ur(:)
       case(2);  read(ip,*)uc(:);  call check_real(uc,size(uc));  ur = uc
       end select
    endif
  end subroutine read_ur

  subroutine nrmlz(ur)
    implicit none
    real*8 ur(n)
    real*8 norm2
    norm2 = sum(abs(ur)**2)*dv 
    if(abs(norm2-1d0)>0.01) then; write(6,*)' ERROR: orbital ',isf,' norm2= ',norm2; stop;
    endif
    ur = ur / sqrt(norm2)
  end subroutine nrmlz
   
end subroutine read_orbj_fromwf
