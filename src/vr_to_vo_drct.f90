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
subroutine vr_to_vo_drct
  use gwm
  use simple_mpi, only : cs=> color_size
  use simple_mpi, only : cr=> color_rank, rank, nodes
  implicit none
  integer j, it, mb, mt
  integer, parameter:: i0=0,i1=1,i2=2
  
  real*4,     allocatable  :: vs0(:,:,:)
  real*4,     allocatable  :: vs1(:,:,:)
  complex*16, allocatable  ::  vd(:,:,:)

  integer irnk

  if(gamflg/=3) return
  
  !
  ! status index
  !
  call vr_open(i0,'old')
  call vr_open(i1,'old')
  call vr_open(i2,'rpl')

  do j=1,nvr/cs
     call mbt_set(j, mb, mt)
     mbmtif: if(mb.le.mt) then

        call alloc_vst
        call vs_read_rmv(  j,mb,mt,i0,vs0)
        call vs_read_rmv(  j,mb,mt,i1,vs1)
        call vd_wg(      mb,mt,vs0,vs1,vd)
        call vos(  (mt-mb+1)*nspv,vd)
        call vd_write
        call dealloc_vst
     else
        call vs_rmv(j,i0)
        call vs_rmv(j,i1)
     end if mbmtif
  end do
  call vr_close(i2)

contains
  subroutine alloc_vst
    implicit none
    integer st
    allocate(vs0(mb:mt,1:nspv,0:nt), stat=st); call check0(st,' vs0 ')
    allocate(vs1(mb:mt,1:nspv,0:nt), stat=st); call check0(st,' vs1 ')
    allocate(vd( mb:mt,1:nspv,0:nt), stat=st); call check0(st,' vd  ')
  end subroutine alloc_vst

  subroutine dealloc_vst
    implicit none
    deallocate(vs0,vs1,vd)
  end subroutine dealloc_vst

  subroutine vd_write
    implicit none
    integer ip
    call ip_make(i2,j,ip)
    write(ip)nt,mb,mt,nspv
    do it=0,nt
       write(ip)it
       ! 
       ! an ugly way of turning complex*16 to complex*8. 
       !
       write(ip)cmplx( real(real(vd(mb:mt,1:nspv,it))),real(aimag(vd(mb:mt,1:nspv,it))))
    enddo
    write(ip)nt
  end subroutine vd_write
end subroutine vr_to_vo_drct
