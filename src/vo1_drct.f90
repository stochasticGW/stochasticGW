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
subroutine vo1_drct(it)
  use gwm
  use simple_mpi, only : cs=> color_size
  use simple_mpi, only : cr=> color_rank
  use simple_mpi, only : color_gather_c8
  use simple_mpi, only : color_bcast_c16, rank
  implicit none
  integer                :: it, j, m, mf, kb, kt, mb, mt
  integer,     parameter :: i2=2
  complex*8, allocatable :: cv(:,:), cvb(:,:)

  if(it==0)call vr_open(i2,'old')  
  call set_bnd
  call alloc_cv

  do j=1,nvr/cs
     call mbt_set(j, mb, mt)!; call check(mt-mb+1,mf,' mtbf ')
     if(mb.le.mt) call vos_read
  enddo
  if(it==nt)call vr_close_rmv(i2)
  call color_gather_c8(cv,cvb,size(cv),0)
  if(cr==0) vo(1:n,1:nspv) = cvb(1:n,1:nspv)
  call color_bcast_c16(vo,n*nspv,0)

  if(.not.flgdyn) vo(:,nsp) = vo(:,1)

  call dealloc_cv
contains

  subroutine set_bnd
    implicit none

    call set_mf(mf)
    m  = mf*(nvr/cs)
    
    kb = 1+m*cr
    kt =   m*(cr+1)
    call check(kt-kb+1,m, ' ktbf ')
    call check_le(kb,  n, ' kb n ')
    call check_lele((cs-1)*m,n,cs*m,' csmn ')
  end subroutine set_bnd

  subroutine alloc_cv
    implicit none
    integer st
    allocate(cv(kb:kt,1:nspv), stat=st)
    cv = 0d0
    call check0(st,  ' cv-vo1 ')

    if(cr==0) then; allocate(cvb(m*cs,1:nspv),stat=st); call check0(st,' cvb-vo1 ')
    else;           allocate(cvb(1,1), stat=st);        call check0(st,' cvb-1 ')
    endif
  end subroutine alloc_cv

  subroutine vos_read
    implicit none
    integer nnt,mmb,mmt,itt,mspv,ip

    if(mb>mt) return
    call ip_make(i2,j,ip)

    if(it==0) then
       read(ip)nnt,mmb,mmt,mspv
       call check(nnt,nt,' nt_beg  ')
       call check(mmb,mb,' mmb     ')
       call check(mmt,mt,' mmt_vos ')
       call check(nspv,mspv,' nmspv_vos ')
    endif

    read(ip)itt; call check(it,itt,' it, itt ')
    read(ip)cv(mb:mt,1:nspv)

    if(it==nt) then
       read(ip)nnt;  call check(nnt, nt,' nt_end ')
    end if
  end subroutine vos_read

  subroutine dealloc_cv
    implicit none
    deallocate(cv,cvb)
  end subroutine dealloc_cv
end subroutine vo1_drct
