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
subroutine vr1_prnt(i,it) ! i: 0 no perturb., 1 perturb
  use gwm
  use simple_mpi, only : cs=> color_size
  use simple_mpi, only : rank
  implicit none
  integer i,j,it,st,mt,mb,ip
  real*4, allocatable :: vs(:,:)  ! s for slice
  call checks

  if(it==0) call vr_open(i,'rpl')

  do j=1,nvr/cs
     call ip_make(i, j, ip)
     call mbt_set(j, mb, mt)
     if(mb.le.mt) then
        call vs_alloc
        call vs_set
        call vs_write
        call vs_dealloc
     end if
  enddo
  if(it==nt) call vr_close(i)
contains
  subroutine checks
    implicit none
    call check_lele(0,i,1,' 0i1 ')
    call check0(mod(nvr,cs),' mod_nvr_cs ')
  end subroutine checks

  subroutine vs_alloc
    implicit none
    allocate(vs(mb:mt,nspv), stat=st); call check0(st,' vs ')
    vs = 0d0
  end subroutine vs_alloc

  subroutine vs_set
    implicit none
    dyn: if(flgdyn) then
       call check(nsp,nspv,' nsp,nspv ')
       vs(mb:mt,:) = real(vt(mb:mt,:) - vks(mb:mt,:)) 
    else
       call check(nspv,1,' nspv,one ')
       vs(mb:mt,1) = real(vt(mb:mt,1) - vks(mb:mt,1))
       call check_vs
    end if dyn
  end subroutine vs_set

  subroutine vs_write
    implicit none
    
    if(it==0)  write(ip)nt,mb,mt,nspv
    
    write(ip)it
    write(ip)vs(mb:mt,1:nspv)
    if(it==nt) write(ip)nt
  end subroutine vs_write

  subroutine check_vs
    implicit none
    real*8 diff
    diff = maxval(abs(vs(mb:mt,1)-(vt(mb:mt,nsp)-vks(mb:mt,nsp))))
    if(diff>1d-5) then
       write(6,*)' ERROR: diff in vs',diff
       stop
    endif
  end subroutine check_vs

  subroutine vs_dealloc
    implicit none
    deallocate(vs)
  end subroutine vs_dealloc
end subroutine vr1_prnt

