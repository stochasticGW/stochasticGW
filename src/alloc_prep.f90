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
subroutine alloc_prep
  use gwm
  implicit none
  call alloc_v
  call alloc_e
  call alloc_gw
  if(nj>0) call alloc_j
contains
  subroutine alloc_v 
    implicit none
    integer st
    allocate(vks(n, nsp),    stat=st); if(st/=0) stop ' vks alloc '
    allocate(vxc(n, nsp),    stat=st); if(st/=0) stop ' vxc alloc '
    allocate(vxc0(n, nsp),   stat=st); if(st/=0) stop ' vxc alloc '
  end subroutine alloc_v

  subroutine alloc_e ! actually, ear set later.
    implicit none
    integer i,ie
    allocate(ear(ne),stat=i); if(i/=0) stop ' ear '
    if(ne/=1) stop ' set ne to 1 for now '
  end subroutine alloc_e

  subroutine alloc_gw
    implicit none
    integer st
    allocate(rho(n, nsp),    stat=st);  if(st/=0) stop ' prob alloc rho   '
    allocate(vt(n, nsp),     stat=st);  if(st/=0) stop ' prob alloc vr    '
    allocate(ct(-nt:nt,ne),  stat=st);  if(st/=0) stop ' prob alloc ct    '
    allocate(cveta(ne),      stat=st);  if(st/=0) stop ' prob alloc cveta '
    allocate(cvxi(ne),       stat=st);  if(st/=0) stop ' prob alloc cvxi  '
    allocate(normpt(ns),     stat=st);  if(st/=0) stop ' prob alloc normpt'
    allocate(ekn(n),         stat=st);  if(st/=0) stop ' prob alloc ek    '
    allocate(expnk(n),       stat=st);  if(st/=0) stop ' prob alloc expnk '
    
    call prnt_mem_or_scratch_req
    if(gamflg==1.or.gamflg==2) then
       call check_lele(1,ngam,10000000,' 1-ngam-10000000 ')
       allocate(ft(0:nt,ngam),stat=st);  
       if(st/=0) then; write(6,*)' ERROR: prob alloc ft; nt,ngam ',nt,ngam 
       endif
       allocate(gamvr(ngam,nsp),    stat=st);  if(st/=0) stop ' prob alloc gamvr '
    end if
  end subroutine alloc_gw

  subroutine alloc_j
    implicit none
    integer st
    allocate(ctj(-nt:nt,ne,nj),  stat=st);  if(st/=0) stop ' prob alloc ctj    '
    allocate(cvetaj(ne,nj),      stat=st);  if(st/=0) stop ' prob alloc cvetaj '
    allocate(cvxij(ne,nj),       stat=st);  if(st/=0) stop ' prob alloc cvxij  '
  end subroutine alloc_j

end subroutine alloc_prep

