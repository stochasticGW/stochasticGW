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
subroutine nvr_prep
  use gwm
  use simple_mpi, only : cs=>color_size, rank
  implicit none
  integer nvr_old
  real*8  mem_slice

  nspv = 1
  if(flgdyn)nspv=nsp

  if(gamflg/=3) return
  nvr_old = nvr

  mem_slice = 8.d0 * dble(nt+1) * dble(nvr) *dble(nspv)

  if(mem_slice>mxmem_vo) then ! 8: c*8 or  2 r*4 var., vo(nvr,nspv,0:nt), vr0(...), vr1(...)
     nvr=max(1,int(mxmem_vo/(8d0*dble(nt+1)*dble(nspv))))
     if(rank==0) then
       write(6,*) ' WARNING: memory vs. mxmem_vo readjustment!'
       call flush(6)
     end if
  endif
  
  if(mod(nvr,cs)/=0) then
     nvr = nvr+cs-mod(nvr,cs)
     call check0(mod(nvr,cs),' modnvrcs ')
  end if

  if(rank==0) then
     write(6,*)
     write(6,*)' Nspan was ',nvr_old, ' ; Now: ',nvr,' so mod(nspan/buffer_size)=0 & mem.(nvr*nt)moderate '
     call flush(6)
  endif
end subroutine nvr_prep
  
