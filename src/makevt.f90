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
!option dyn    v(t)= vht+vxct
!option nondyn v(t)= vks+ (vh-vh0) 

subroutine makevt(it)
  use gwm
  implicit none
  integer st, it
  if(.not.allocated(vt)) then
     allocate(vt(n, nsp),stat=st)
     call check0(st,' vt in makevt ')
  endif

  if(flgdyn) then
     call make_vt_dyn
  else
     call make_vt_h(it)
  end if
end subroutine makevt

subroutine make_vt_dyn
  use gwm
  implicit none  
  call dens_to_v_spn(rho,vt,n,nsp)
end subroutine make_vt_dyn

subroutine make_vt_h(it)
  use gwm
  implicit none
  integer st, it, sp
  real*8, allocatable :: vh(:)

  if(.not.allocated(vh0)) then
     if(it>1) stop ' vh0 should have been allocated once it>1 '
     allocate(vh0(n),stat=st)
     call check0(st,' vh0 ')
  endif

  allocate(vh(n),stat=st)
  call check0(st,' vh ')

  call vh_sub(rho_p, vh, n)
  if(it==0) vh0 = vh

  do sp=1,nsp
     vt(:,sp)= vks(:,sp)+ vh-vh0
  enddo

  deallocate(vh)
end subroutine make_vt_h

