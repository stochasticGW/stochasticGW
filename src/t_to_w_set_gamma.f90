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
subroutine set_gamatw
  use simple_mpi, only : rank
  use gwm
  implicit none
  real*8 tmax
  real*8 pi
  integer iw, iwmin, iwmax
  integer st
  if(rank==0) write(17,*) "setgamatw: gama",gama; call flush(17)
  pi = dacos(-1d0)
  tmax  = 3d0/gama
  nt   = tmax/dt
  nw   = 8* 2**ceiling(log(dble(nt))/log(2d0))
  dw   = 2d0*pi/(nw*dt) 
  if(rank==0) then
     write(6,*) " ############ DETERMINED VARIABLES ############ "
     write(6,*) " tmax         =",real(tmax)
     write(6,*) " nw           =",nw
     write(6,*) " dw           =",real(dw)
     write(6,*)
     write(17,*)' gama, tmax, dt, nt, nw, dw '
     write(17,*)  gama, tmax, dt, nt, nw, dw
     call flush(6)
     call flush(17)
  end if

  iwmin = -nw/2
  iwmax =  nw/2-1

  allocate( war(iwmin:iwmax), stat=st);   call check0(st,' prepew ')

  do iw=iwmin,iwmax
     war(iw)=iw*dw
  end do
  
end subroutine set_gamatw

