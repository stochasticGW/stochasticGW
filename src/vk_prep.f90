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
subroutine prep_vk( nx,ny,nz,dx,dy,dz,vk)
  use simple_mpi, only : rank
  use gwm, only : dim_periodic
  implicit none
  integer nx,ny,nz
  real*8 dx,dy,dz
  real*8 vk(nx,ny,nz)
  select case(dim_periodic)
  case(0)
     call prep_vk_MT( nx,ny,nz,dx,dy,dz,vk)
  case(2)
     call prep_vk_2d( nx,ny,nz,dx,dy,dz,vk)
  case(3)
     call prep_vk_3d( nx,ny,nz,dx,dy,dz,vk)
  case default
     stop ' ERROR: dim_periodic not 0,2,3 '
  end select

  if(rank==0)then
     write(17,*)' minval, maxval vk ',real(minval(vk)), real(maxval(vk))
     call flush(17)
  endif

  end subroutine prep_vk

