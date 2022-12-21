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
subroutine vh_sub(dens,vh,n)
  use simple_mpi, only : rank
  use gwm, only : nx,ny,nz,vk, dv,scale_vh
  implicit none
  integer n, nxb,nyb,nzb,i
  real*8 dens(n)
  real*8   vh(n)
  real*8, allocatable :: densb(:)
  real*8, allocatable ::   vhb(:)

  if(scale_vh/=1.and.scale_vh/=2.and.scale_vh/=4) then
     write(6,*)'scale_vh should be 1 or 2 (recc.) or 4 '
     stop
  endif
  nxb = scale_vh*nx
  nyb = scale_vh*ny
  nzb = scale_vh*nz

  call alloc_densb_vhb
  call dens_to_densb(dens,densb)
  call get_vcon(densb,nxb,nyb,nzb,vk,vhb,rank)
  call vhb_to_vh(vhb,vh)
  call dealloc_densb_vhb
contains
  subroutine alloc_densb_vhb
    implicit none
    integer nb
    nb = n*(scale_vh)**3
    allocate(densb(nb), stat=i); if(i/=0) stop ' densb '
    allocate(vhb(nb),   stat=i); if(i/=0) stop ' vhb '
  end subroutine alloc_densb_vhb

  subroutine dealloc_densb_vhb
    implicit none
    deallocate(densb, vhb)
  end subroutine dealloc_densb_vhb

  subroutine dens_to_densb(dens,densb)
    implicit none
    real*8 dens( nx, ny, nz)
    real*8 densb(nxb,nyb,nzb)

    integer mx,my,mz,px,py,pz
    densb = 0d0
    
    mx= 1 + nxb/2 - nx/2
    px=     mx+nx-1

    my= 1 + nyb/2 - ny/2
    py=     my+ny-1

    mz= 1 + nzb/2 - nz/2
    pz=     mz+nz-1

    densb(mx:px,my:py,mz:pz) = dens
  end subroutine dens_to_densb

  subroutine vhb_to_vh(vhb,vh)
    implicit none
    real*8 vh( nx, ny, nz)
    real*8 vhb(nxb,nyb,nzb)
    integer mx,my,mz,px,py,pz

    mx= 1 + nxb/2 - nx/2
    px=     mx+nx-1

    my= 1 + nyb/2 - ny/2
    py=     my+ny-1

    mz= 1 + nzb/2 - nz/2
    pz=     mz+nz-1

    vh = vhb(mx:px,my:py,mz:pz) 

  end subroutine vhb_to_vh

end subroutine vh_sub

