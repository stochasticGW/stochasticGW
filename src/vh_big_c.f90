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
subroutine cvh_sub(cdens,cvh,nin)
  use simple_mpi, only : rank
  use gwm
  implicit none
  integer nin, nxb,nyb,nzb,i
  complex*16 cdens(nin)
  complex*16   cvh(nin)
  complex*16, allocatable :: cdensb(:)
  complex*16, allocatable :: cvhb(:)

  call check(n,nin,' n, nin  ')

  if(scale_vh/=1.and.scale_vh/=2.and.scale_vh/=4) then
     write(6,*)'scale_vh should be 1 or 2 (recc.) or 4 '
     stop
  endif
  nxb = scale_vh*nx
  nyb = scale_vh*ny
  nzb = scale_vh*nz

  ! NEW -- modify to handle non-even grids
  if(mod(nxb,2)/=0.or.mod(nyb,2)/=0.or.mod(nzb,2)/=0) stop ' make nxb,nyb,nzb even '

  call alloc_cdensb_cvhb
  call cdens_to_cdensb(cdens,cdensb)
  call get_cvco(cdensb,nxb,nyb,nzb,vk,cvhb,rank)
  call cvhb_to_cvh(cvhb,cvh)
  call dealloc_cdensb_cvhb
contains
  subroutine alloc_cdensb_cvhb
    implicit none
    integer nb
    nb = n*(scale_vh)**3
    allocate(cdensb(nb), stat=i); if(i/=0) stop ' cdensb '
    allocate(cvhb(nb),   stat=i); if(i/=0) stop ' cvhb '
  end subroutine alloc_cdensb_cvhb

  subroutine dealloc_cdensb_cvhb
    implicit none
    deallocate(cdensb, cvhb)
  end subroutine dealloc_cdensb_cvhb


  subroutine cdens_to_cdensb(cdens,cdensb)
    implicit none
    complex*16 cdens( nx, ny, nz)
    complex*16 cdensb(nxb,nyb,nzb)

    integer mx,my,mz,px,py,pz
    cdensb = 0d0
    
    mx= 1 + nxb/2 - nx/2
    px=     mx+nx-1

    my= 1 + nyb/2 - ny/2
    py=     my+ny-1

    mz= 1 + nzb/2 - nz/2
    pz=     mz+nz-1

    cdensb(mx:px,my:py,mz:pz) = cdens
  end subroutine cdens_to_cdensb

  subroutine cvhb_to_cvh(cvhb,cvh)
    implicit none
    complex*16 cvh( nx, ny, nz)
    complex*16 cvhb(nxb,nyb,nzb)
    integer mx,my,mz,px,py,pz

    mx= 1 + nxb/2 - nx/2
    px=     mx+nx-1

    my= 1 + nyb/2 - ny/2
    py=     my+ny-1

    mz= 1 + nzb/2 - nz/2
    pz=     mz+nz-1

    cvh = cvhb(mx:px,my:py,mz:pz) 

  end subroutine cvhb_to_cvh

end subroutine cvh_sub

