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
subroutine read_xyz(ip)
  implicit none
  integer  ip
  logical  binrd
  binrd = .false.
  call read_xyz_either(ip, binrd)
end subroutine read_xyz

subroutine read_xyz_wf(ip)
  use gwm, only : binwf
  implicit none
  integer  ip
  call read_xyz_either(ip, binwf)
end subroutine read_xyz_wf

subroutine read_xyz_either(ip, binrd)
  use simple_mpi, only : rank, nodes, bcast_scalar_i, bcast_scalar_r8, sync_mpi
  use gwm,        only : nx, ny, nz, dx, dy, dz, n, dv, nsp
  implicit none
  integer ip, msp, ierr
  integer mx, my, mz
  real*8  ddx, ddy, ddz
  logical binrd
  character*9  :: ch
  character*14 :: a
  character*14 :: aold=""
  
  call sync_mpi

  select case(ip)
  case(111)
     a=' from orb.txt '
  case(440)
     if(binrd)      a=' from wf.bin '
     if(.not.binrd) a=' from wf.txt '
  case(900)
     a=' from dens.txt '
  case(901)
     a=' from vxc.txt ' 
  case default
     write(6,*)' ip should be 111 or 440 or 900, not ',ip; stop
  end select
  
  rnk0: if(rank==0) then 
   if(binrd) then
     write(6,*)' reading ',a
     call flush(6)
     read(ip)ch, mx; if(ch/='nx') then;write(6,*)' ERROR, expect nx,got ',ch,a;stop;endif
     read(ip)ch, my; if(ch/='ny') then;write(6,*)' ERROR, expect ny,got ',ch,a;stop;endif
     read(ip)ch, mz; if(ch/='nz') then;write(6,*)' ERROR, expect nz,got ',ch,a;stop;endif
     read(ip)ch, ddx; if(ch/='dx')then;write(6,*)' ERROR, expect dx,got ',ch,a;stop;endif
     read(ip)ch, ddy; if(ch/='dy')then;write(6,*)' ERROR, expect dy,got ',ch,a;stop;endif
     read(ip)ch, ddz; if(ch/='dz')then;write(6,*)' ERROR, expect dz,got ',ch,a;stop;endif
     read(ip)ch, msp;if(ch/='nsp')then;write(6,*)' ERROR, expect nsp,got ',ch,a;stop;endif
    else
     read(ip,*)ch, mx; if(ch/='nx') then;write(6,*)' ERROR, expect nx,got ',ch,a;stop;endif
     read(ip,*)ch, my; if(ch/='ny') then;write(6,*)' ERROR, expect ny,got ',ch,a;stop;endif
     read(ip,*)ch, mz; if(ch/='nz') then;write(6,*)' ERROR, expect nz,got ',ch,a;stop;endif
     read(ip,*)ch, ddx; if(ch/='dx')then;write(6,*)' ERROR, expect dx,got ',ch,a;stop;endif
     read(ip,*)ch, ddy; if(ch/='dy')then;write(6,*)' ERROR, expect dy,got ',ch,a;stop;endif
     read(ip,*)ch, ddz; if(ch/='dz')then;write(6,*)' ERROR, expect dz,got ',ch,a;stop;endif
     read(ip,*)ch, msp;if(ch/='nsp')then;write(6,*)' ERROR, expect nsp,got ',ch,a;stop;endif
    end if

    call check(msp,nsp,' spin: read vs. from INPUT ')
    if(nx>-1) call check(  nx,mx,' nx not matching: wf.txt or dens.txt or orb.txt or vxc.txt ')
    if(ny>-1) call check(  ny,my,' ny not matching: wf.txt or dens.txt or orb.txt or vxc.txt ')
    if(nz>-1) call check(  nz,mz,' nz not matching: wf.txt or dens.txt or orb.txt or vxc.txt ')
    if(dx>0d0)call check_r(dx,ddx,' dx not matching: wf.txt or dens.txt or orb.txt or vxc.txt ')
    if(dy>0d0)call check_r(dy,ddy,' dy not matching: wf.txt or dens.txt or orb.txt or vxc.txt ')
    if(dz>0d0)call check_r(dz,ddz,' dz not matching: wf.txt or dens.txt or orb.txt or vxc.txt ')
    nx=mx
    ny=my
    nz=mz
    dx=ddx
    dy=ddy
    dz=ddz


    n =nx*ny*nz
    dv=dx*dy*dz
    call check_r(dble(n),dble(nx)*dble(ny)*dble(nz),' n-nxyz ')

    ierr = 0
    if(nx<1.or.ny<1.or.nz<1.or.dx.le.0.d0.or.dy.le.0.d0.or.dz.le.0d0) ierr=1
    if(max(mod(nx,2),mod(ny,2),mod(nz,2))/=0) ierr=1
    
    if(aold/=a.or.ierr==1) then
       if(ierr==1) write(6,*)' ERROR: problem with grids.  Specifically '
       write(6,*)" ############## Grid parameters as read ",a," ############ "
       write(6,*)" nx           ",nx       ;  if(nx<1) stop ' ERROR: nx should be > 0 '
       write(6,*)" ny           ",ny       ;  if(ny<1) stop ' ERROR: ny should be > 0 '
       write(6,*)" nz           ",nz       ;  if(nz<1) stop ' ERROR: nz should be > 0 '
       write(6,*)" n            ",n           
       write(6,'(1x,A,10X,F12.6)') " dx ",dx  ;  if(dx.le.0d0) stop ' ERROR: dx should be > 0 '
       write(6,'(1x,A,10X,F12.6)') " dy ",dy  ;  if(dy.le.0d0) stop ' ERROR: dy should be > 0 '
       write(6,'(1x,A,10X,F12.6)') " dz ",dz  ;  if(dz.le.0d0) stop ' ERROR: dz should be > 0 '
       write(6,'(1x,A,10X,F12.6)') " dv ",dv  
       
       if(max(mod(nx,2),mod(ny,2),mod(nz,2))/=0) stop ' ERROR: nx, ny, or nz not even ! '
    endif

    aold = a

  end if rnk0

  call bcast_scalar_i(nx)
  call bcast_scalar_i(ny)
  call bcast_scalar_i(nz)
  call bcast_scalar_i(n)
  call bcast_scalar_r8(dx)
  call bcast_scalar_r8(dy)
  call bcast_scalar_r8(dz)
  call bcast_scalar_r8(dv)
  call check_r(dble(n),dble(nx)*dble(ny)*dble(nz),' n-nxyz post ')    

end subroutine read_xyz_either
