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
subroutine read_vxc
  use gwm, only : xc_type
  implicit none
  integer ip
  ip = 901
  if(xc_type/=1) stop ' ERROR: read_vxc should be done only if xc_type=1 '
  call open_vxc_file(ip)
  call read_xyz(ip)
  call read_vxc_in(ip)
  call close_vxc_file(ip)
end subroutine read_vxc

subroutine open_vxc_file(ip)
  use simple_mpi, only : rank, sync_mpi
  implicit none
  integer ip, ierr
  if(rank==0) then
     open(ip,file='vxc.txt',status='old',iostat=ierr)
     call check0(ierr,' ierr in vxc.txt opening ')
     rewind(ip)
  end if
  call sync_mpi
end subroutine open_vxc_file

subroutine read_vxc_in(ip)
  use simple_mpi, only : rank, nodes, bcast_scalar_i, bcast_scalar_r8, bcast_r8
  use gwm
  implicit none
  integer ip
  character*40 ch

  if(.not.allocated(vxc0)) stop ' ERROR: vxc0 not allocated in read_vxc_in '

  rnk0b : if(rank==0) then
     read(ip,*)ch
     if(ch.ne.'vxc')then
        write(6,*)' ERROR: in vxc.txt read ',ch,', not vxc'
     endif
     read(ip,*)vxc0
     
     write( 6,*) " vxc was read. "
     write( 6,*) 
     call flush(6) 
     write(17,*) " vxc was read "
     call flush(17) 
  end if rnk0b
  
  call bcast_r8(vxc0,size(vxc0),0)
  if(.not.allocated(vxc)) stop ' ERROR: vxc not allocated '
  vxc = vxc0
end subroutine read_vxc_in

subroutine close_vxc_file(ip)
  use simple_mpi, only : rank
  implicit none
  integer ip, ierr
  if(rank==0) then
     close(ip,iostat=ierr)
     call check0(ierr,' ierr in vxc.txt closing ')
     rewind(ip)
  end if
end subroutine close_vxc_file



  
