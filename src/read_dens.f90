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
subroutine open_dens_file(ip)
  use simple_mpi, only : rank, sync_mpi
  implicit none
  integer ip, ir
  if(rank==0) then
     open(ip,file='dens.txt',status='old',iostat=ir)
     call check0(ir,' error in dens.txt opening ')
     rewind(ip)
  end if
  call sync_mpi
end subroutine open_dens_file

subroutine read_dens(ip)
  use simple_mpi, only : rank, nodes, bcast_scalar_i, bcast_scalar_r8, bcast_r8
  use gwm
  implicit none
  integer ip
  character*40 ch
  rnk0b : if(rank==0) then
     read(ip,*)ch
     if(ch.ne.'dens')then
        write(6,*)' ERROR: in dens.txt read ',ch,', not dens'
     endif
     read(ip,*)dens0
     rnel = sum(dens0)*dv


     write( 6,*) " Density was read.  Number of electrons == integral(dens*d3r)= ",real(rnel)
     call flush(6) 
     write(17,*)
     write(17,*) " Number of electrons == integral(dens*d3r)= ",real(rnel)
     write(17,*) " everything read successfully"
     call flush(17) 
  end if rnk0b  
  call bcast_r8(dens0,size(dens0),0)
end subroutine read_dens

subroutine make_chrg_from_dens
  use atoms,      only : valch_a, charge_a => ch_a
  use gwm
  use simple_mpi, only : rank
  implicit none
  rnel         = sum(dens0*dv) 
  rnel_orig    = rnel
  rnel_neutral = sum(valch_a)
  chrg_net     = rnel_neutral-rnel
  if(rank==0) then
     write(6 ,*)' recalculated net charge (from density): ',real(chrg_net)
     write(17,*)' recalculated net charge (from density): ',real(chrg_net)
     if(abs(chrg_net)>1d0) then
        write(6,*)' WARNING: net charge more than +-1 ; is that what you want? Stopping  '
        call flush(6)
        stop
     end if
     write(6,*)
     call flush(6)
  end if
end subroutine make_chrg_from_dens

subroutine close_dens_file(ip)
  use simple_mpi, only : rank
  implicit none
  integer ip, ir
  if(rank==0) then
     close(ip,iostat=ir)
     call check0(ir,' error in dens.txt closing ')
     rewind(ip)
  end if
end subroutine close_dens_file



  
