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
subroutine color_make
  use simple_mpi, only : rank, nodes, bcast_scalar_i, sync_mpi
  use simple_mpi, only : color_size, color_rank, icolor
  use simple_mpi, only : prepare_mpi_colors
  implicit none
  integer nds
  integer ncolors
  nds = max(nodes,1)

  !
  ! checking, reseting color_size
  !
  if(mod(nds,color_size)/=0) then
     if(rank==0) then
        write(6,*)
        write(6,*)' ERROR: color-blk-size= ',color_size,',should divide cores number =',nodes
        call flush(6)
     endif

     call sync_mpi
     stop

  end if

  call check_le(1,color_size,' one vs. color_size ')

  if(color_size>nds) then
     color_size=nds
     if(rank==0) then
        write(6,*)
        write(6,*)' NOTE: color_size was reduced to be =nodes=',nds
        call flush(6)
     endif
  end if
  
  ncolors = nds/color_size
  call check(ncolors*color_size,nds,' ncolors*color_size, nds ')
  
  call prepare_mpi_colors(color_size)
end subroutine color_make
