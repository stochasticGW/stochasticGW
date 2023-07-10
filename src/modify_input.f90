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
subroutine modify_nmctot
  use simple_mpi, only : rank, nodes, bcast_scalar_i, color_size
  use gwm, only : ncolors, nmctot, nmc
  implicit none
  integer  nds
  nds = max(nodes,1)
  call check(ncolors*color_size, nds,' ncolor*color_size, nds ')

  if(rank==0) then
     if(mod(nmctot,ncolors)/=0) then
        write(6,*)
        nmctot = nmctot - mod(nmctot,ncolors)+ncolors
        write(6,*)' NOTE: modifying nmctot to be divisble by(#nodes/color_size)=',ncolors
        write(6,*)' NOTE: modified nmctot = ',nmctot
        write(17,*)' NOTE: modifying nmctot to be divisble by(#nodes/color_size)=',ncolors
        write(17,*)' NOTE: modified nmctot = ',nmctot
     endif
     nmc    = nmctot/ncolors
     if(mod(nmctot,nmc)/=0.or.mod(nmctot,ncolors)/=0) stop ' nmc, nmctot problems '
  endif

  call bcast_scalar_i(nmctot)
  call bcast_scalar_i(nmc)
end subroutine modify_nmctot
