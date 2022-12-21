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
subroutine color_indx
  use simple_mpi, only : color_size, rank, nodes
  use gwm,        only : ncolors, icolor
  implicit none
  integer nds
  nds=max(nodes,1)
  call check0(mod(nds,color_size),' mod(nds, color_size) ')
  ncolors = nds/color_size
  icolor  = rank/color_size
end subroutine color_indx
