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
subroutine set_ge
  use gwm
  use gwm, only        : ernk=> eta_rank
  use simple_mpi, only : crnk=> color_rank
  use simple_mpi, only : cs  => color_size
  use simple_mpi, only : color_bcast_r8
  use simple_mpi, only : color_gather_r8
  use simple_mpi, only : color_scatter_r8
  use simple_mpi, only : color_allsum_r8
  use simple_mpi, only : sync_mpi
  use simple_mpi, only : rank
  implicit none
  integer             :: i
  integer             :: st
  real*8, allocatable :: ao(:)
  real*8, allocatable :: gb(:)
  real*8, save        :: one=1d0
  real*8, parameter   :: tollao = 1d-8
  ! checks

  if(.not.allocated(ge))  stop " Error: ge not allocated "

  if(ne/=1) stop " ERROR: at present set for ne=1 only "

  if(crnk==ernk) call rand_r(ge(:,1),n, dv) ! dont move this statement after the "if trace"

  if(.not.trace) then
     ge(:,1) = orb_rd(:) * sqrt(dv)
     if(nj>0) gej(:,1,:)=rdorbj(:,:)*sqrt(dv)
     return
  endif

  stop ' fix set_ge for trace when not using strips '
end subroutine set_ge

