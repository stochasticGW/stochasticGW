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
subroutine pt_stoch
  use simple_mpi, only : rank
  use rand_seed_mod
  use gwm
  implicit none
  integer   :: line_seed_old, chunk, ib, isp, slc, isf
  slc = 1

  if(det_tddft) stop ' ERROR: need stoch tddft in stoc_pt '

  chunk = 1
  if(rank==0) then
     write(17,'(A,I5,A)')' entering stochastic projection for  = ',ns_blk,' states '
     call flush(17)
  end if

  !
  ! make the functions
  !

  do ib=1,ns
     call set_seed_pt(ib)
     call rand_c(pt(1,ib), n, dv)   ! note that 2 spin functions will get same (spatial) pt.
  enddo

end subroutine pt_stoch

