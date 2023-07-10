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
subroutine make_chrg_for_wf
  use atoms, only : valch_a, charge_a => ch_a
  use gwm
  use simple_mpi, only : rank
  implicit none
  rnel_neutral = sum(valch_a) 
  rnel         = rnel_neutral - chrg_net
  rnel_orig    = rnel
  if(rank==0) then
     write(6 ,*)' Total charge: ',real(rnel)
     write(17,*)' Total charge: ',real(rnel)
     call flush(6)
     call flush(17)
  endif
end subroutine make_chrg_for_wf
