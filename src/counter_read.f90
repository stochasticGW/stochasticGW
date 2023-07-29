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
subroutine counter_read
  use gwm,        only : icounter
  use simple_mpi, only : rank , bcast_scalar_i
  implicit none
  logical fex
  if(rank==0) then
     inquire(file='counter.inp',exist=fex)
     if(.not.fex) stop ' ERROR: counter.inp missing '
     open(3,file='counter.inp',status='old')   
     rewind(3)
     icounter = 0
     read(3,*,end=33)icounter
33   continue
     close(3)

     write(6,*)" icounter =    ",icounter
  end if
  
  call bcast_scalar_i(icounter)
end subroutine counter_read
