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
subroutine set_nrm_trace
  use simple_mpi, only : rank
  use gwm, only : nrm_trace, filter_trace, evl, ntop, eorb, de0
  implicit none
  integer i
  real*8  f2
  real*8, parameter :: toll = 1d-8

  call check(ntop, size(evl), ' ntop sz_evl ')
  call check_e

  f2 = 0d0
  do i=1, ntop
     f2 =  f2 + filter_trace(evl(i),eorb, de0)**2
  enddo

  if(sqrt(f2)<toll) then
     write(6,*)' ERROR: eorb, min. dist. to evl ',eorb,minval(abs(evl)-eorb)
     write(6,*)' ERROR: min. dist. smaller than de0 = ',de0
     stop
  endif

  nrm_trace = 1d0/sqrt(f2)

  if(rank==0) write(17,*)" nrm_trace for filter_trace ",nrm_trace


contains
  subroutine check_e
    implicit none

    rnk0 : if(rank==0) then
       
       orbif: if(eorb<minval(evl)) then
          write(6,*)" WARNING: eorb=",eorb, " < minval(eval) = ",minval(evl)
       elseif(eorb>maxval(evl)) then
          write(6,*)" WARNING: eorb=",eorb, " > maxval(eval) = ",maxval(evl)
       else
          iloop : do i= 1,ntop-1
             if(evl(i).le.eorb.and.eorb.le.evl(i+1)) then
                write(17,*)" eorb & eigenvalues below/above ",eorb,evl(i),evl(i+1)
                exit iloop
             endif
          enddo iloop
       end if orbif
    end if rnk0

  end subroutine check_e
end subroutine set_nrm_trace
