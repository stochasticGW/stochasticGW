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
subroutine prnt(j,a)
  use mpi_lib_ours, only : rank,nodes
  implicit none
  character(*) a
  integer j,ip
  integer, save::  i1=1
  real*8,  save :: tstart, t0, t1
  logical, save :: first=.true.
  character(len=1024) :: filename
  character(len=1024) :: fr
  integer, save :: count_p(20)

  if(rank/=0) return 

  if(j>20.or.j<1) stop ' prnt: 1st arguments should be 1-20 '

  ip = 17

  if(first) then;
     count_p(1:20) = 0
     first = .false.
     call cpu_time(tstart)
     t0 = tstart
  endif

  if(count_p(j)>100.and.j/=1) return

  count_p(j) = count_p(j)+1
  call cpu_time(t1)
  write(ip,*)'stage: ',j,' ',a,real(t1-t0),' acc. ',real(t1-tstart); 
  call flush(ip)

  t0=t1
  i1 = -1
end subroutine prnt
