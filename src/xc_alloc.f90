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
subroutine xc_alloc
  use gwm
  implicit none
  integer st
  allocate(exce(ne), stat=st); call check0(st,' exce ')
  allocate(vxce(ne), stat=st); call check0(st,' vxce ')
  allocate( wge(ne), stat=st); call check0(st,' wge  ')

  if(nj/=0) then
    allocate(excej(ne,nj), stat=st); if(st/=0) stop ' excej'
    allocate(vxcej(ne,nj), stat=st); if(st/=0) stop ' vxcej'
    allocate( wgej(ne,nj), stat=st); if(st/=0) stop '  wgej'
  endif

  exce = 0d0
  vxce = 0d0
  wge  = 0d0

  if(nj/=0) then
    excej = 0d0
    vxcej = 0d0
    wgej = 0d0
  endif

  if(rdexch) then
     exce(:) = exchange_value
     !!! PT corresponding line for excej ???
  endif

end subroutine xc_alloc
  
  
