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
subroutine vr_make
  use gwm
  use simple_mpi, only : sync_mpi
  implicit none
  integer i,it
  integer, external :: is_it_power2 
  complex*16, parameter :: ci=(0.d0,1.d0)

  if(gamflg.eq.1.or.gamflg.eq.2) then
     if(ngam>0) ft = 0d0
     if(ngam<1) stop ' ERROR: ngam<1 '
  end if

  do it=0,nt ! whole part move to gpu
     call prop_pert(it)
     call vr1_process(0,it)
     call prntt(0)
  enddo
  
  call read_pt 
  do i=1,ns
     pt(:,i)=exp( -ci*sm*del(:) )*pt(:,i)     ! for both spins
  enddo
  
  do it=0,nt
     call prop_pert(it)
     call vr1_process(1,it)
     call prntt(1)
  enddo

contains
  subroutine prntt(i)
    implicit none
    integer i
    if(i==0.and.mod(it,100)==0) call prnt(1,'DEBUG: W prop 0 after: 0,100,200,... steps')
    if(i==1.and.mod(it,100)==0) call prnt(1,'DEBUG: W prop 1 after: 0,100,200,... steps')
  end subroutine prntt

end subroutine vr_make

