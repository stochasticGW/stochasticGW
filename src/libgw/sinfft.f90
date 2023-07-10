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
subroutine sinfft_1d(g, n) 
  ! g(r)-->sum(sin(kr) g(r)), k=0,dk,..,(nf-1)*dk, and dk= pi/(nf*dr). Not 2pi!
  ! r and k start from 0
  implicit none
  integer n,st,i,m,nb
  real*8 g(0:n-1), d
  complex*16, allocatable :: ca(:)
  
  allocate(ca(0:2*n-1),stat=st); if(st/=0) stop ' ca in sinfft '
  ca(0:n-1)=g
  ca(0) = 0d0 ! note
  ca(n)=0d0
  do i=1,n-1   
     ca(2*n-i) = -ca(i)   !  0 1 2 3 4->zero 5=-3 6=-2 7=-1 
  enddo

  nb=2*n
  m=nint(dlog(dble(nb))/dlog(2d0))
  call check(2**m,nb,' 2^m,2*n ')

  call fftsa(nb,ca,m)
  d = maxval(abs(dble(ca)))
  g = -aimag(ca(0:n-1))/2d0

  if(d>1d-9) then
     write(6,*) ' potential problem in ca sinfft; d= ',d
     write(6,*) ' vs imag part ',maxval(abs(g))
     write(6,*)' ca(0)       ',ca(0)
     write(6,*)' ca(1)       ',ca(1)
     write(6,*)' ca(n-1)     ',ca(n-1)
     write(6,*)' ca(n)       ',ca(n)
     write(6,*)' ca(2*n-1)   ',ca(2*n-1)
     stop
  endif

  deallocate(ca)
end subroutine sinfft_1d

