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
subroutine prop_nl_a(p,q,n,ms,ia) 
  use gwm,   only : dt
  implicit none
  integer                 ::  n, ms, ia, is
  complex*16              ::  p(n,ms)
  complex*16              ::  q(n,ms)
  do is=1,ms
     call prop1_nl_a(p(1,is),q(1,is),n,ia)
  enddo
end subroutine prop_nl_a

subroutine prop1_nl_a(p,q,n,ia)
  use gwm,    only : dt
  use kb_mod, only : dv
  use kb_mod, only : mapai
  use kb_mod, only : ngs, mapkbg
  use kb_mod, only : vp, pvp
  use kb_mod, only : iub, iut
  use kb_mod, only : lpptop, lpploc
  implicit none

  integer    :: n, ia
  integer    :: ma, igg, l, st
  integer    :: iu, iukeep
  integer    :: ib, it
  real*8, parameter       :: toll_o = 1d-8
  real*8, allocatable     :: ov(:)
  complex*16              ::  p(n)
  complex*16              ::  q(n)
  complex*16, parameter   :: ci=(0d0,1d0)
  complex*16, allocatable :: ce(:) 

  
  ma = mapai(ia)
  iu=iub(ia)
  lloop : do l= 0,lpptop(ma)
     lif : if(l/=lpploc(ma)) then
        iukeep = iu
        ib=l**2     ! 0,1,4,9
        it=ib+2*l   ! 0,3,8,15
        
        ! ov=<vphi|vphi> in 3d.        
        allocate(ov(ib:it), stat=st); call check0(st,' ov-propnl ')
        ov= 0d0
        do igg=1,ngs(ia)
           ov(ib:it)= ov(ib:it) + vp(iu:iu+it-ib)**2*dv
           iu = iu+it-ib+1  
        enddo
        iu = iukeep

        where (ov<toll_o) ov = toll_o
           
        ! <vphi| psi>        
        allocate(ce(ib:it), stat=st); call check0(st,' ov-propnl ')
        ce = 0d0
        do igg=1,ngs(ia)
           ce(ib:it)= ce(ib:it) +  p( mapkbg(igg,ia))*vp(iu:iu+it-ib)
           iu = iu+it-ib+1  
        enddo
        iu = iukeep
        
        ce(:) = ce(:) * dv /ov(:) 
        
        ce(:) = ce(:) * (exp(-ci*dt/2d0*ov(:)/pvp(l,ma))-1d0) 
        
        do igg=1,ngs(ia)
           q(mapkbg(igg,ia)) = &
           q(mapkbg(igg,ia)) + &
           sum(ce(ib:it)*vp(iu:iu+it-ib))

           iu = iu+it-ib+1  !
        enddo

        deallocate(ce, ov)
     end if lif
  end do lloop

  call check(iu-1,iut(ia),' iu-1, iut(ia) ')

end subroutine prop1_nl_a

