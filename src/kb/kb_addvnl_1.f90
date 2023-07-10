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
subroutine addvnl_c_kb(p,hp,n,ms,cnst,ia)
  use kb_mod, only : nx,ny,nz,dv
  use kb_mod, only : mapai,ngs,mapkbg
  use kb_mod, only : vp, pvp
  use kb_mod, only : iub, iut
  use kb_mod, only : lpptop, lpploc
  implicit none

  integer    :: n,ms,ia,ma
  integer    :: igg, l, is, st
  integer    :: ib, it
  integer    :: iu, iukeep
  real*8     :: cnst
  complex*16 :: p(n,ms), hp(n,ms)
  complex*16 :: cs
  complex*16, allocatable :: ce(:)

  ma = mapai(ia)

  sloop : do is=1,ms 
     iu=iub(ia)
     lloop : do l= 0,lpptop(ma)
        lif : if(l/=lpploc(ma)) then
           ib=l**2      ! 0,1,4,9
           it=ib+2*l   ! 0,3,8,15
           
           allocate(ce(ib:it), stat=st);  call check0(st,' ce ')
           ce = 0d0 

           iukeep = iu
           do igg=1,ngs(ia)
              ce(ib:it)=&
              ce(ib:it) +  p( mapkbg(igg,ia),is)*vp(iu:iu+it-ib)
              iu = iu+it-ib+1  
           enddo
           iu = iukeep
           
           ce(:) = ce(:) * dv * cnst / pvp(l,ma)           

           do igg=1,ngs(ia)
              hp(mapkbg(igg,ia),is) = &
              hp(mapkbg(igg,ia),is)+ &
              sum(ce(ib:it)*vp(iu:iu+it-ib))
              iu = iu+it-ib+1  !
           enddo
           deallocate(ce)
        end if lif
     end do lloop
     call check(iu-1,iut(ia),' iu-1, iut(ia) ')
  end do sloop

end subroutine addvnl_c_kb

