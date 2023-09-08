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
subroutine prop_nl_hamann(p,q,n,ms) 
  !
  ! Assumes that VNL is unpolarized.  Could be extended.
  ! 
  use gwm,    only : dt
  use kb_mod, only : dv 
  use kb_mod, only : mapai
  use kb_mod, only : ngs
  use kb_mod, only : mapkbg
  use kb_mod, only : dij_diag
  use kb_mod, only : nsuper_ianl
  use kb_mod, only : indx_ianl
  use kb_mod, only : start_ianl
  use kb_mod, only : vp_hamann

  implicit none
  integer, intent(in) :: n, ms
  complex*16  :: p(n,ms), q(n,ms)
  integer     :: j, i, igg, is, ia, ma
  complex*16  :: ce 
  complex*16, parameter :: ci=(0d0,1d0)
  
  spn: do is=1,ms
    super: do j=1,nsuper_ianl
    
      ia = indx_ianl(j,1)
      i  = indx_ianl(j,2)
      ma = mapai(ia)
            
      ce = 0d0
      do igg=1,ngs(ia)
         ce = ce + vp_hamann(start_ianl(j)+igg) * p(mapkbg(igg,ia),is)
      enddo

      ce = ce * dv * (exp(-ci*dt/2d0*dij_diag(i,ma))-1d0) 
        
      do igg=1,ngs(ia)
         q(mapkbg(igg,ia),is) = q(mapkbg(igg,ia),is) + ce * vp_hamann(start_ianl(j)+igg)
      end do

    enddo super
  enddo spn
end subroutine prop_nl_hamann
