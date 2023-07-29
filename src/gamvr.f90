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
subroutine gamvrmake
  use gwm
  implicit none
  integer j, sp, st
  real*8, allocatable ::   vtd(:)

  allocate(vtd(n), stat=st); call check0(st,' vtd ')

  dyn: if(flgdyn) then
     do sp=1,nsp
        vtd(:) = vt(:,sp)-vks(:,sp)
        select case(gamflg)
        case(1); do j=1,ngam; gamvr(j,sp) = sum(gam(:,j)*vtd(:))*dv; enddo
        case(2)
           do j=1,ngam
              call gam_dot( vtd, n, gam(1,j), seg, j, dv, gamvr(j,sp))
           enddo
        case default; write(6,*)' gamflg ',gamflg; stop
        end select
     enddo
  else
     vtd(:) = vt(:,1)-vks(:,1)
     select case(gamflg)
     case(1); do j=1,ngam; gamvr(j,1) = sum(gam(:,j)*vtd(:))*dv; enddo
     case(2)
        do j=1,ngam
           call gam_dot( vtd, n, gam(1,j), seg, j, dv, gamvr(j,1))
           if(nsp>1) gamvr(j,nsp) = gamvr(j,1)
        enddo
     case default; write(6,*)' gamflg ',gamflg; stop
     end select
  end if dyn
  deallocate(vtd)

end subroutine gamvrmake
