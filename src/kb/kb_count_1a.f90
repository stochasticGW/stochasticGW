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
subroutine kb_count_1a(ia,iu)  ! iu is a super index over: small-grid, m,l,atom
  use kb_mod
  implicit none
  integer       :: ia, iu, ma, l
  integer       :: igg
  integer       :: ix, iy, iz
  real*8        :: r0(3), ri, x, y, z

  r0 = (/ -nx/2d0*dx, -ny/2d0*dy, -nz/2d0*dz /)

  ma = mapai(ia)

  lloop : do l=0,lpptop(ma)
     lif : if(l/=lpploc(ma)) then
        iggloop : do igg=1,ngs(ia)
           !ig = mapkbg(igg,ia)
           !call convert_super(ig,ix,iy,iz)
           x = rgn(1,igg,ia)
           y = rgn(2,igg,ia)
           z = rgn(3,igg,ia)
           ri = max(sqrt(x**2+y**2+z**2),1.d-10)
           
           rcut : if(ri<racut_a(ma)) then
              select case(l)
              case(0)
                 iu=iu+1
              case(1)
                 iu=iu+3
              case(2)
                 iu=iu+5
              case(3)
                 iu=iu+7
              case default
                 stop ' l too high '
              end select
           end if rcut
        end do iggloop
     end if lif
  end do lloop
end subroutine kb_count_1a
