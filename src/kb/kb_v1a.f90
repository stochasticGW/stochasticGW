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
subroutine kb_vphi
  use kb_mod, only : na
  implicit none
  integer  ia
  do ia=1,na
     call kb_v1a(ia)
  end do
end subroutine kb_vphi

! based on https://en.wikipedia.org/wiki/Table_of_spherical_harmonics

subroutine kb_v1a(ia)  ! iu is a super index.
  use kb_mod, only : mapai, iub, iut
  use kb_mod, only : cnt
  use kb_mod, only : vp, ngs
  use kb_mod, only : racut_a, mapkbg, rgn
  use kb_mod, only : lpptop, lpploc
  use kb_mod, only : vphipp, nrpp, rrpp
  use kb_mod, only : dx, dy, dz, nx, ny, nz
  implicit none
  integer       :: ia, ma
  integer       :: iu, l
  integer       :: ig, igg, ix, iy, iz
  real*8        :: r0(3), ri, x, y, z
  real*8        :: pi, vi, as, ap, ap0, ad, ad0, af, af0, dc(5), fc(7)

  pi = dacos(-1d0)
  as  = 1d0/sqrt(4d0*pi) 
  ap0 = sqrt( 3d0)/sqrt(4d0*pi) 
  ad0 = sqrt(15d0)/2d0/sqrt(pi)
  af0 = 1d0/4d0/sqrt(pi)

  dc(:) = (/ 1d0, 1d0, 1d0/2d0/sqrt(3d0), 1d0, 1d0/2d0 /)
  fc(:) = &
  (/ sqrt(35d0/2d0),2*sqrt(105d0),sqrt(21d0/2d0),sqrt(7d0),&
                                  sqrt(21d0/2d0),sqrt(105d0),sqrt(35d0/2d0)/)

  ma = mapai(ia)
  r0 = (/ -(nx)/2d0*dx, -(ny)/2d0*dy, -(nz)/2d0*dz /)
  vp(iub(ia):iut(ia))=0d0
  
  iu=iub(ia)
  lloop: do l=0,lpptop(ma)
     lif : if(l/=lpploc(ma)) then
        iggloop : do igg=1,ngs(ia)
           !ig = mapkbg(igg,ia)
           !call convert_super(ig,ix,iy,iz)
           
           !x = r0(1)+ (ix-1)*dx - cnt(1,ia)
           !y = r0(2)+ (iy-1)*dy - cnt(2,ia)
           !z = r0(3)+ (iz-1)*dz - cnt(3,ia)     
           x = rgn(1,igg,ia)
           y = rgn(2,igg,ia)
           z = rgn(3,igg,ia)
           ri = max(sqrt(x**2+y**2+z**2),1.d-10)
 
           cutoff: if(ri<racut_a(ma)) then
              call interp_lin(ri, vi, rrpp(1,ma), vphipp(1,l,ma), nrpp(ma))
              lcase : select case(l)
              case(0)
                 vp(iu) = as * vi;                           iu=iu+1
              case(1)
                 ap = ap0/ri * vi
                 vp(iu) = ap*x;                              iu=iu+1
                 vp(iu) = ap*y;                              iu=iu+1
                 vp(iu) = ap*z;                              iu=iu+1
              case(2)
                 ad = ad0/ri**2 *vi
                 vp(iu) = ad*dc(1)*x*y;                      iu=iu+1
                 vp(iu) = ad*dc(2)*y*z ;                     iu=iu+1
                 vp(iu) = ad*dc(3)*(2*z**2-x**2-y**2);       iu=iu+1
                 vp(iu) = ad*dc(4)*x*z;                      iu=iu+1
                 vp(iu) = ad*dc(5)*(x**2-y**2);              iu=iu+1
              case(3)
                 af = af0/ri**3 *vi
                 vp(iu) = af*fc(1)*y*(3*x**2-y**2);          iu=iu+1
                 vp(iu) = af*fc(2)*x*y*z;                    iu=iu+1
                 vp(iu) = af*fc(3)*y*(4*z**2-x**2-y**2);     iu=iu+1
                 vp(iu) = af*fc(4)*z*(2*z**2-3*x**2-3*y**2); iu=iu+1
                 vp(iu) = af*fc(5)*x*(4*z**2-x**2-y**2);     iu=iu+1
                 vp(iu) = af*fc(6)*z*(x**2-y**2);            iu=iu+1
                 vp(iu) = af*fc(7)*x*(x**2-3*y**2);          iu=iu+1
              case default
                 stop ' l too high '
              end select lcase
           end if cutoff
        end do iggloop
     end if lif
  end do lloop
  iu = iu-1
  call check(iu,iut(ia),' iu iut ')
end subroutine kb_v1a

subroutine convert_super(ig,ix,iy,iz)
  use kb_mod, only : nx,ny,nz
  implicit none
  integer ig,ix,iy,iz,ig2,jg

  ix = mod((ig-1),nx)+1
  ig2 = (ig-ix)/nx+1 
  iy = mod(ig2-1,ny)+1
  iz = (ig2-iy)/ny+1

  jg = ix+(iy-1)*nx + (iz-1)*nx*ny
  if(jg/=ig) stop ' jg ig '
end subroutine convert_super
