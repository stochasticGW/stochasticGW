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
subroutine kb_vp_allatoms_extended
  use gwm, only : lmax=>lmx
  use kb_mod, only : nproj_m
  use kb_mod, only : lpp_m
  use kb_mod, only : lpptop
  use kb_mod, only : lpploc
  use kb_mod, only : nsuper_ianl
  use kb_mod, only : indx_ianl
  use kb_mod, only : start_ianl
  use kb_mod, only : vp_hamann
  use kb_mod, only : na
  use kb_mod, only : racut_a
  use kb_mod, only : rgn
  use kb_mod, only : ngs
  use kb_mod, only : mapai

  implicit none
  integer ia, ma, lh, st, np, nilm
  integer igg, ix,iy,iz
  integer i, l, m, j, k, in
  integer, allocatable :: at_ilms(:), at_nilm(:)
  real*8  ra, x, y, z, pi
  real*8  as, ap, ap0, dc(5), fc(7)
  real*8  ad0, ad, af0, af
  real*8, allocatable  :: fm(:,:), pb(:)

  pi = dacos(-1d0)
  call const_angular
  call make_super_ianl

  DO ia=1,na

     ma = mapai(ia)
     np = nproj_m(ma)
     nilm=at_nilm(ia)
     j=at_ilms(ia)

     lh = maxval(lpp_m(1:np,ma))
     call check_lele(0,lh,3,        ' 2-lh-3        ')
     call check_lele(0,lpptop(ma),3,' 2-lpptop_ma-3 ')
     call check_lele(2,lmax,3,      ' 2-lmax-3      ')

     allocate(fm(2*lh+1,0:lh), stat=st)
     call check0(st,' fm alloc ')

     allocate(pb(np), stat=st)
     call check0(st,' pb alloc ')
  
     iggloop : do igg=1,ngs(ia)
        x = rgn(1,igg,ia)
        y = rgn(2,igg,ia)
        z = rgn(3,igg,ia)
        ra = max(sqrt(x*x+y*y+z*z),1.d-8)

        ra_le_racut : if(ra<racut_a(ma)) then

           call make_fm_angular

           call interpp(ra, pb(1:np), ma, np)

           super: do in=1,nilm
               k=j+in
               i=indx_ianl(k,2)
               l=indx_ianl(k,3)
               m=indx_ianl(k,4)
               vp_hamann(igg+start_ianl(k)) = pb(i)*fm(m,l)
           enddo super

        end if ra_le_racut

     enddo iggloop
     deallocate(fm)
     deallocate(pb)

  ENDDO ! ia loop

  deallocate(at_nilm,at_ilms)

  contains
  subroutine make_super_ianl
    implicit none
    integer ia, in, l, m, j, count_ianl

!   indices for i-l-m index numbers and starting points
    allocate(at_nilm(na),at_ilms(na))

!   Count the number of ia-n-l-m combinations
    j=0
    do ia=1,na
       ma=mapai(ia)
       np=nproj_m(ma)
       at_ilms(ia)=j
       at_nilm(ia)=0
       do in=1,np
          l=lpp_m(in,ma)
          if (l/=lpploc(ma).and.l.le.lmax) then
             do m=1,2*l+1
                j=j+1
                at_nilm(ia)=at_nilm(ia)+1
             enddo
          endif
       end do
    end do

    nsuper_ianl = j

    allocate(indx_ianl(nsuper_ianl,4),start_ianl(nsuper_ianl))
    
    j=0
    count_ianl=0
    do ia=1,na
       ma=mapai(ia)
       np=nproj_m(ma)
       do in=1,np
          l=lpp_m(in,ma)
          if (l/=lpploc(ma).and.l.le.lmax) then
             do m=1,2*l+1
                j=j+1
                indx_ianl(j,1)=ia
                indx_ianl(j,2)=in
                indx_ianl(j,3)=l
                indx_ianl(j,4)=m
                start_ianl(j)=count_ianl
                count_ianl=count_ianl+ngs(ia)
             enddo
          endif
       enddo
    enddo

    allocate(vp_hamann(count_ianl))

  end subroutine make_super_ianl

  subroutine const_angular
    implicit none
    
    dc(:) = (/ 1d0, 1d0, 1d0/2d0/sqrt(3d0), 1d0, 1d0/2d0 /)
    fc(:) = &
         (/ sqrt(35d0/2d0),2*sqrt(105d0),sqrt(21d0/2d0),sqrt(7d0),&
            sqrt(21d0/2d0),sqrt(105d0),sqrt(35d0/2d0)/)
    as  = 1d0 / sqrt(4d0*pi) 
    ap0 = as  * sqrt(3d0)
    ad0 = sqrt(15d0)/2d0/sqrt(pi)
    af0 =        1d0/4d0/sqrt(pi)
  end subroutine const_angular

  subroutine make_fm_angular
    implicit none
    fm = 0d0
    fm(1,0) = as
    
    if(lh.ge.1) then
       ap = ap0/ra 
       fm( 1,1) = ap* x 
       fm( 2,1) = ap* y 
       fm( 3,1) = ap* z            
    end if
    
    if(lh.ge.2) then
       ad = ad0/ra**2 
       fm( 1,2) = ad*dc(1)*x*y
       fm( 2,2) = ad*dc(2)*y*z 
       fm( 3,2) = ad*dc(3)*(2*z**2-x**2-y**2)
       fm( 4,2) = ad*dc(4)*x*z
       fm( 5,2) = ad*dc(5)*(x**2-y**2)
    endif
    
    if(lh.ge.3) then
       af       = af0/ra**3 
       fm( 1,3) = af*fc(1)*y*(3*x**2-y**2)
       fm( 2,3) = af*fc(2)*x*y*z           
       fm( 3,3) = af*fc(3)*y*(4*z**2-x**2-y**2)
       fm( 4,3) = af*fc(4)*z*(2*z**2-3*x**2-3*y**2) 
       fm( 5,3) = af*fc(5)*x*(4*z**2-x**2-y**2)  
       fm( 6,3) = af*fc(6)*z*(x**2-y**2)      
       fm( 7,3) = af*fc(7)*x*(x**2-3*y**2)      
    end if
    
    if(lh.ge.4) stop ' error: lpp_high > 3, need to put g states '
  end subroutine make_fm_angular

end subroutine kb_vp_allatoms_extended
