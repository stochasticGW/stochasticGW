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
subroutine vloc_tuma_prep
  use simple_mpi, only : rank, bcast_r8
  use kb_mod

  implicit none
  real*8 :: vloc_c

  call vloc_tot_prep
  call vloc_cnst
  vloc_tot = vloc_tot + vloc_c
  call bcast_r8(vloc_tot, ng, 0)

contains

subroutine vloc_cnst

  ! Payne et al., RMP 1992 -- Eq. 2.14.  Payne has a mistake, it is Z/r + vloc
  !
  implicit none

  logical :: expspc
  integer :: ma, ia, st, ir, mr
  real*8  :: e1, z, pi, ratio, r, dr
  real*8, allocatable :: ec(:)

  vloc_c = 0d0
  if(dim_periodic.eq.0) return
  pi = dacos(-1d0)

  allocate(ec(matop), stat=st)
  call check0(st,' ec ')


  do ma=1,matop
     if(rrpp(1,ma).le.1d-10) rrpp(1,ma) = 1d-10 ! corrected 5/1/19

     mr = nrpp(ma)
     z = -vpploc(mr,ma)*rrpp(mr,ma)

     e1 = 0d0

     expspc= .true.
     ratio = rrpp(2,ma)/rrpp(1,ma)
     do ir=1,mr-1 ! last point unimportant
        ! first determine if exp. spaced
        if(abs(rrpp(ir+1,ma)/rrpp(ir,ma)-ratio)>1d-5) expspc=.false.
     enddo

     e1 = 0d0
     do ir=1,mr-1
        r = rrpp(ir,ma)
        if(expspc) then
           dr = r*log(ratio) ! dr = dexp(s) = r ds.  ds = log(ratio)
        else
           if(ir==1) then; dr = rrpp(2,ma)-rrpp(1,ma)
           else;           dr =(rrpp(ir+1,ma)-rrpp(ir-1,ma))/2d0
           end if
        endif
        e1 = e1 + r*r*dr* (z/r + vpploc(ir,ma)) ! note signs
     enddo
     ec(ma) = e1*4d0*pi

     if(rank==0) write(6,*)' For atom type ',ma,' intgrl[(-Z/r-vloc(r))d3r] = ',ec(ma)
  end do

  vloc_c = 0d0
  do ia=1,na
     ma = mapai(ia)
     vloc_c = vloc_c + ec(ma)
  end do

  vloc_c = vloc_c / (ng*dv)
  if(rank==0)write(6,*)' shift due to pseudopotential ',vloc_c

  deallocate(ec)
end subroutine vloc_cnst

end subroutine vloc_tuma_prep

