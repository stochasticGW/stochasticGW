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
 subroutine debug_analyze_h(pa,vks,n)
  use gwm, only : na, nx, ny, nz, ekn
  implicit none
  integer n, st, ia
  real*8     :: vks(n), e_nl, e_nlb, e_k, e_ks
  complex*16 :: pa(n)
  complex*16, allocatable :: pb(:)

  allocate(pb(n),stat=st)
  call check0(st,' pb_alloc ')
  pb = 0d0
  do ia=1,na
     call debug_addvnl_c_kb(pa,pb,n,ia)
  enddo
  e_nl = sum(conjg(pa)*pb)/sum(conjg(pa)*pa)

  ! 
  ! kin
  ! 

  pb = pa
  call fft3d_forward_many(nx, ny, nz, 1, pb)
  pb = pb * ekn
  call fft3d_backward_many(nx, ny, nz, 1, pb)
  e_k = sum(conjg(pa)*pb)/sum(conjg(pa)*pa)

  !
  ! kohn-sham
  !

  pb = vks(:)*pa
  e_ks = sum(conjg(pa)*pb)/sum(conjg(pa)*pa)
  
  pb = 0d0
  call addvnl_c(pa,pb,n,1,1d0)   ! pb=chp+ vnlscaled p1
  e_nlb = sum(conjg(pa)*pb)/sum(conjg(pa)*pa)

  write(17,*)' <p|vnl|p> ',e_nl, e_nlb
  if(abs(e_nl-e_nlb)>1d-6) then
     write(6,*)' ERROR:  <p|vnl|p>  issues ',e_nl, e_nlb
     stop
  endif
  
  write(17,*)' eorb: kin, ks, nl, tot ',real(e_k), real(e_ks), real(e_nlb),real(e_k+e_ks+e_nlb)
  call flush(17)

  deallocate(pb)
end subroutine debug_analyze_h

subroutine debug_addvnl_c_kb(p,hp,n,ia)
  use kb_mod, only : nx,ny,nz,dv
  use kb_mod, only : mapai,ngs,mapkbg
  use kb_mod, only : vp, pvp
  use kb_mod, only : iub, iut
  use kb_mod, only : lpptop, lpploc
  implicit none

  integer    :: n,ia,ma,ms
  integer    :: igg, l, st
  integer    :: ib, it
  integer    :: iu, iukeep
  complex*16 :: p(n), hp(n)
  complex*16 :: cs
  complex*16, allocatable :: ce(:)

  ma = mapai(ia)

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
                ce(ib:it) +  p( mapkbg(igg,ia))*vp(iu:iu+it-ib)
           iu = iu+it-ib+1  
        enddo
        iu = iukeep
        
        ce(:) = ce(:) * dv / pvp(l,ma)           
        
        do igg=1,ngs(ia)
           hp(mapkbg(igg,ia)) = &
           hp(mapkbg(igg,ia))+ &
                sum(ce(ib:it)*vp(iu:iu+it-ib))
           iu = iu+it-ib+1  !
        enddo
        deallocate(ce)
     end if lif
  end do lloop
  call check(iu-1,iut(ia),' iu-1, iut(ia) ')

end subroutine debug_addvnl_c_kb

