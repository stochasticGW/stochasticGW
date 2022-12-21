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
subroutine kb_ovrlp
  use simple_mpi, only : rank
  use kb_mod, only     : mapat, racut_a, matop
  use kb_mod, only     : nrpp, lpptop, lpploc
  use kb_mod, only     : rr=>rrpp, vl=>vlpp, pl => phipp
  use kb_mod, only     : pvp, vp
  implicit none
  integer  ma, l, mr, st  
  real*8,  allocatable :: dr(:)
  real*8,  parameter   :: toll_vr=1d-10

  maloop : do ma=1,matop

     !
     ! maximum r point
     !
     mrmake : do mr=nrpp(ma),1,-1
        do l=0,lpptop(ma)
           if(l/=lpploc(ma)) then
              if (abs(vl(mr,l,ma))>toll_vr) exit mrmake
           endif
        enddo
     enddo mrmake
     mr = min(nrpp(ma),mr+4)
     racut_a(ma) = rr(mr,ma)
     if(rank==0) write(17,*)' atomindex(ma), mr, racut ',ma,mr,racut_a(ma)
     if(rank==0.and.racut_a(ma)>4.d0) write(17,*)' WARNING: racut_a(ma) too high '
     if(rank==0.and.racut_a(ma)>4.d0) write( 6,*)' WARNING: racut_a(ma) too high '

     
     ! 
     ! set dr
     ! 
     
     allocate(dr(1:mr), stat=st); call check0(st,' dr ' )
     call dr_make(rr(1,ma),dr,mr)
     if(rank==0) write(17,*)' ma, minmax dr ',ma,minval(dr),maxval(dr)

     
     !
     ! make ovrlp
     !
     do l=0,lpptop(ma)
        if(l/=lpploc(ma)) then

           pvp(l,ma) =   &
           sum( pl(1:mr,l,ma)**2*rr(1:mr,ma)**2*dr(1:mr)* vl(1:mr,l,ma)) ! no 4 pi--y0 integ.

           if(rank==0) write(17,*)' ma, l, min pl, min rr, mindr, min vl '
           if(rank==0) write(17,*)ma,l,minval(pl(1:mr,l,ma)),minval(rr(1:mr,ma)),&
                                       minval(dr(1:mr)),minval(vl(1:mr,l,ma))
        end if
     end do
     
     !
     ! print
     !
     if(rank==0)write(17,*)' atom charge ',mapat(ma),' maximum grid index ',mr,' dist ',rr(mr,ma)
     do l=0,lpptop(ma)
        if(rank==0.and.l/=lpploc(ma)) &
             write(17,*)' l,atom,cutoff,<p|v|p>',l,mapat(ma),racut_a(ma),pvp(l,ma)
     enddo
     if(rank==0) then
        write(17,*)
        call flush(17)
     endif
     !
     ! deallocate
     !
     deallocate(dr)

  end do maloop
end subroutine kb_ovrlp

