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
subroutine pp_check
  use simple_mpi, only : rank
  use ppm,        only : nrppmx, nrpp, rrpp, vpploc, phipp, lpploc, lpptop, matop, mapat
  implicit none
  integer             :: ma, l, nr, st
  real*8              :: ovrlp
  real*8, allocatable :: dr(:)

  rnk0: if(rank==0) then
     
     ! 
     ! checking that all were read
     ! 

     do ma=1,matop
        if(nrpp(ma).le.0) then
           write(6,*)" Error.  Did not read PP for atom with charge= ",mapat(ma)
           stop
        endif
     enddo

     !
     ! checking nrpp
     !
     do ma=1,matop
        call check_lele(50,nrpp(ma),nrppmx,' 50-nrpp(ma)-nrppmx ')
     enddo

     !
     ! making dr
     !
     maloop : do ma=1, matop

        nr = nrpp(ma)
        allocate(dr(nr), stat=st);    call check0(st,' dr_inpp ')

        call dr_make(rrpp(1,ma),dr,nr)
        lloop: do l=0,lpptop(ma)

           lnonloc : if(l/=lpploc(ma)) then

              ovrlp = sum(rrpp(1:nr,ma)**2*dr(1:nr)*phipp(1:nr,l,ma)**2)
              write(17,*)  ' atomcharge ',mapat(ma),' L ',l,' <phi_L|phi_L> ',ovrlp

              if(abs(ovrlp-1d0)>0.01) then
                 write(6,*)' WARNING: potential problem in pseudopotential: <phi_L|phi_L> isnt 1:'
                 write(6,'(X,2(A,I0),A,F12.6)')' -> atom-charge ',mapat(ma),' and L = ',l,': <phi_L|phi_L> = ',real(ovrlp)
                 write(6,*)
                 call flush(6)
              endif
           end if lnonloc
        end do lloop

        deallocate(dr)
     end do maloop

  end if rnk0
end subroutine pp_check
