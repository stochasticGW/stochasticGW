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
subroutine kb_makeracut
  use simple_mpi, only : rank
  use kb_mod, only : nrpp, rrpp, phipp, racut_a, matop
  implicit none
  
  if(rank==0)  call make_racut_a

contains

  subroutine make_racut_a
    implicit none
    integer mr, ma
    real*8, parameter :: toll_phipp=1d-7
    maloop : do ma=1,matop

       !
       ! maximum r point
       !
       mrmake : do mr=nrpp(ma),1,-1
          if (maxval(abs(phipp(mr,:,ma)))>toll_phipp) exit mrmake
       enddo mrmake
       mr = min(nrpp(ma),mr+1) 
       racut_a(ma) = rrpp(mr,ma)
       if(rank==0) write(17,*)' atomindex(ma), mr, racut ',ma,mr,racut_a(ma)
       if(rank==0.and.racut_a(ma)>4.d0) write(17,*)' WARNING: racut_a(ma) too high '
       if(rank==0.and.racut_a(ma)>4.d0) write( 6,*)' WARNING: racut_a(ma) too high '
    end do maloop

  end subroutine make_racut_a

end subroutine kb_makeracut
