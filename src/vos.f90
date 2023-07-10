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
subroutine vos(km,vd)
  use gwm, only : nw
  use gwm, only : ft
  use gwm, only : nt
  implicit none
  integer :: km,many, st, it, kl, kh

  integer, parameter :: kf=1000
  complex*16 vd(km,0:nt)
  complex*16, allocatable :: fw(:,:)
  
  call check_le(nt*2+2,  nw,' nt2nw        ')

  !bottom, top. low, high.

  do kl=1,km,kf
     kh=min(kl+kf-1,km)
     many=kh-kl+1

     allocate(fw(0:nw-1,kl:kh), stat=st);  call check0(st,' fw_d ')

     fw = 0d0
     do it=0,nt ! transposing
        fw(it,kl:kh)=vd(kl:kh,it)
     enddo

     !fw(nw-nt:nw-1,:)= fr(-nt:-1,:)  -- zero!!!                  ! expand grid
     fw(0,:)  = fw(0,:)/2d0                                   ! theta(t=0)= 1/2

     fw = conjg(fw);        
     call fftw_multiple(fw,nw,many);     
     fw = conjg(fw)/dble(nw)                                      ! fft with exp(iwt)
     fw(nw/2:nw-1,:) = conjg(fw( nw/2:nw-1,:))                    ! ret  to non ret
     call fftw_multiple(fw,nw,many)                                ! fft with exp(-iwt)

     do it=0,nt
        vd( kl:kh,it) = fw(it,kl:kh)
     enddo
     call check_symm_fw
     !vd(-nt:-1,kl:kh) = fw(nw-nt:nw-1,:)  ! not needed, symm    ! compress grid
     deallocate(fw)
  enddo
contains
  subroutine check_symm_fw
    implicit none
    real*8 diff
    diff = maxval(abs(fw(nw-1:nw-nt:-1,kl:kh)-fw(1:nt,kl:kh)))
    if(diff>1d-7) then; write(6,*)' ERROR: fw_drct non symmetrical, stopping; diff= ',diff; stop
    endif
  end subroutine check_symm_fw
end subroutine vos
