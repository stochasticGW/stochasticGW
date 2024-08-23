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
subroutine xc_expect
  use gwm
  use simple_mpi
  implicit none
  integer i
  call check_alloc_vxc                    
  if(.not.det_dft.and..not.rdexch) then
     call exchange_stoch
     call color_reduce_sum_r8(exce, size(exce), ptheader_rank)  ! exchange should be summed
     if(nj/=0) call color_reduce_sum_c16(excej,size(excej),ptheader_rank)
  endif

  call xc_vw
end subroutine xc_expect

! this subroutine is extremely inefficient; but that may be OK
!  To optimize it -- replace the two vh calculations (4 ffts) by a single fft on a big grid
! Also, put it all on the kernel replacing ne*ne transfers by ne+n
! also, reduce ne.
! At present, speed is not an issue.

subroutine check_alloc_vxc
  use gwm, only : vxc
  implicit none
  if(.not.allocated(vxc)) stop ' vxc should have been prepared by now '
end subroutine check_alloc_vxc

subroutine xc_vw
  use gwm
  use simple_mpi, only : rank
  implicit none
  integer  ie, i, is, ij
  real*8 ex
  real*8, allocatable :: p(:), vp(:)
  allocate(p(n), vp(n), stat=i); if(i/=0) stop ' p problems '
  do ie=1,ne
     !wge( ie) = wge(ie)  + 1d0/dv * sum(ge(:,ie)**2)*dv             ! dont average within color
     !vxce(ie) = vxce(ie) + 1d0/dv * sum(ge(:,ie)**2*vxc(:,sp0))*dv  ! dont average within color
      wge( ie) =            1d0/dv * sum(ge(:,ie)**2)*dv             ! dont average within color
      vxce(ie) =            1d0/dv * sum(ge(:,ie)**2*vxc(:,sp0))*dv  ! dont average within color

      if(nj/=0) then
        do ij=1,nj
          wgej( ie,ij) =  1d0/dv * sum(gej(:,ie,ij)**2)*dv ! no need to average within color
          vxcej(ie,ij) =  1d0/dv * sum(gej(:,ie,ij)*ge(:,ie)*vxc(:,sp0))*dv ! no need to average within color
        enddo
      endif
  enddo
  deallocate(p, vp)
end subroutine xc_vw

  
  
