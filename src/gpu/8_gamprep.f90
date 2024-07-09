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
subroutine gam_prep_gpu
  use gwm
  use device_mem_module
  use curand

  implicit none
  integer(4) :: cerr
  real*8  :: fac
  integer :: i,j

! Error handling
  if (.not.rand_dev_setup) stop ' gam_prep_gpu(): curand not set up'

  fac=2.d0*dsqrt(3d0/dv)

  !$acc data copyout(gam)

  !$acc host_data use_device(gam)
  cerr = curandGenerate(curand_plan_gam, gam, seg*ngam)
  !$acc end host_data

  !$acc parallel loop gang vector collapse(2) async
  do j=1,ngam
     do i=1,seg
        gam(i,j)=fac*(dble(gam(i,j))-0.5d0)
     enddo
  enddo
  !$acc wait
  !$acc end data

!  write(*,*) 'Gamma test'
!  do i=1,50
!     write(*,*) i,gam(i,1),gam(i,2),gam(i,3),gam(i,4)
!  enddo

end subroutine gam_prep_gpu
