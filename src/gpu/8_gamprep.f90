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
subroutine gam_prep_curand_gpu
  use gwm
  use device_mem_module
  use curand

  implicit none
  integer(4) :: cerr
  real*8  :: fac
  integer :: i,j

  if (.not.rand_dev_setup) stop ' gam_prep_curand_gpu(): curand not set up'
  if (.not.gam_dev_setup) call init_gam_device

  fac=2.d0*dsqrt(3d0/dv)

  ! acc data copyout(gam)
  !$acc data present(gam)

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

end subroutine gam_prep_curand_gpu
