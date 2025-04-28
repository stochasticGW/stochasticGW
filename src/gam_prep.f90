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
subroutine gam_prep
  use gwm
#if GPU_ENABLED
  use device_mem_module
  use seed_gw_modu, only : seed_array, line_seed_gam
  use curand
#endif
  implicit none
  integer :: st,i
  integer(4) :: cerr

  if(ngam<1) return
  if(gamflg==2) call gam_seg
  if(.not.allocated(gam)) then
     allocate(gam(seg,ngam),stat=st); call check0(st,' gam ')
  endif

#if GPU_ENABLED
  if (disable_gpu_gam) then
     call gam_prep_cpu() ! CPU version
  else
!    Error handling
     if (.not.rand_dev_setup) stop ' gam_prep(): curand not set up'
!    Get the seed from the kiss generator and pass to cuRAND
!    so that each sample has a different seed
     call set_seed_gam(1)
     cerr = curandSetPseudoRandomGeneratorSeed(curand_plan_gam,seed_array(1,line_seed_gam))

     if (use_host_curand) then
        call gam_prep_curand_cpu()
     else
        call gam_prep_curand_gpu()
     endif
  endif
#else
  call gam_prep_cpu()
#endif

  if (rank==0) then
#if GPU_ENABLED
     if (disable_gpu_gam) then
        write(17,*) 'Gamma test (intrinsic RNG on host)'
     else
        if (use_host_curand) then
           write(17,*) 'Gamma test (curand on host)'
        else
           write(17,*) 'Gamma test (curand on device)'
        endif
     endif
#else
     write(17,*) 'Gamma test (intrinsic RNG on host)'
#endif
!     do i=1,50
!        write(*,*) i,gam(i,1),gam(i,2),gam(i,3),gam(i,4)
!     enddo
  endif

end subroutine gam_prep

subroutine gam_prep_cpu
  use gwm
  implicit none
  integer igam

  do igam=1,ngam
     call set_seed_gam( igam)
     call rand_r48(gam(:,igam),size(gam,1), dv)
  end do

end subroutine gam_prep_cpu
