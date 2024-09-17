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

subroutine hmin_hmax_gpu
  use gwm
  use simple_mpi, only : rank
  use device_mem_module

  implicit none
  integer :: st
  real*8  :: normhm(2),rayqt(2)

! Error handling
  if (.not.device_setup) stop ' hmin_hmax_gpu(): device not setup'
  if (.not.ffts_setup) stop ' hmin_hmax_gpu(): ffts not setup'
  if (.not.hmin_hmax_setup) call init_hmin_hmax_device

  !$acc data create(normhm) copyout(rayqt)

  call set_seed_general()
  call prep_hmin_hmax_gpu
  call run_hmin_hmax_gpu

  !$acc wait
  !$acc end data

  hmin = min(rayqt(1),rayqt(1)+rayqt(2))
  hmax = max(rayqt(1),rayqt(1)+rayqt(2))

  call flush_hmin_hmax_device

contains

  subroutine prep_hmin_hmax_gpu

    implicit none
    integer :: i,j

    call rand_c(po(:,:),2*n,dv)
    !$acc update device(po)

!   Set qcr,qci to zero before propagating to ensure that the correct 
!   subset of the grid is processed by chbflt_nl()
    !$acc data present(qcr,qci)
    !$acc parallel loop gang vector collapse(2) async
    do j=1,2
       do i=1,n
          qcr(i,j)=0.d0 ; qci(i,j)=0.d0
       enddo
    enddo
    !$acc end data

  end subroutine prep_hmin_hmax_gpu

  subroutine run_hmin_hmax_gpu

    implicit none
    integer :: i,j,k
    integer, parameter :: niter=200

    if (rank==0) then
       write(17,*) 'Hamiltonian spectral estimation (power iterations):'
    endif

!   Normalize po
    call renrmlz_hmin_hmax

!   Main Chebyshev loop
    !$acc data present(pe,po,rayqt)
    do k=0,niter

!      Both pe and po must contain Pn(H)*psi, so copy here
       !$acc parallel loop gang vector collapse(2) async
       do j=1,2
          do i=1,n
            pe(i,j)=po(i,j)
          enddo
       enddo

!      Iteration: calc 2*H*Tn(H) - Tn-1(H) using Chebyshev code,
!      but where Pn(H) is passed as both Tn ("pe") and Tn-1 ("po")
       call chebfilt_pteta(n,2,.TRUE.,.TRUE.)

!      Pn+1(H) =  0.5*(2*H*Pn(H) - Pn(H) + Pn(H))
       !$acc parallel loop gang vector collapse(2) async
       do j=1,2
          do i=1,n
             po(i,j) = 0.5d0*(po(i,j)+pe(i,j))
          enddo
       enddo

!      Second limit: apply shift after estimating first limit
       if (k>0) then
          !$acc parallel loop gang vector async
          do i=1,n
             po(i,2) = po(i,2)-rayqt(1)*pe(i,2)
          enddo
       endif

!      Compute Rayleigh coefficient to test convergence
       if (mod(k,10)==0) then
          call rayleigh_quotient
          !$acc wait
          !$acc update host (rayqt)
          if (rank==0) then
             write(17,*) k, rayqt(1),rayqt(2)
          endif
       endif

!      Normalize po
       call renrmlz_hmin_hmax

    enddo
    !$acc end data

  end subroutine run_hmin_hmax_gpu

  subroutine renrmlz_hmin_hmax
    implicit none
    integer :: i,j,iout
    real*8  :: nrm

    !$acc data present(normhm,pe,po)
    !$acc parallel loop gang async
    do j=1,2
       normhm(j) = 0.d0
    enddo

    !$acc parallel loop gang private(nrm) collapse(2) vector_length(256) async
    do j=1,2
       do iout = 1,n,256
          nrm=0.d0
          !$acc loop vector reduction(+:nrm) shortloop
          do i=iout,min(iout+255,n)
             nrm = nrm + abs(po(i,j))**2
          enddo
          !$acc atomic
          normhm(j) = normhm(j) + nrm
       enddo
    enddo

    !$acc parallel loop gang async
    do j=1,2
       normhm(j)=1.d0/sqrt(normhm(j)*dv)
    enddo

    !$acc parallel loop gang vector collapse(2) async
    do j=1,2
      do i=1,n
         po(i,j)=po(i,j)*normhm(j)
      enddo
    enddo

    !$acc end data

  end subroutine renrmlz_hmin_hmax

  subroutine rayleigh_quotient

    implicit none
    integer :: i,j,iout
    real*8  :: nrm

    !$acc data present(pe,po,rayqt)
    !$acc parallel loop gang async
    do j=1,2
       rayqt(j) = 0.d0
    enddo

    !$acc parallel loop gang private(nrm) collapse(2) vector_length(256) async
    do j=1,2
       do iout = 1,n,256
          nrm=0.d0
          !$acc loop vector reduction(+:nrm) shortloop
          do i=iout,min(iout+255,n)
             nrm = nrm + conjg(pe(i,j))*po(i,j)
          enddo
          !$acc atomic
          rayqt(j) = rayqt(j) + nrm
       enddo
    enddo

    !$acc parallel loop gang async
    do j=1,2
       rayqt(j) = rayqt(j)*dv
    enddo

    !$acc end data

  end subroutine rayleigh_quotient

end subroutine hmin_hmax_gpu
