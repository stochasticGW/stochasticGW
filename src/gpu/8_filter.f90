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

subroutine pt_eta_fltr_gpu
  use gwm
  use simple_mpi, only : rank, color_rank
  use device_mem_module

  implicit none
  integer :: st

! Error handling
  if (.not.device_setup) stop ' pt_eta_gpu(): device not setup'
  if (.not.ffts_setup) stop ' pt_eta_gpu(): ffts not setup'
  if (.not.filter_setup) call init_filter_device

  if (.not.allocated(pt)) stop ' pt_eta_gpu(): pt not allocated'
  if (.not.allocated(zeta)) stop ' pt_eta_gpu(): zeta not allocated'
  if (.not.allocated(eta)) stop ' pt_eta_gpu(): eta not allocated'
  if (.not.allocated(th_co)) stop ' pt_eta_gpu(): th_co not allocated'

  if (SIZE(pt,1).ne.n .or. SIZE(pt,2).ne.ns) then
     write(*,'(X,4(A,I0),A)') 'pt must be (',n,' x ',ns,&
                              ') but is (',SIZE(pt,1),' x ',SIZE(pt,2),')'
     stop ' pt_eta_gpu(): wrong size: pt'
  endif
  if (SIZE(zeta).ne.n) then
     write(*,'(X,2(A,I0),A)') 'zeta must be (',n,') but is (',SIZE(zeta),')'
     stop ' pt_eta_gpu(): wrong size: zeta'
  endif
  if (SIZE(eta).ne.n) then
     write(*,'(X,2(A,I0),A)') 'eta must be (',n,') but is (',SIZE(eta),')'
     stop ' pt_eta_gpu(): wrong size: eta'
  endif
  if (SIZE(th_co).ne.nchbc) then
     write(*,'(X,2(A,I0),A)') 'th_co must be (',nchbc,') but is (',SIZE(th_co),')'
     stop ' pt_eta_gpu(): wrong size: th_co'
  endif
  if (SIZE(normpt).ne.ns) then
     write(*,'(X,2(A,I0),A)') 'normpt must be (',ns,') but is (',SIZE(normpt),')'
     stop ' pt_eta_gpu(): wrong size: normpt'
  endif

  allocate(pa(n,ns+1), stat=st); &
           if(st/=0) stop ' init_filter_device(): trouble allocating pa '

! Pack 'pt' and 'zeta' arrays into module-array 'pa' since we do the same
! filtering operation to both in parallel
  pe(:,1:ns)=pt(:,1:ns)
  pe(:,ns+1)=zeta(:)
  !$acc update device(pe)
  !$acc data copyout(pa,normpt)

  call prep_chebyshev_gpu
  call filter_chebyshev_gpu
  call renrmlz_cheb
  !$acc wait
  !$acc end data

! Unpack accumulator on host
  pt(:,1:ns)=pa(:,1:ns)
  if(color_rank==eta_rank) eta(:)=pa(:,ns+1)

  deallocate(pa, stat=st); &
             if(st/=0) stop ' pt_eta_gpu(): trouble deallocating pa '

  call flush_filter_device

contains

  subroutine prep_chebyshev_gpu

    implicit none
    integer :: i,is

!   Initialize arrays for Chebyshev recursion
    !$acc data present(pe,po,pa)
    !$acc parallel loop gang vector collapse(2) async
    do is=1,ns+1
       do i=1,n
          pa(i,is)=0.d0; po(i,is)=pe(i,is)
       enddo
    enddo
    !$acc end data

!   Set qcr,qci to zero before propagating to ensure that the correct 
!   subset of the grid is processed by chbflt_nl()
    !$acc data present(qcr,qci)
    !$acc parallel loop gang vector collapse(2) async
    do is=1,ns+1
       do i=1,n
          qcr(i,is)=0.d0 ; qci(i,is)=0.d0
       enddo
    enddo
    !$acc end data

  end subroutine prep_chebyshev_gpu

  subroutine filter_chebyshev_gpu

    implicit none
    integer :: i,is,k
    real*8  :: coef

!   Main Chebyshev loop
    !$acc data present(pe,po,pa)
    do k=0,nchbc-1
       coef=th_co(k)

       if (mod(k,2)==0) then

          !$acc parallel loop gang vector collapse(2) async
          do is=1,ns+1
             do i=1,n
                pa(i,is) = pa(i,is)+coef*pe(i,is)
             enddo
          enddo
          call chebfilt_pteta(n,ns+1,.TRUE.,.FALSE.)

          if (k.eq.0) then
             !$acc parallel loop gang vector collapse(2) async
             do is=1,ns+1
                do i=1,n
                   po(i,is) = 0.5d0*(po(i,is)+pe(i,is))
                enddo
             enddo
          endif

       else

          !$acc parallel loop gang vector collapse(2) async
          do is=1,ns+1
             do i=1,n
                pa(i,is) = pa(i,is)+coef*po(i,is)
             enddo
          enddo
          call chebfilt_pteta(n,ns+1,.FALSE.,.FALSE.)

       endif

    enddo
    !$acc end data

  end subroutine filter_chebyshev_gpu

  subroutine renrmlz_cheb
    implicit none
    integer :: i,is,iout
    real*8  :: nrm

!   Note that 'pt' in cols (1:ns) of 'pa' should be renormalized but
!   'eta' in column ns+1 should not be, so 'is' loop runs only (1:ns)

    !$acc data present(pa,normpt)
    !$acc parallel loop gang async
    do is=1,ns
       normpt(is) = 0.d0
    enddo

    !$acc parallel loop gang private(nrm) collapse(2) vector_length(256) async
    do is=1,ns
       do iout = 1,n,256
          nrm=0.d0
          !$acc loop vector reduction(+:nrm) shortloop
          do i=iout,min(iout+255,n)
             nrm = nrm + abs(pa(i,is))**2
          enddo
          !$acc atomic
          normpt(is) = normpt(is) + nrm
       enddo
    enddo

    !$acc parallel loop gang async
    do is=1,ns
       normpt(is)=sqrt(normpt(is)*dv)
    enddo

    !$acc parallel loop gang vector collapse(2) async
    do is=1,ns
      do i=1,n
         pa(i,is)=pa(i,is)/normpt(is)
      enddo
    enddo
    !$acc end data

  end subroutine renrmlz_cheb

end subroutine pt_eta_fltr_gpu
