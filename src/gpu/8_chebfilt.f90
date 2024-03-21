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

  subroutine chebfilt_pteta(n,ms,odd)

  use gwm,    only : na, ekn
  use kb_mod, only : ngs, mapkbg, vp_hamann, nsuper_ianl, indx_ianl, start_ianl
  use device_mem_module, only : map_sp => map_sp_pteta
  use device_mem_module, only : pe, po, ptmp, qcr, qci, shscvks, cufft_plan_pteta, two_ov_dh
  use device_mem_module, only : device_setup, ffts_setup, filter_setup
  use device_mem_module, only : scpsi
  implicit none
  integer, intent(in) :: n,ms
  logical, intent(in) :: odd

  ! The 'odd' flag replaces pe (.F.) or po (.T.); this is a workaround which is
  ! necessary since dynamic pointer mapping causes unwanted host->device copies

  if (.not.device_setup) stop ' chebfilt_pteta(): device not setup'
  if (.not.ffts_setup) stop ' chebfilt_pteta(): ffts not setup'
  if (.not.filter_setup) stop ' chebfilt_pteta(): filter module not setup'

  !$acc data present(pe,po,ptmp,ekn), &
  !$acc present(map_sp), &
  !$acc present(ngs,vp_hamann,mapkbg), &
  !$acc present(qcr,qci,shscvks,scpsi,indx_ianl,start_ianl)
 
  call chbflt_l
  call chbflt_nl
  call chbflt_k
  !$acc end data

contains

  subroutine chbflt_l
    implicit none
    integer :: i,is

    !$acc data present(pe,po,shscvks,map_sp)
    if (odd) then
       !$acc parallel loop gang vector async
       do i=1,n
          !$acc loop seq
          do is=1,ms
             po(i,is)=shscvks(i,map_sp(is))*pe(i,is)-po(i,is)
          enddo
       enddo
    else
       !$acc parallel loop gang vector async
       do i=1,n
          !$acc loop seq
          do is=1,ms
             pe(i,is)=shscvks(i,map_sp(is))*po(i,is)-pe(i,is)
          enddo
       enddo
    endif
    !$acc end data

  end subroutine chbflt_l

  subroutine chbflt_k
    use openacc
    use cudafor
    use cufft

    implicit none
    integer :: i,is
    integer(4) :: cerr
    real*8     :: ekntmp

    !$acc data present(pe,po,ptmp)
    if (odd) then
       !$acc host_data use_device(pe,ptmp)
       cerr = cufftExecZ2Z(cufft_plan_pteta,pe,ptmp,CUFFT_FORWARD)
       !$acc end host_data
       if (cerr/=0) stop ' chbflt_k(): problem with cutfft_forward (odd)'
    else
       !$acc host_data use_device(po,ptmp)
       cerr = cufftExecZ2Z(cufft_plan_pteta,po,ptmp,CUFFT_FORWARD)
       !$acc end host_data
       if (cerr/=0) stop ' chbflt_k(): problem with cutfft_forward (even)'
    endif

    !$acc data present(ekn)
    !$acc parallel loop gang vector private(ekntmp) async
    do i=1,n
       ekntmp=ekn(i)
       do is=1,ms
          ptmp(i,is)=ekntmp*ptmp(i,is)
       end do
    end do
    !$acc end data

    !$acc host_data use_device(ptmp)
    cerr = cufftExecZ2Z(cufft_plan_pteta,ptmp,ptmp,CUFFT_INVERSE)
    !$acc end host_data
    if (cerr/=0) stop ' chbflt_k(): problem with cutfft_inverse'

    if (odd) then
       !$acc parallel loop gang vector collapse(2) async
       do i=1,n
          do is=1,ms
             po(i,is)=po(i,is)+two_ov_dh*ptmp(i,is)
          end do
       end do
    else
       !$acc parallel loop gang vector collapse(2) async
       do i=1,n
          do is=1,ms
             pe(i,is)=pe(i,is)+two_ov_dh*ptmp(i,is)
          end do
       end do
    endif
    !$acc end data

  end subroutine chbflt_k

  subroutine chbflt_nl
    implicit none
    integer :: i, igg, l, is, it, ia, j
    real*8                  :: cer,cei
    complex*16, parameter   :: ci=(0d0,1d0)
    complex*16              :: ce

    !$acc data present(qcr,qci,pe,po,ngs,mapkbg)
    !$acc parallel loop gang collapse(2) async
    do is=1,ms
       do ia=1,na
          !$acc loop vector
          do igg=1,ngs(ia)
             it=mapkbg(igg,ia)
             qcr(it,is) = 0.0d0; qci(it,is) = 0.0d0
          end do
       end do
    end do

    !$acc data present(scpsi,vp_hamann,indx_ianl,start_ianl)
    if (odd) then
       !$acc parallel loop gang collapse(2) vector_length(160) async
       isoloop: do is=1,ms
          loloop: do l=1,nsuper_ianl
             ia = indx_ianl(l,1)

             ce=0d0
             !$acc loop vector reduction(+:ce)
             do igg=1,ngs(ia)
                ce=ce+pe(mapkbg(igg,ia),is)*vp_hamann(start_ianl(l)+igg)
             enddo
             ce=ce*scpsi(l)
             cer=dble(ce)
             cei=aimag(ce)

             !$acc loop vector
             do igg=1,ngs(ia)
                it=mapkbg(igg,ia)
                j=start_ianl(l)+igg
                !$acc atomic
                qcr(it,is) = qcr(it,is) + cer*vp_hamann(j)
                !$acc atomic
                qci(it,is) = qci(it,is) + cei*vp_hamann(j)
             end do
          end do loloop
       end do isoloop
    else
       !$acc parallel loop gang collapse(2) vector_length(160) async
       iseloop: do is=1,ms
          leloop: do l=1,nsuper_ianl
             ia = indx_ianl(l,1)

             ce=0d0
             !$acc loop vector reduction(+:ce)
             do igg=1,ngs(ia)
                ce=ce+po(mapkbg(igg,ia),is)*vp_hamann(start_ianl(l)+igg)
             enddo
             ce=ce*scpsi(l)
             cer=dble(ce)
             cei=aimag(ce)

             !$acc loop vector
             do igg=1,ngs(ia)
                it=mapkbg(igg,ia)
                j=start_ianl(l)+igg
                !$acc atomic
                qcr(it,is) = qcr(it,is) + cer*vp_hamann(j)
                !$acc atomic
                qci(it,is) = qci(it,is) + cei*vp_hamann(j)
             end do
          end do leloop
       end do iseloop
    endif
    !$acc end data

    if (odd) then
       !$acc parallel loop gang vector async
       do i=1,n
          do is=1,ms
             if ((qcr(i,is).ne.0.0d0).or.(qci(i,is).ne.0.0d0)) then
                po(i,is)=po(i,is)+qcr(i,is)+ci*qci(i,is)
             end if
          end do
       end do
    else
       !$acc parallel loop gang vector async
       do i=1,n
          do is=1,ms
             if ((qcr(i,is).ne.0.0d0).or.(qci(i,is).ne.0.0d0)) then
                pe(i,is)=pe(i,is)+qcr(i,is)+ci*qci(i,is)
             end if
          end do
       end do
    endif
    !$acc end data

  end subroutine chbflt_nl

end subroutine chebfilt_pteta
