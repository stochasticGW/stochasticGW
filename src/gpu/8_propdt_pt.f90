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

  subroutine propdt_pt(n,ms,short)

  use gwm,    only : dv, dt, na, expnk
  use gwm,    only : map_sp => map_sp_pt
  use kb_mod, only : ngs, mapkbg, vp_hamann, nsuper_ianl, indx_ianl, start_ianl
  use device_mem_module, only : pt0, expvks, qr, qi, nrmarray, cufft_plan
  use device_mem_module, only : device_setup, ffts_setup, propdtpt_setup
  use device_mem_module, only : ovdev
  implicit none
  integer, intent(in) :: n,ms
  logical, intent(in) :: short

  if (.not.device_setup) stop ' propdt_pt(): device not setup'
  if (.not.ffts_setup) stop ' propdt_pt(): ffts not setup'
  if (.not.propdtpt_setup) stop ' propdt_pt(): propdtpt not setup'

  !$acc data present(pt0,expnk), &
  !$acc present(map_sp), &
  !$acc present(ngs,vp_hamann,mapkbg), &
  !$acc present(qr,qi,nrmarray,expvks,ovdev,indx_ianl,start_ianl)
  
  call propl
  call propnl
  call propk
  call propnl
  call propl
  call renrmlz
  !$acc end data

contains

  subroutine propl
    implicit none
    integer :: j,is,i
    complex*16 :: exptmp

!   Shortcut potential step for etaxi single-spin (sp0) propagation
!   Here, ms=1 and the same expvks is used for eta and xi, which are in 
!   pt(:,1,1) and pt(:,1,2), respectively
    if (short) then
       !$acc data present(pt0,expvks)
       !$acc parallel loop gang vector collapse(2) async
       do i=1,n
          do j=1,2
             pt0(i,1,j)=pt0(i,1,j)*expvks(i,1,1)
          enddo
       enddo
       !$acc end data

!   Full potential step for pt propagation
    else
       !$acc data present(pt0,expvks,map_sp)
       !$acc parallel loop gang vector async
       do i=1,n
          !$acc loop seq
          do j=1,2
             !$acc loop seq
             do is=1,ms
                pt0(i,is,j)= pt0(i,is,j)*expvks(i,map_sp(is),j)
             enddo
          enddo
       enddo
       !$acc end data
    endif

  end subroutine propl

  subroutine propk
    use openacc
    use cudafor
    use cufft

    implicit none
    integer :: i,j,is,plan,iout,ii
    integer(4) :: cerr
    complex*16 :: exptmp

!   Select the cufft plan
    if (short) then
       plan=2
    else
       plan=1
    endif

    !$acc data present(pt0)
    !$acc host_data use_device(pt0)
    cerr = cufftExecZ2Z(cufft_plan(plan),pt0,pt0,CUFFT_FORWARD)
    !$acc end host_data
    if (cerr/=0) stop ' propk(): problem with cutfft_forward'
    !$acc data present(expnk)
    !$acc parallel loop gang vector private(exptmp) async
    do i=1,n
       exptmp=expnk(i)
       do is=1,ms
          do j=1,2
             pt0(i,is,j)=exptmp*pt0(i,is,j)
          end do
       end do
    end do
    !$acc end data
    !$acc host_data use_device(pt0)
    cerr = cufftExecZ2Z(cufft_plan(plan),pt0,pt0,CUFFT_INVERSE)
    !$acc end host_data
    if (cerr/=0) stop ' propk(): problem with cutfft_inverse'

    !$acc end data
  end subroutine propk

  subroutine propnl
    implicit none
    integer :: i, igg, l, j, is, it, ia
    real*8                  :: cer,cei
    complex*16, parameter   :: ci=(0d0,1d0)
    complex*16              :: ce

    !$acc data present(qr,qi,pt0,ngs,mapkbg)
    !$acc parallel loop gang collapse(3) async
    do i=1,2
       do is=1,ms
          do ia=1,na
             !$acc loop vector
             do igg=1,ngs(ia)
                it=mapkbg(igg,ia)
                qr(it,is,i) = 0.0d0; qi(it,is,i) = 0.0d0
             end do
          end do
       end do
    end do

    !$acc data present(ovdev,vp_hamann,indx_ianl,start_ianl)
    !$acc parallel loop gang collapse(3) vector_length(160) async
    do i=1,2
       isloop : do is=1,ms
          lloop: do l=1,nsuper_ianl
             ia = indx_ianl(l,1)

             ce=0d0
             !$acc loop vector reduction(+:ce)
             do igg=1,ngs(ia)
                 ce=ce+pt0(mapkbg(igg,ia),is,i)*vp_hamann(start_ianl(l)+igg)
             enddo
             ce=ce*ovdev(l)
             cer=dble(ce)
             cei=aimag(ce)

             !$acc loop vector
             do igg=1,ngs(ia)
                it=mapkbg(igg,ia)
                j=start_ianl(l)+igg
                !$acc atomic
                qr(it,is,i) = qr(it,is,i) + cer*vp_hamann(j)
                !$acc atomic
                qi(it,is,i) = qi(it,is,i) + cei*vp_hamann(j)
             end do
          end do lloop
       end do isloop
    end do ! i loop

    !$acc end data

    !$acc parallel loop gang vector async
    do i=1,n
       do is=1,ms
          do j=1,2
             if ((qr(i,is,j).ne.0.0d0).or.(qi(i,is,j).ne.0.0d0)) then
                pt0(i,is,j)=pt0(i,is,j)+qr(i,is,j)+ci*qi(i,is,j)
             end if
          end do
       end do
    end do
    !$acc end data

  end subroutine propnl

  subroutine renrmlz
    implicit none
    integer :: i,j,is,iout
    real*8  :: nrm

    !$acc data present(pt0,nrmarray)
    !$acc parallel loop gang collapse(2) async
    do j=1,2
      do is=1,ms
         nrmarray(is,j) = 0.d0
      enddo
    enddo

    !$acc parallel loop gang private(nrm) collapse(3) vector_length(256) async
    do j=1,2
       do is=1,ms
          do iout = 1,n,256
             nrm=0.d0
             !$acc loop vector reduction(+:nrm) shortloop
             do i=iout,min(iout+255,n)
                nrm = nrm + abs(pt0(i,is,j))**2
             enddo
             !$acc atomic
             nrmarray(is,j) = nrmarray(is,j) + nrm
          enddo
       enddo
    enddo

    !$acc parallel loop gang collapse(2) async
    do j=1,2
      do is=1,ms
         nrmarray(is,j)=1.d0/sqrt(nrmarray(is,j)*dv)
      enddo
    enddo

    !$acc parallel loop gang vector collapse(3) async
    do j=1,2
      do is=1,ms
        do i=1,n
           pt0(i,is,j)=pt0(i,is,j)*nrmarray(is,j)
        enddo
      enddo
    enddo
    !$acc end data

  end subroutine renrmlz

end subroutine propdt_pt
