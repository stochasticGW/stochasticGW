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
subroutine gam_seg1(shft, wdth)
  use gwm, only : seg, n
  implicit none
  integer shft, wdth, i, ib, ih
  call check_le(seg,n,' seg, n ')

  call rand_i(i, n+seg-1)
  shft= max(0,i-seg)
  ib = shft+1
  ih = min(i,      n)
  wdth = ih-ib+1
  call check_lelele(1,ib,ih,n,' 1-ib-ih-n ')

  ! Example: n=100, seg=10;  
  !  i = 1,...,10,11,12,... 99,100,101,......108,109
  !shft= 0,....0, 1,........89, 90, 91,.......98,99
  ! ib = 1,...,1, 2, 3,... ,90 ,91, 92,.......99,100
  ! ih = 1,...10, 11,12,....99,100,100,......100,100
  !wdth= 1,...10, 10,.......10, 10, 9,........2, 1

  ! note that an arbitrary site k is included (i.e.,  ib.le.k.le.ih ) seg times out of n+seg-1 
  ! (i.e., here 10 times out of 109).  For example, site 11 is included in i=11...20,
  ! site 1 included for i=1,..,10, site 99 for i=99...108, site 100 for i=100,..,109
  !
  ! so seg_fctr= 109/10 = (n+seg-1)/seg
end subroutine gam_seg1

subroutine gam_seg1_block
  use gwm, only : seg, n, ngam, seg_frctn, seg_sh, seg_w
  implicit none
  integer :: j,shft
  integer :: nseggrp,nsamples

  nseggrp=(n+seg-1)/seg
  nsamples=ngam/nseggrp

  do j=1,ngam
     shft=min((j-1)/nsamples*seg,n)
     seg_sh(j)=shft
     seg_w(j)=min(seg, max(n-shft,0))
!     write(*,*) j,seg_sh(j),seg_w(j)
  enddo

end subroutine gam_seg1_block
