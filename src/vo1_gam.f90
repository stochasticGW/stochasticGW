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
subroutine vo1_gam(it)
  use gwm
  implicit none
  complex*16, allocatable:: u(:)
  complex*16 :: c1
  integer i, j, it, st, sh
  if(ngam<1) stop ' ngam<1 in prepvo '

  vo = 0d0

  select case(gamflg)
  case(1)
     do i=1,n
        vo(i,sp0) = sum(gam(i,:)*ft(it,:))/ dble(ngam_nzero_blk)
     enddo
  case(2)
     allocate( u(n), stat=st)
     u = 0d0
     do j=1,ngam
        sh = seg_sh(j)
        c1 = ft(it,j)
        do i=1,seg_w(j)
           u(i+sh) = u(i+sh) + c1*gam(i,j)
        enddo
     enddo
     vo(:,sp0) = u(:) *  (seg_fctr/dble(ngam_nzero_blk))
     deallocate(u)
  case default
     stop ' ERROR: gamflg not 1 nor 2 '
  end select
end subroutine vo1_gam
