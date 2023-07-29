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
subroutine set_ngam_nzero_blk
  use simple_mpi, only : cs=>color_size
  use simple_mpi, only : rank, bcast_scalar_i
  use gwm, only : n,seg,ngam_blk,ngam_nzero_blk,block_gam_alg
  implicit none
  integer :: i,nspanning_segments,ng,ng_blk,nsamples
  integer :: ng_blk_next,ngam_nzero_blk_next

! Number of segments to fully span the spatial grid
  nspanning_segments=(n+seg-1)/seg

  if(rank==0) then

    if (block_gam_alg) then

!      Compute the number of gammas for each color
       ngam_nzero_blk=0
       ngam_nzero_blk_next=0
       do i=1,cs
          ng=ngam_blk/cs
          if(i<mod(ngam_blk, cs)) ng=ng+1
!         Number of sampling gammas for each segment on this processor
          nsamples=ng/nspanning_segments
          ng_blk=nspanning_segments*nsamples
          ng_blk_next=nspanning_segments*(nsamples+1)
          ngam_nzero_blk = ngam_nzero_blk + ng_blk
          ngam_nzero_blk_next = ngam_nzero_blk_next + ng_blk_next
       enddo

       if (ngam_nzero_blk.lt.ngam_blk) then
          write(6,*) 'NOTE: block gamma algorithm quantizes nxi: '
          write(6,'(X,3(A,X,I0,X),A/)') '[',ngam_nzero_blk,&
          '(SELECTED) --',ngam_blk,'(input) --',ngam_nzero_blk_next,'(next)]'
          write(17,*) 'NOTE: block gamma algorithm quantizes nxi: '
          write(17,'(X,3(A,X,I0,X),A/)') '[',ngam_nzero_blk,&
          '(SELECTED) --',ngam_blk,'(input) --',ngam_nzero_blk_next,'(next)]'
       endif
    else
       ngam_nzero_blk = ngam_blk
    endif

  endif

  call bcast_scalar_i(ngam_nzero_blk)

end subroutine set_ngam_nzero_blk

