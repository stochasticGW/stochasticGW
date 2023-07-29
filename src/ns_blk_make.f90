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
subroutine ns_blk_make
  use simple_mpi, only : rank,bcast_scalar_i, color_size
  use gwm
  implicit none
  integer nold
  rnk0 : if(rank==0) then
     write(6,*)
     dettddft : if(det_tddft) then

        if(rnel<0.5d0) stop ' problem, rnel too small '
        nold = ns_blk
        ns_blk = ceiling(rnel*dble(nsp)/2d0-1d-8)
        if(mod(ns_blk,nsp)/=0)ns_blk=ns_blk+1

        write(6,*)" NOTE: deterministic tddft, # states modifed from ",nold," to: ",ns_blk
        if(color_size>ns_blk+1) then
           write(6,*)' buffer_size in input (=',color_size,') is too large, '
           write(6,*)' since only ',ns_blk,' deterministic states are used '
           write(6,*)' Therefore, stopping calculation. Repeat with smaller buffer_size (e.g., 1)'
           call flush(6)
           stop
        end if

     else

        if(ns_blk<color_size-1) then
           nold = ns_blk
           ns_blk=color_size-1
           if(mod(ns_blk,nsp)/=0)ns_blk =ns_blk+1
           write(6,*)' NOTE: number_tddft_states was ',nold,&
                ' < buffer_size-1 =',color_size-1,' so increased to  ',ns_blk
        endif
     end if dettddft
     call flush(6)
  end if rnk0
  call bcast_scalar_i(color_size)
  call bcast_scalar_i(ns_blk)
end subroutine ns_blk_make
