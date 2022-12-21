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
subroutine makerho
  use gwm
  use simple_mpi, only : color_allsum_r8
  implicit none
  integer i, is, sp
  real*8  a
  dyn: if(flgdyn) then ! need full rho
     call makerho_full
     call rho_sum
  else
     call makerho_p  ! note that then rho is not defined
  end if dyn

contains

  subroutine makerho_full
    implicit none

    rho = 0d0
    
    do is=1,ns
       if(det_tddft) then
          a = occ_det(is)            ! no 1/ns, no 2/nsp, since absorbed in occ_det
       else
          a = 2d0/nsp*normpt(is)**2d0/ns_blk 
       end if
       sp = map_sp_pt(is)
       
       do i=1,n
          rho(i,sp) = rho(i,sp)+a*( (dble(pt(i,is)))**2+(aimag(pt(i,is)))**2 )
       end do
    enddo
    
    call color_allsum_r8(rho, size(rho))
    rho = rho /sum(rho*dv) * rnel ! normalizing rho!!!!
  end subroutine makerho_full

  subroutine rho_sum
    implicit none
    integer st
    if(.not.allocated(rho_p)) then
       allocate(rho_p(n), stat=st)
       call check0(st,' rho_p ')
    end if

    rho_p = 0d0
    do sp=1,nsp
       rho_p = rho_p+ rho(:,sp)
    enddo
  end subroutine rho_sum

  subroutine makerho_p
    implicit none
    integer st
    if(.not.allocated(rho_p)) then
       allocate(rho_p(n), stat=st)
       call check0(st,' rho_p ')
    end if
    rho_p = 0d0
    
    do is=1,ns
       if(det_tddft) then
          a = occ_det(is)            ! no 1/ns, no 2.
       else
          a = 2d0/nsp*normpt(is)**2d0/ns_blk  ! spin: 2.  
       end if
       
       do i=1,n
          rho_p(i) = rho_p(i)+a*( (dble(pt(i,is)))**2+(aimag(pt(i,is)))**2 )
       end do
    enddo
    
    call color_allsum_r8(rho_p, size(rho_p))
    rho_p = rho_p /sum(rho_p*dv) * rnel ! normalizing rho!!!!
  end subroutine makerho_p
    
end subroutine makerho

