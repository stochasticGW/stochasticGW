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
subroutine cnt_check_shift
  use atoms,     only : cnt
  use gwm,       only : na
  use gwm,       only : nx, ny, nz
  use gwm,       only : dx, dy, dz
  use gwm,       only : boxpos
  use simple_mpi,only : rank,  bcast_r8

  implicit none
  integer     :: ia
  real*8      :: x,y,z

  rank_if : if(rank==0) then
     atom_do : do ia=1,na
        x=cnt(1,ia)
        y=cnt(2,ia)
        z=cnt(3,ia)
        box: if(boxpos) then
           if(x<0d0.or.x>nx*dx.or.y<0d0.or.y>ny*dy.or.z<0d0.or.z>nz*dz)then
              write(6,*)' ERROR: Problem in atom on line ',ia,' of cnt.ini '
              write(6,*)'  Specifically, (possibly scaled) coordinates ',real(x),real(y),real(z) 
              write(6,*)'  are outside box of zero to nx*dx,ny*dy,nz*dz ',nx*dx,ny*dy,nz*dz
              call flush(6)
              stop
           endif
           x=x-dx*(nx/2)
           y=y-dy*(ny/2)
           z=z-dz*(nz/2)
           cnt(:,ia)=(/x,y,z/)
           write(17,*)' atom ',ia,' ; modified cnt(:,ia) ',real(cnt(:,ia))
        else
           if(abs(x)>nx*dx/2d0.or.abs(y)>ny*dy/2d0.or.abs(z)>nz*dz/2d0)then
              write(6,*)' ERROR: Problem in atom on line ',ia,' of cnt.ini '
              write(6,*)'  Specifically,(possibly scaled) coordinates ',real(x),real(y),real(z) 
              write(6,*)'  are outside 0-centered box of half-lengths ',&
                   real(nx*dx/2d0),real(ny*dy/2d0),real(nz*dz/2d0)
              call flush(6)
              stop
           endif
        end if box
     end do atom_do
  end if rank_if

  call bcast_r8(cnt,size(cnt),0)
end subroutine cnt_check_shift
