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
subroutine vk_prep
  use gwm
  implicit none
  integer nxb, nyb, nzb

  nxb = scale_vh * nx
  nyb = scale_vh * ny
  nzb = scale_vh * nz

  call vk_prep_b(nxb,nyb,nzb)

end subroutine vk_prep

subroutine vk_prep_b(nxb,nyb,nzb)
  use gwm
  implicit none
  integer nxb, nyb, nzb, ngb,nd,ix,iy,iz,i
  real*8, allocatable:: grdb(:,:,:,:)
  nd = 3
  ngb = nxb*nyb*nzb; if(abs(ngb-dble(nxb)*dble(nyb)*dble(nzb))>1d-4) stop ' ngb not proper ; maybe need integer*8 '

  if(allocated(vk)) deallocate(vk);  allocate(vk(ngb), stat=i); if(i/=0) stop ' vkb '

  allocate(grdb(nxb,nyb,nzb,3), stat=i); if(i/=0) stop ' grdb '

  do iz=1,nzb
     do iy=1,nyb
        do ix=1,nxb
           grdb(ix,iy,iz,:) = (/ (ix-1-nxb/2)*dx, (iy-1-nyb/2)*dy, (iz-1-nzb/2)*dz /); 
        enddo
     enddo
  enddo
  
  call prep_vk( nxb,nyb,nzb,dx,dy,dz,vk)
  vk_scale_cnst = 1d0/dble(nxb)/dble(nyb)/dble(nzb)
  vk= vk * vk_scale_cnst
  deallocate(grdb)
end subroutine vk_prep_b

