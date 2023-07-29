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
subroutine dealloc_gw
  use gwm

  if(allocated(vks))     deallocate(vks)
  if(allocated(vxc))     deallocate(vxc)
  if(allocated(vxc0))    deallocate(vxc0)
  if(allocated(ear))     deallocate(ear)
  if(allocated(rho))     deallocate(rho)
  if(allocated(vt))      deallocate(vt)
  if(allocated(ct))      deallocate(ct)
  if(allocated(cveta))   deallocate(cveta)
  if(allocated(cvxi))    deallocate(cvxi)
  if(allocated(normpt))  deallocate(normpt)
  if(allocated(ekn))     deallocate(ekn)
  if(allocated(expnk))   deallocate(expnk)
  if(allocated(ft))      deallocate(ft)
  if(allocated(gamvr))   deallocate(gamvr)
end subroutine dealloc_gw
