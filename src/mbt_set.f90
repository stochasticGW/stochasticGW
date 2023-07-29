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
subroutine set_mf(mf)
  use gwm, only : n, nvr
  implicit none
  integer mf
  mf  = n/nvr
  if(mf*nvr<n)mf=mf+1
end subroutine set_mf

subroutine mbt_set(j, mb, mt)
  use simple_mpi, only : cr=> color_rank
  use simple_mpi, only : cs=> color_size
  use gwm
  implicit none
  integer :: j, mb, mt
  integer :: mf, ivr

  ivr = j + nvr/cs* cr   

  call checks
  call set_mf(mf)
  mb  = (ivr-1)*mf+1
  mt  = min(ivr*mf,n)
!  call check(mt-mb+1,mf,' mtbf ')
contains
  subroutine checks
    implicit none
    call check0(mod(nvr,cs),' mod_nvr_cs ')
    call check_le(ivr,nvr,  ' ivr, nvr   ')
  end subroutine checks
end subroutine mbt_set
  
