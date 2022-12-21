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
subroutine vs_read_rmv(  j,mb,mt,i,vs)
  use gwm
  use simple_mpi, only : rank
  implicit none
  integer j, mb, mt, i, it
  integer mmb, mmt, nnt, itt, mspv, ip
  real*4  vs(mb:mt,1:nspv,0:nt)

  if(mb>n)return
  call check_le(mb,mt,' mb,mt ')
  call ip_make(i,j,ip)

  read(ip)nnt,mmb,mmt,mspv
  call check(nnt, nt,  ' nt_beg ')
  call check(mmb,mb,   ' mmb    ')
  call check(mmt,mt,   ' mmt    ')
  call check(nspv,mspv,' nmspv  ')

  do it=0,nt
     read(ip)itt
     call check(it,itt,' it-itt-vs ')
     read(ip)vs(mb:mt,:,it)
  enddo
  read(ip)nnt
  call check(nnt, nt,' nt_end ')

  close(ip,status='delete')
end subroutine vs_read_rmv
  
  
subroutine vs_rmv(  j,i)
  implicit none
  integer j, i, ip
  call ip_make(i,j,ip)
  close(ip,status='delete')
end subroutine vs_rmv
  
  

