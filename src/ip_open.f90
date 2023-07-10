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
subroutine ip_open(i,ch,j,ip)
  use gwm
  use simple_mpi, only : rank
  implicit none
  integer ip, j, i
  integer ist
  character(len=20)  :: fr0,fj
  character(len=60)  :: fc
  character(len=100) :: fn
  character(len=3)   :: ch

  call cnvrt_number_to_string_20(rank,fr0)
  call cnvrt_number_to_string_20(j,fj)

  select case(i)
  case(0); fc='/v0_'//trim(adjustl(fj))//'_'//trim(adjustl(fr0))//'.bin'
  case(1); fc='/v1_'//trim(adjustl(fj))//'_'//trim(adjustl(fr0))//'.bin'
  case(2); fc='/v2_'//trim(adjustl(fj))//'_'//trim(adjustl(fr0))//'.bin'
  case default; stop ' ip_open '
  end select
  fn = trim(adjustl(scratch_dir))//trim(adjustl(fc))
  select case(ch)
  case('rpl')
     open(ip,file=trim(adjustl(fn)),status='replace',form='unformatted')
  case('old')
     open(ip,file=trim(adjustl(fn)),status='old',form='unformatted')
  case default
     stop ' ip_open '
  end select

  rewind(ip)

end subroutine ip_open
