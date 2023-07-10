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
subroutine open_wf_file(ip)
  use simple_mpi, only : rank
  use gwm,        only : binwf
  implicit none
  integer ip, ios
  rnk0 : if(rank==0) then
     if(binwf) then
        open ( ip,file='wf.bin',status='old',form='unformatted',iostat=ios)
        call check0(ios,' ios in opening wf.bin ')
     else
        open ( ip,file='wf.txt',status='old',iostat=ios)
        call check0(ios,' ios in opening wf.txt ')
     end if
     rewind(ip)
  end if rnk0
end subroutine open_wf_file

subroutine close_wf_file(ip)
  use simple_mpi, only : rank
  implicit none
  integer ip, ios
  rnk0 : if(rank==0) then
     close ( ip,iostat=ios)
     call check0(ios,' ios in closing wf.txt ')
  end if rnk0
end subroutine close_wf_file
