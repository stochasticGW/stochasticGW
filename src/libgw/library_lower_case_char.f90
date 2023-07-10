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
subroutine lower_the_case(ch)
  implicit none
  character*(*) ch
  character*1 b
  integer i
  do i=1,len(ch)
     b=ch(i:i)
     if(iachar(b).ge.iachar('A').and.iachar(b).le.iachar('Z')) then
        b=achar(iachar(b)+iachar('a')-iachar('A'))
     endif
     ch(i:i)=b
  enddo
end subroutine lower_the_case

subroutine drv_lower_the_case
  implicit none
  character*30 a
  a=' aRpBBAAaa 32bB5D'
  write(6,*)' a pre lowering: ',a
  call lower_the_case(a)
  write(6,*)' a pst lowering: ',a
end subroutine drv_lower_the_case

!program lowrcase
!  call drv_lower_the_case
!end program lowrcase
