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
subroutine cnvrt_number_1K_string(i,fr)
  implicit none
  character(len=1024) :: fr
  integer i

  fr= ""

  select case(i)
     case(0  :9 );       write(fr,"(I1)")i
     case(10 :99);       write(fr,"(I2)")i
     case(100:999);      write(fr,"(I3)")i
     case(1000:9999);    write(fr,"(I4)")i
     case(10000:99999);  write(fr,"(I5)")i
     case(100000:999999);write(fr,"(I6)")i
     case default; write(6,*)' problem with i in cnvrt_number_1K_string, =',i; stop
  end select

end subroutine cnvrt_number_1K_string
