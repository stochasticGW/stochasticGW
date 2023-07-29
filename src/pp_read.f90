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
subroutine pp_read
  use gwm
  use simple_mpi, only : rank, bcast_scalar_i
  implicit none
  character(lngth_char) :: varn
  integer               :: inpstat
  logical               :: rdpps =.false.

  rnk0: if(rank==0) then
     open(101,file=trim(inputfname),iostat=inpstat)
     if(inpstat /= 0) stop "PROBLEM READING INPUT FILE"
     rewind(101)

     do 
        read(101,*,iostat=inpstat) varn
        if(inpstat/=0) exit
        select case(varn)
        case('pps')
           if(rdpps)        stop 'ERROR: ALREADY READ pps '; rdpps =.true.
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: pps"
           call pp_read_in
        end select
     end do
     close(101)
     if(.not.rdpps) stop "ERROR: NEED TO READ PSEUDOPOTENTIAL!"
  end if rnk0

end subroutine pp_read


