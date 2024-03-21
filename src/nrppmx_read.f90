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
subroutine nrppmx_read
  use gwm
  use simple_mpi, only : rank, bcast_scalar_i
  use ppm,        only : nrppmx
  implicit none
  character(lngth_char) :: varn
  integer               :: inpstat
  logical               :: rdnrppmx       =.false.

  rnk0: if(rank==0) then
     nrppmx=2000
     
     open(101,file=trim(inputfname),iostat=inpstat)
     if(inpstat /= 0) stop "PROBLEM READING INPUT FILE"
     rewind(101)
     
     do 
        read(101,*,iostat=inpstat) varn
        if(inpstat /= 0) exit
        select case(varn)
        case('nrppmx')
           if(rdnrppmx)stop 'ERROR: ALREADY READ nrppmx '; rdnrppmx =.true.
           read(101,*,iostat=inpstat) nrppmx
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: nrppmx"
           if(nrppmx<1000.or.nrppmx>20000) stop "ERROR: nrppmx should be in 1,000-20,000 range" 
        end select
     end do

     close(101)
     end if rnk0
  call bcast_scalar_i(nrppmx)
end subroutine nrppmx_read
