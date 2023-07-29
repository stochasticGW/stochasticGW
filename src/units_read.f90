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
subroutine read_units !(inputfname)
  use gwm
  use simple_mpi, only : rank
  
  implicit none
  character(10) :: varn
  integer       :: inpstat
  logical       :: rdunits_nuc    =.false.
  
  rnk0 : if(rank==0) then
     open(101,file=trim(inputfname),iostat=inpstat)
     rewind(101)
     if(inpstat /= 0) stop "PROBLEM READING INPUT FILE"

     units_nuc    = 'Bohr'
     do 
        read(101,*,iostat=inpstat) varn
        if(inpstat /= 0) exit
        call lower_the_case(varn)
        select case(varn)
        case('units_nuc')
           if(rdunits_nuc)stop ' ERROR: ALREADY READ units_nuc '; rdunits_nuc =.true.
           read(101,*,iostat=inpstat) units_nuc
           if(inpstat /= 0) stop "VALUE ABSENT FOR VAR: units_nuc"
        case default
        end select
     end do

     close(101)
     call assign_flg_units_nuc
  end if rnk0
  call bcast_units
end subroutine read_units

subroutine bcast_units
  use gwm
  use simple_mpi, only : bcast_scalar_i
  implicit none
  call bcast_scalar_i(flg_units_nuc)
end subroutine
