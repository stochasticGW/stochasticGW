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
subroutine scratch_dir_name_read
  use gwm
  use simple_mpi, only : rank, bcast_scalar_char, sync_mpi
  implicit none
  character*1          :: a
  character(40)        :: varn
  integer              :: inpstat  
  character*20 ch
  call cnvrt_number_to_string_20(icounter,ch)

  rnk0 : if(rank==0) then

     scratch_dir = './GW_SCRATCH.'//trim(adjustl(ch))
     open(101,file=trim(inputfname),iostat=inpstat)
     if(inpstat /= 0) stop "PROBLEM READING INPUT FILE"
     rewind(101)
     do 
        read(101,*,iostat=inpstat) varn
        if(inpstat /= 0) exit
        call lower_the_case(varn)
        select case(trim(adjustl(varn)))
        case('scratch_path')
           if(rdscratch_path)stop ' ERROR: ALREADY READ scratch_path '; rdscratch_path =.true.
           read(101,"(A)",iostat=inpstat) scratch_path
           scratch_path=trim(adjustl(scratch_path))
           write(6,*)' scratch_path=',scratch_path
           if(inpstat /= 0) stop "ERROR: VALUE ABSENT FOR VAR: scratch_path"
           scratch_dir = trim(adjustl(scratch_path))//'/GW_SCRATCH.'//trim(adjustl(ch))
           write(6,*)' scratch_dir=',scratch_dir
        case default
        end select
     end do
     close(101)
     
  end if rnk0
  call sync_mpi
  call bcast_scratchdir
contains 
  subroutine bcast_scratchdir
    use simple_mpi, only : rank
    implicit none
    integer i
    do i=1,lngth_char
       if(rank==0) then
          a=scratch_dir(i:i)
       else
          a=' '
       endif
       call bcast_scalar_char(a)
       scratch_dir(i:i)=a
    enddo
  end subroutine bcast_scratchdir
end subroutine scratch_dir_name_read
