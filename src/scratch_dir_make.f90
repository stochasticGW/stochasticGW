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
subroutine scratch_dir_make
  use simple_mpi, only : nodes,rank, sync_mpi
  use gwm,        only : scratch_dir
  implicit none
  integer ir, stat
  logical fex
  character(len=20)  :: fr
  character(len=100) :: ch
  do ir=0,max(nodes-1,0)
     call sync_mpi
     call system(" sleep 0.01")
     if(rank==ir) then 
        inquire(file=trim(scratch_dir),exist=fex)
        if(.not.fex) then
           call system("mkdir "//trim(scratch_dir)//" > /dev/null 2>&1 ") 
        end if
     end if
     if(mod(ir,100)==0) call sleep(1)
  end do
  call sync_mpi

  call cnvrt_number_to_string_20(rank,fr)
  ch=trim(scratch_dir)//'/tmp.'//trim(fr)
  open(unit=99, iostat=stat, file=trim(ch), status='replace')
  if (stat == 0) then
     close(99, status='delete')
  else
     write(6,*) ' ERROR: problem opening '//trim(ch),' rank= ',rank
     stop
  endif
end subroutine scratch_dir_make
  
