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
subroutine open_ge_file_replace
  use simple_mpi
  use gwm, only : scratch_dir, icounter
  implicit none
  character(len=100) :: filename
  character(len=20)  :: fr,fi
  character(len=49)  :: fc

  call cnvrt_number_to_string_20(rank,fr)
  call cnvrt_number_to_string_20(icounter,fi)
  fc = '/ge_'//trim(adjustl(fr))//'.'//trim(adjustl(fi))//'.bin'
  filename = trim(adjustl(scratch_dir))//trim(adjustl(fc))
  open(unit=310000,file=trim(filename),status='replace',form='unformatted')
end subroutine open_ge_file_replace

subroutine open_ge_file_old
  use simple_mpi
  use gwm, only : scratch_dir, icounter
  implicit none
  character(len=100) :: filename
  character(len=20)  :: fr,fi
  character(len=49)  :: fc

  call cnvrt_number_to_string_20(rank,fr)
  call cnvrt_number_to_string_20(icounter,fi)
  fc = '/ge_'//trim(adjustl(fr))//'.'//trim(adjustl(fi))//'.bin'
  filename = trim(adjustl(scratch_dir))//trim(adjustl(fc))
  open(unit=310000,file=trim(filename),status='old',form='unformatted')
end subroutine open_ge_file_old



subroutine open_pt_file_replace
  use simple_mpi
  use gwm, only : scratch_dir, icounter
  implicit none
  character(len=100) :: filename
  character(len=20)  :: fr,fi
  character(len=49)  :: fc

  call cnvrt_number_to_string_20(rank,fr)
  call cnvrt_number_to_string_20(icounter,fi)
  fc = '/pt_'//trim(adjustl(fr))//'.'//trim(adjustl(fi))//'.bin'
  filename = trim(adjustl(scratch_dir))//trim(adjustl(fc))
  open(unit=320000,file=trim(filename),status='replace',form='unformatted')
end subroutine open_pt_file_replace

subroutine open_pt_file_old
  use simple_mpi
  use gwm, only : scratch_dir, icounter
  implicit none
  character(len=100) :: filename
  character(len=20)  :: fr,fi
  character(len=49)  :: fc

  call cnvrt_number_to_string_20(rank,fr)
  call cnvrt_number_to_string_20(icounter,fi)
  fc = '/pt_'//trim(adjustl(fr))//'.'//trim(adjustl(fi))//'.bin'
  filename = trim(adjustl(scratch_dir))//trim(adjustl(fc))
  open(unit=320000,file=trim(filename),status='old',form='unformatted')
end subroutine open_pt_file_old

