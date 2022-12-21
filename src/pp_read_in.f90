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
subroutine pp_read_in
  use simple_mpi, only : rank
  use atoms, only : matop, mapat
  implicit none
  integer       :: ich, ival
  integer       :: inpstat
  logical       :: fex 
  character(50) :: ppfname

  call check0(rank,' rank-0 ')

  write(6,*) " ############ READING PSEUDOPOTENTIAL FILES ############"
  write(6,*)

  atomloop : do 
     ! 
     ! name 
     !

     read(101,*,iostat=inpstat) ppfname
     if(scan(ppfname,'#')/=0) exit atomloop
     write(6,*)' Reading pseudopotential file: ',trim(ppfname); call flush(6)
     inquire(file='PP/'//trim(ppfname),exist=fex)
     if(.not.fex) then
        write(6,*) "ERROR: FILE: ",trim(ppfname)," is missing!"; stop
     end if
     
     !
     ! convert, and put in v1bTM and rrTM
     !
     call pp_convert(ppfname,ich,ival)
     
     !
     ! check and assign charges
     !
     call check_charge(ich)
     call assign_valence(ich,ival)
  end do atomloop
  call check_all_charges_assigned_pps
end subroutine pp_read_in

subroutine check_charge(ich)
  use atoms, only : found_charge,n118
  implicit none
  integer ich

  call check_lele(1,ich,n118,' 1-ich-118; charge should be 1 to 118 ')

  if(found_charge(ich)) then
     write(6,*)' ERROR: 2 or more PP for same Z,',ich
     stop
  endif
  found_charge(ich)=.true.
end subroutine check_charge

subroutine assign_valence(ich,ival)
  use atoms, only : valence_list
  implicit none
  integer ich, ival
  call check_lele(1,ival,ich,' 1-ival-Zch; valence should be < Z ')
  valence_list(ich)=ival
end subroutine assign_valence

subroutine check_all_charges_assigned_pps
  use atoms, only : matop, mapat, found_charge
  implicit none
  integer ma,ic
  do ma=1,matop
     ic= mapat(ma)
     if(.not.found_charge(ic)) then
        write(6,*)' ERROR: no PPs found for atom with charge ',ic
     end if
  end do
end subroutine check_all_charges_assigned_pps
