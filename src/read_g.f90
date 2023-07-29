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
subroutine read_g  ! grid, dens, possibly exchange too
  use gwm
  implicit none
  integer ip

  !
  !  read the orbitals if (i) trace not used; and  (ii.a) not read in wf_det nor (ii.b) stochastic-dft used  
  !
  ! now either reading wf_det and making density, or just reading density

  if (use_unified_input) then
     call make_chrg_from_dens
  else
     if (det_dft) then ! deterministic orbitals
        ip = 440
        call open_wf_file(ip)
        call read_xyz_wf(ip)
        call cnt_check_shift
        call alloc_dens0
        call read_wf_det_calc_dens0(ip)
        call close_wf_file(ip)
     else
        ip = 900
        call open_dens_file(ip)
        call read_xyz(ip)
        call cnt_check_shift
        call alloc_dens0
        call read_dens(ip)
        call close_dens_file(ip)
        call make_chrg_from_dens
     endif
  endif

end subroutine read_g

subroutine alloc_dens0
  use gwm, only : n, dens0,nsp
  implicit none
  integer st
  call check_le(1,n,' 1 n ')
  call check_lele(1,nsp,2,' 1-nsp-2 ')
  allocate(dens0(n,nsp), stat=st)
  call check0(st,' dens0 ')
end subroutine alloc_dens0

