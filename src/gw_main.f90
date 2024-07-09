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
program gw
  call main
end program gw

subroutine main
  use simple_mpi
  use gwm
#if GPU_ENABLED
  use device_mem_module
#endif
  implicit none
  real*8 :: tstart,tend
  call prepare_mpi_lib
  call prep_gw
  call run_gw
  call dealloc_gw
  call prnt_time_mpi(" program end ",tstart)
  call finalize_mpi_lib
contains
  subroutine prep_gw
    implicit none
    call cpu_time(tstart)
    call write_header
    call check_input_name
    call counter_read

    call work_dir_make
    call scratch_dir_name_read
    call scratch_dir_make
    call output_dir_set
    call output_dir_make

    call set_details_file_name
    call make_atoms_list
    call read_units
    call read_input_vars
    call read_unified
    call read_cnt
    call map_cnt
    call prep_rand
    call nrppmx_read
    call pp_prep 
    call bcast_valence_atoms
    call make_chrg_for_wf
    call ns_blk_make
    call color  
    call modify_nmctot
    call prnt_settings 
    call set_gamatw
    call xc_alloc
    call read_orb
    call vk_prep
    call vk_set_exch
    call vk_possibly_reset
    call read_g    ! dns, possibly exchange
    call prep_map_sp_pt
    call alloc_prep
    call ek_prep
    call make_kb
    call make_rho_core
    call dens_to_v_spn(dens0,vks,n, nsp)
    call prep_vxc
    call nvr_prep
    call seg_prep
    call set_ngam_nzero_blk
    if(.not.trace)   call calc_eorb
#if GPU_ENABLED
    if (usegpu) call init_device
#endif
    if(filter_cheby .and. .not.rdhavg)  call improve_hmin_hmax
    if(filter_cheby) call set_theta
    !call exchange_set 
    call prnt_time_mpi(" preparation ",tstart)
  end subroutine prep_gw

  subroutine run_gw
    implicit none
    mcloop : do imc=1,nmc
       call check_prnt
       call write_mc_header
       call gw_core
       call write_mc_timing
    end do mcloop
  end subroutine run_gw

  subroutine write_mc_header
    use simple_mpi, only : nodes
    implicit none
    integer nds
    nds = max(nodes,1)
    imctot = imc*(nds/color_size)
    if(rank==0.and.print_now) then; 
       write(6,*)
       write(6,'(A,I10,A)')" ######## Accumulative # of Monte Carlo steps :  MC= ",imctot,&
                           " ######## "
       write(6,'(A,I10)')  " Completed Monte Carlo steps in each core : MC_STEP= ",imc  
       write(6,*)
       call flush(6)
    end if
  end subroutine write_mc_header
  
  subroutine write_mc_timing
    implicit none
    if(print_now)call prnt_time_i_mpi(" MC_STEP: ",imc,tstart)
  end subroutine write_mc_timing
end subroutine main


