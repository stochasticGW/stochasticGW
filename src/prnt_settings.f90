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
subroutine prnt_settings
  use gwm
  use simple_mpi, only : rank,nodes,color_size
  use ppm, only : nrppmx
  implicit none
  integer nds
  character(40)axc
  character(25) homotag,lumotag

  nds = max(nodes,1)

  select case(xc_type)
  case(1)
     axc=" read from vxc.txt"
  case(2)
#if LIBXC_ENABLED
     if (funct_x.eq.0 .and. funct_c.eq.0) then
        axc=" LDA (PW91) "
     else
        write(axc,'(X,A,2(X,I0),A)') ' (LIBXC#',funct_x,funct_c,')'
     endif
#else
     axc=" LDA (PW91) "
#endif
  case default
     stop " xc_type at present should be 1 or 2 "
  end select

  if(rank==0) then
     write(6,*) " ############# ECHO OF INPUT VARIABLES #############" 
     write(6,*) " MC (total)            =",nmctot
     write(6,*) " # MC_STEP (each core) =",nmctot/ncolors
     write(6,*) " multiplicity (=nsp)   =",nsp
     write(6,*) " work_dir:       = ",trim(work_dir)
     write(6,*) " scratch_dir     = ",trim(scratch_dir)
     write(6,*) " centers units   = ",trim(units_nuc)
     if (flgdyn) then
        write(6,*) " flgdyn          =",flgdyn,' i.e., W-TDDFT '
     else
        write(6,*) " flgdyn          =",flgdyn,' i.e., W-RPA   '
     endif
     write(6,*) " det_tddft       =",det_tddft
     write(6,*) " binary inp.     =",binwf
     write(6,*) " periodic        =",(dim_periodic.gt.0)
     write(6,*) " block_gam_alg   =",block_gam_alg
     write(6,*) " scale_vh        =",scale_vh
     write(6,*) " N_cores (mpi)   =",nds
     write(6,*) " buffer_size     =",color_size
     write(6,*) " ntddft*nsp      =",ns_blk
     write(6,*) " orb_indx        =",orb_indx
     write(6,*) " PP:nrppmx       =",nrppmx
     write(6,*) " xc_type         =",xc_type,axc
     write(6,*) " projection      =",.not.filter_cheby
     if (filter_cheby) then
        write(6,*) " gapped-filtering=",gapped_flg
        if (gapped_flg) then
           if (read_homo_input) then
              homotag=' (read from INPUT)'
           else
              homotag=' (read from "sgwinp.txt")'
           endif
           if (read_lumo_input) then
              lumotag=' (read from INPUT)'
           else
              lumotag=' (read from "sgwinp.txt")'
           endif
           write(6,*) " homo            =",homo,homotag
           write(6,*) " lumo            =",lumo,lumotag
        else
           write(6,*) " mu              =",mu
           write(6,*) " Tp              =",Tp 
        endif
#if GPU_ENABLED
    write(6,*) " hminmax on GPU  =",(.not.disable_gpu_hminmax)
    write(6,*) " filter  on GPU  =",(.not.disable_gpu_filter)
#endif
     endif
#if GPU_ENABLED
    write(6,*) " gam     on GPU  =",(.not.(disable_gpu_gam.or.use_host_curand))
    write(6,*) " prop    on GPU  =",(.not.disable_gpu_prop)
#endif
     write(6,'(1X,A,F12.5)') " ekcut           =",ekcut
     write(6,'(1X,A,F12.5)') " dt              =",dt
     write(6,'(1X,A,F12.5)') " gamma           =",gama
     write(6,'(1X,A,F12.5)') " sm              =",sm
     write(6,'(1X,A,F12.5)') " charge_net      =",chrg_net
     !write(6,*) " det_dft         =",det_dft
     !write(6,*) " trace           =",trace
     if(trace) then
     write(6,'(1X,A,F12.5)') " eorb (trace)    =",eorb
     write(6,'(1X,A,F12.5)') " de0  (trace)    =",de0
     endif
     select case(gamflg)
     case(1)
        write(6,*) " using old xi    ="
        write(6,'(A,I10)') " nxi             =",ngam_blk
     case(2)
        write(6,'(1X,A,I10)') " nxi             =",ngam_blk
        write(6,'(1X,A,F12.5)') " segment_fraction=",seg_frctn
     case(3)
        write(6,*) " direct vr->vo   ="
        write(6,*) " writing strp sz =",nvr
     case default
        stop " ERROR: gamflg problem  "
     end select
     write(6,*)
     call flush(6)
  end if !rank=0 only
end subroutine prnt_settings

