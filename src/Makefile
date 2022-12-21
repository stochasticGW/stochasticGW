## compiler options
FCMPI = mpifort
MPIFLG = -DMPI -O3
PREPROCFLG = -Mpreprocess

## Set below to use GPU
USEGPU = yes

# Libraries
FFTFLG = -L/usr/bin/fftw-3.3.9/lib -lfftw3
LIBXC = -I/home/filip/Documents/LBL/libxc/nv-5.2.2/include -L/home/filip/Documents/LBL/libxc/nv-5.2.2/lib -lxcf03 -lxc
#LIBXC = -I/usr/bin/LIBXC/include -L/usr/bin/LIBXC/lib -lxcf03 -lxc


# GPU libraries
CUFLG = -L/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/cuda/lib64 -lcudart
CULIB = -L/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/math_libs/lib64 -lcufft -lcublas

MAINFILES1 = \
         1_gwm.f90 \
         2_atoms.f90 \
         2_kb_top_module.f90 \
         3_ppm.f90

GPU_FILES = \
        7_device_mem_module.f90 \
        8_chebfilt.f90 \
        8_ct.f90 \
        8_filter.f90 \
        8_propdt_pt.f90 \
        8_vr.f90

MAINFILES2 = \
        allocate_calc_del_zeta.f90 \
        alloc_prep.f90 \
        atoms_bcast.f90 \
        atoms_list.f90 \
        atoms_valence_bcast.f90 \
        check_chrg_det.f90 \
        check_prnt.f90 \
        chrg_make.f90 \
        cnt_check_shift.f90 \
        cnt_map.f90 \
        cnt_read.f90 \
        color.f90 \
        color_flgs.f90 \
        color_indx.f90 \
        color_make.f90 \
        counter_read.f90 \
        dealloc_gw.f90 \
        debug_analyze_h.f90 \
        debug_mpi.f90 \
        ek.f90 \
        eorb.f90 \
        eorbj.f90 \
        eta_det_tddft.f90 \
        eta_fltr.f90 \
        eta_xi.f90 \
        exch_accum.f90 \
        exchange_corr_lda_spn.f90 \
        exchange_stoch.f90 \
        exch_sum.f90 \
        fft.f90 \
        fft_fff.f90 \
        fftsa.f90 \
        fftw_good.f90 \
        filter_chebyshev.f90 \
        filter_trace.f90 \
        fr1_gama.f90 \
        gam_dot.f90 \
        gam_prep.f90 \
        gam_seg1.f90 \
        gam_seg.f90 \
        gamvr.f90 \
        get_ma.f90 \
        gw_core.f90 \
        gw_main.f90 \
        gwplt.f90 \
        header.f90 \
        hmin_hmax.f90 \
        hmin_max_improve.f90 \
        hpsi_addvnl.f90 \
        hscld_psi.f90 \
        h_psi.f90 \
        ip_make.f90 \
        ip_open.f90 \
        make_rho.f90 \
        makevt.f90 \
        map_sp_pt.f90 \
        mbt_set.f90 \
        misc.f90 \
        modify_input.f90 \
        name_check_input.f90 \
        name_file_set.f90 \
        nrm_trace.f90 \
        nrppmx_read.f90 \
        ns_blk_make.f90 \
        nvr_prep.f90 \
        open_close_wf.f90 \
        openfile.f90 \
        output_dir_make.f90 \
        output_dir_set.f90 \
        pp_allocate.f90 \
        pp_bcast.f90 \
        pp_check.f90 \
        pp_convert.f90 \
        pp_prep.f90 \
        pp_read.f90 \
        pp_read_in.f90 \
        prnt.f90 \
        prnt_req.f90 \
        prnt_settings.f90 \
        prnt_time_mpi.f90 \
        propdt.f90 \
        prop_pert.f90 \
        pt_det_tddft.f90 \
        pt_eta_effic.f90 \
        pt_eta.f90 \
        pt_fltr.f90 \
        pt_stoch.f90 \
        rand.f90 \
        read_dens.f90 \
        read_det_prep.f90 \
        read_evlocc.f90 \
        read_g.f90 \
        read_input_vars.f90 \
        read_orb_1.f90 \
        read_orb.f90 \
        read_orb_fromwf.f90 \
        read_orbj_1.f90 \
        read_orbj_fromwf.f90 \
        read_unified.f90 \
        read_vxc.f90 \
        read_wf1.f90 \
        read_wf.f90 \
        read_xyz.f90 \
        rho_core.f90 \
        rho_det.f90 \
        rw_ge.f90 \
        rw_pt.f90 \
        scratch_dir_make.f90 \
        scratch_dir_name_read.f90 \
        scratch_dir_reset.f90 \
        seed_gw_modu.f90 \
        seed_set.f90 \
        seed_set_ran_ps.f90 \
        seg_prep.f90 \
        set_ge.f90 \
        set_ngam_nzero_blk.f90 \
        setocc.f90 \
        start_prt_version.f90 \
        theta.f90 \
        t_to_w.f90 \
        t_to_w_set_gamma.f90 \
        units_assign.f90 \
        units_read.f90 \
        vd_wg.f90 \
        vfuzzy_mc_2d.f90 \
        vfuzzy_mc_3d.f90 \
        vh_big_c.f90 \
        vh_big.f90 \
        vh_sub_exch.f90 \
        vk_prep.f90 \
        vh_sum.f90 \
        vk_scale.f90 \
        vks_spn.f90 \
        vk_2d.f90 \
        vk_3d.f90 \
        vk_possibly_reset.f90 \
        vk_set_exch.f90 \
        vk_TuMa.f90 \
        vo1_drct.f90 \
        vo1.f90 \
        vo1_gam.f90 \
        vo.f90 \
        vos.f90 \
        vr1_prnt.f90 \
        vr1_process.f90 \
        vr_close.f90 \
        vr.f90 \
        vr_open.f90 \
        vr_to_vo_drct.f90 \
        vr_to_vo_gam.f90 \
        vs_read_rmv.f90 \
        vxc_libxc.f90 \
        vxc_pbe.f90 \
        vxc_prep.f90 \
        vxc_spn.f90 \
        wgf.f90 \
        work_dir_make.f90 \
        work_dir_reset.f90 \
        xc_alloc.f90 \
        xc_expect.f90

# Process flags
USEGPU := $(strip ${USEGPU})
SRC_FILES = ${MAINFILES1}
GPULIB =
ifeq (${USEGPU},yes)
   GPULIB += $(CUFLG) $(CULIB)
   FCMPI += -cuda -acc=gpu -fast -Minfo=accel -DGPU_ENABLED=1
   SRC_FILES += ${GPU_FILES}
else
   FCMPI += -DGPU_ENABLED=0
endif
SRC_FILES += ${MAINFILES2}
FCMPI += ${PREPROCFLG}

# compile
all :  clean sgw clnomod 

sgw: simmpi gwm atoms kbtop ppm kbmod
	$(FCMPI) $(MPIFLG) -o sgw.x kb/*f90 libgw/*f libgw/*f90 $(SRC_FILES) XCI/*.f90 $(LIBXC) $(FFTFLG) $(GPULIB)

simmpi:
	$(FCMPI) $(MPIFLG) -c ./libgw/0_library_mpi_module.f90 ./libgw/0_simple_mpi.f90

gwm:
	$(FCMPI) $(MPIFLG) -c 1_gwm.f90

atoms:
	$(FCMPI) $(MPIFLG) -c 2_atoms.f90

kbtop:
	$(FCMPI) $(MPIFLG) -c 2_kb_top_module.f90

ppm:
	$(FCMPI) $(MPIFLG) -c 3_ppm.f90

kbmod:
	$(FCMPI) $(MPIFLG) -c kb/1_kb_mod.f90

clnomod:
	rm -rf *o *mod kb/*o libgw/*o 

clean:
	rm -rf *o *mod kb/*o libgw/*o sgw.x 

