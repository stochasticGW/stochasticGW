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
module gwm
  implicit none; save
  integer,    parameter :: lmx=3  ! 2 or 3
  integer,    parameter :: nvrl=1000
  integer,    parameter :: ipb = 40000
  real*8,     parameter :: angstrom_in_au = 1.8897161646d0
  real*8,     parameter :: tollocc = 1d-9
  real*8,     parameter :: tp_dns  = 1d-5
  integer,    parameter :: lngth_char=100
  integer               :: nproj_max = 20 ! increase if needed.  Modify to determine by reading.
  character*30          :: method_exch
  character(lngth_char) :: inputfname='INPUT'
  character(lngth_char) :: work_dir    
  character(lngth_char) :: work_dir_prmtv
  character(lngth_char) :: output_dir    
  character(lngth_char) :: scratch_dir
  character(lngth_char) :: scratch_path
  character(lngth_char) :: units_nuc

  integer :: funct_x,funct_c
  integer :: icounter
  integer :: color_size_pt
  integer :: ncolors, icolor
  integer :: nmctot 
  integer :: ngam   
  integer :: ngam_blk   
  integer :: nvr
  integer :: ns     
  integer :: ns_blk     
  integer :: ne     
  integer :: orb_kind
  integer :: orb_indx
  integer :: orb_det_kind
  integer :: nmc   
  integer :: scale_vh
  integer :: dim_periodic ! 0, 2, 3
  integer :: na
  integer :: nx=-1
  integer :: ny=-1
  integer :: nz=-1
  integer :: n
  integer :: nc_x
  integer :: nt
  integer :: nw
  integer :: nct
  integer :: mamx
  integer :: imc
  integer :: imctot
  integer :: imax_nl
  integer :: ncE
  integer :: flg_units_nuc
  integer :: xc_type
  integer :: ptheader_rank
  integer :: eta_rank
  integer :: nsp
  integer :: nspv 
  integer :: sp0
  integer :: ntop
  integer :: gamflg
  integer :: seg
  integer :: nchbmx=-1
  integer :: nchbc
  integer :: ngam_nzero_blk
  integer :: nj

  real*8 :: dv
  real*8 :: dx=-1d0
  real*8 :: dy=-1d0
  real*8 :: dz=-1d0
  real*8 :: dw
  real*8 :: ekcut
  real*8 :: normeta
  real*8 :: normxi
  real*8 :: rnel_orig
  real*8 :: rnel_neutral
  real*8 :: rnel  ! real # of electrons, not necc. integer.
  real*8 :: vk_scale_cnst=-1d0 ! normalization factor of vk.  Needs to be reset.
  real*8 :: dt    
  real*8 :: gama 
  real*8 :: sm    
  real*8 :: toll  
  real*8 :: eorb 
  real*8 :: de0
  real*8 :: nrm_trace
  real*8 :: chrg_net
  real*8 :: exchange_value
  real*8 :: mxmem_vo
  real*8 :: seg_fctr
  real*8 :: seg_frctn
  real*8 :: hmin
  real*8 :: hmax
  real*8 :: dhscl
  real*8 :: dh
  real*8 :: havg
  real*8 :: Tp
  real*8 :: mu
  real*8 :: homo
  real*8 :: lumo
  real*8, external :: filter_trace
  real*8, external :: wgf

  logical              :: trace
  logical              :: det_dft
  logical              :: det_tddft
  logical              :: print_now 
  logical              :: flgdyn    !true for chi_tddft, false for chi_hartree.
  logical              :: ptheader_flg
  logical              :: eta_flg
  logical              :: binwf 
  logical              :: read_eorb
  logical              :: stoch_x
  logical              :: rdexch         = .false.
  logical              :: rdscratch_path = .false.
  logical              :: flg_rd_orb     = .false. ! overwritten automatically if reading orbital.
  logical              :: boxpos
  logical              :: fuzzy_vk
  logical              :: rddh           = .false.
  logical              :: rdhavg         = .false.
  logical              :: filter_cheby   = .false.
  logical              :: rdmu           = .false.
  logical              :: gapped_flg     = .true.
  logical              :: usegpu
  logical              :: block_gam_alg
  logical              :: read_homo_input = .false.
  logical              :: read_lumo_input = .false.
  logical              :: read_jidx_input = .false.
  logical              :: use_unified_input = .false.

  integer, allocatable :: mapai(:)
  !integer, allocatable :: mapkbg(:,:)
  integer, allocatable :: is_map(:,:)
  integer, allocatable :: igam_map(:,:)
  integer, allocatable :: is_start(:)
  integer, allocatable :: is_end(:)
  integer, allocatable :: ns_each(:)
  integer, allocatable :: map_sa(:)
  integer, allocatable :: map_sp_det(:)
  integer, allocatable :: map_sp_pt(:)
  integer, allocatable :: map_x(:)
  integer, allocatable :: seg_w(:)
  integer, allocatable :: seg_sh(:)
  integer, allocatable :: orbj_indx(:)

  real*4, allocatable :: gam(:,:)    ! note real*4

  real*8, allocatable :: vk(:)   ! on big grid
  real*8, allocatable :: vk_exch(:) ! for periodic systems
  real*8, allocatable :: vk_old(:)  ! for periodic systems
  real*8, allocatable :: vks(:,:)
  real*8, allocatable :: vkst(:,:)
  real*8, allocatable :: dens0(:,:)
  real*8, allocatable :: vxc(:,:)
  real*8, allocatable :: vxc0(:,:)
  real*8, allocatable :: gamvr(:,:)
  real*8, allocatable :: ge(:,:)
  real*8, allocatable :: g(:)
  real*8, allocatable :: zeta(:)
  real*8, allocatable :: den(:)
  real*8, allocatable :: del(:)
  real*8, allocatable :: war(:)
  real*8, allocatable :: ear(:)
  real*8, allocatable :: ekn(:)
  !real*8, allocatable :: fr(:,:)
  real*8, allocatable :: th_co(:)
  real*8, allocatable :: gaus_co(:)
  real*8, allocatable :: normpt(:)
  real*8, allocatable :: exce(:), wge(:), vxce(:)
  real*8, allocatable :: phi_det(:,:)! det.orbitals, used only if det_dft
  real*8, allocatable :: occ_det(:)  !occ.det.orbts, used only if det_dt
  real*8, allocatable :: occ_det_full(:)
  real*8, allocatable :: evl(:)
  real*8, allocatable :: orb_rd(:)
  real*8, allocatable :: vloc_tot(:)
  real*8, allocatable :: gge(:)
  real*8, allocatable :: vnl_evl(:,:,:)
  real*8, allocatable :: vkb_many(:,:,:,:)
  real*8, allocatable :: rho_core(:)
  real*8, allocatable :: rho(:,:)
  real*8, allocatable :: rho_p(:)
  real*8, allocatable :: vh0(:)
  real*8, allocatable :: vt(:,:)
  real*8, allocatable :: phi_x(:,:)
  real*8, allocatable :: gej(:,:,:)
  real*8, allocatable :: wgej(:,:), vxcej(:,:)
  real*8, allocatable :: rdorbj(:,:)
  real*8, allocatable :: eorbj(:)

  complex*16, allocatable :: vo(:,:)
  complex*16, allocatable :: xi(:)
  complex*16, allocatable :: eta(:)
  complex*16, allocatable :: ct(:,:)
  complex*16, allocatable :: pt(:,:)
  complex*16, allocatable :: expnk(:)
  complex*16, allocatable :: c_stoch_weight(:,:)
  complex*16, allocatable :: sg(:)
  complex*16, allocatable :: sq(:)
  complex*16, allocatable :: expvks0(:)
  complex*16, allocatable :: cgEco(:,:)
  complex*16, allocatable :: etaxi(:,:)
  complex*16, allocatable :: cveta(:)
  complex*16, allocatable :: cvxi(:)
  complex*8,  allocatable :: ft(:,:)
  complex*16, allocatable :: excej(:,:)
  complex*16, allocatable :: ctj(:,:,:)
  complex*16, allocatable :: cvetaj(:,:)
  complex*16, allocatable :: cvxij(:,:)

end module gwm
