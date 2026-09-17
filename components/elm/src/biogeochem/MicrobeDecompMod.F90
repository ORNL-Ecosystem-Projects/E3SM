module MicrobeDecompMod

  !-----------------------------------------------------------------------
  ! Microbial CTC decomposition cascade.
  !
  ! This module contains only the DOM/bacteria/fungi extension of the ELM
  ! decomposition cascade. Revised methane reactions are intentionally kept
  ! in MicrobeMethaneMod and are not part of this Phase 2 implementation.
  !-----------------------------------------------------------------------

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use abortutils             , only : endrun
  use decompMod              , only : bounds_type
  use elm_varpar             , only : nlevdecomp, i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use elm_varpar             , only : i_bacteria, i_fungi, i_dom
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType            , only : cnstate_type
  use ColumnType             , only : col_pp
  use VegetationType         , only : veg_pp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

  implicit none
  save
  private

  integer, parameter, public :: microbe_decomp_npools = 11
  integer, parameter, public :: microbe_decomp_ntransitions = 47

  type, public :: MicrobeDecompParamsType
     ! PFT-indexed daily turnover probabilities and pathway parameters.
     real(r8), allocatable :: k_dom(:)
     real(r8), allocatable :: k_bacteria(:)
     real(r8), allocatable :: k_fungi(:)
     real(r8), allocatable :: m_rf_s1m(:)
     real(r8), allocatable :: m_rf_s2m(:)
     real(r8), allocatable :: m_rf_s3m(:)
     real(r8), allocatable :: m_rf_s4m(:)
     real(r8), allocatable :: m_batm_f(:)
     real(r8), allocatable :: m_fatm_f(:)
     real(r8), allocatable :: m_bdom_f(:)
     real(r8), allocatable :: m_fdom_f(:)
     real(r8), allocatable :: m_bs1_f(:)
     real(r8), allocatable :: m_bs2_f(:)
     real(r8), allocatable :: m_bs3_f(:)
     real(r8), allocatable :: m_fs1_f(:)
     real(r8), allocatable :: m_fs2_f(:)
     real(r8), allocatable :: m_fs3_f(:)
     real(r8), allocatable :: m_domb_f(:)
     real(r8), allocatable :: m_domf_f(:)
     real(r8), allocatable :: m_doms1_f(:)
     real(r8), allocatable :: m_doms2_f(:)
     real(r8), allocatable :: m_doms3_f(:)
     real(r8), allocatable :: cn_bacteria(:)
     real(r8), allocatable :: cn_fungi(:)

     ! Fixed pool stoichiometry and cascade controls.
     real(r8) :: cn_dom
     real(r8) :: bacteria_pool_cn
     real(r8) :: fungi_pool_cn
     real(r8) :: cp_bacteria
     real(r8) :: cp_fungi
     real(r8) :: cp_dom
     real(r8) :: cue_max
     real(r8) :: cue_cn_target
     real(r8) :: allocation_cn_exponent
     real(r8) :: bacteria_initial_c
     real(r8) :: fungi_initial_c
     real(r8) :: dom_initial_c
     real(r8) :: dom_som_diffusion_multiplier
     real(r8) :: microbe_som2_q10
     real(r8) :: microbe_som3_q10
     real(r8) :: microbe_som4_q10
     real(r8) :: microbe_dom_q10

     ! Litter/SOM solubilization and direct stabilization fractions.
     real(r8) :: l1dom_f
     real(r8) :: l2dom_f
     real(r8) :: l3dom_f
     real(r8) :: s1dom_f
     real(r8) :: s2dom_f
     real(r8) :: s3dom_f
     real(r8) :: s4dom_f
     real(r8) :: l1s1_f
     real(r8) :: l2s2_f
     real(r8) :: l3s3_f
     real(r8) :: s1s2_f
     real(r8) :: s2s3_f
     real(r8) :: s3s4_f
  end type MicrobeDecompParamsType

  type(MicrobeDecompParamsType), public :: MicrobeDecompParamsInst

  public :: readMicrobeDecompParams
  public :: initMicrobeDecompCascade
  public :: getMicrobeDecompTimestepRates
  public :: validateMicrobeDecompDonorClosure

contains

  !-----------------------------------------------------------------------
  subroutine readMicrobeDecompParams(ncid)
    use ncdio_pio, only : file_desc_t, ncd_io, ncd_inqdid, ncd_inqdlen

    type(file_desc_t), intent(inout) :: ncid

    integer :: dimid
    integer :: npft

    call ncd_inqdid(ncid, 'pft', dimid)
    call ncd_inqdlen(ncid, dimid, npft)

    call read_pft_parameter(ncid, 'k_dom', MicrobeDecompParamsInst%k_dom, npft)
    call read_pft_parameter(ncid, 'k_bacteria', MicrobeDecompParamsInst%k_bacteria, npft)
    call read_pft_parameter(ncid, 'k_fungi', MicrobeDecompParamsInst%k_fungi, npft)
    call read_pft_parameter(ncid, 'm_rf_s1m', MicrobeDecompParamsInst%m_rf_s1m, npft)
    call read_pft_parameter(ncid, 'm_rf_s2m', MicrobeDecompParamsInst%m_rf_s2m, npft)
    call read_pft_parameter(ncid, 'm_rf_s3m', MicrobeDecompParamsInst%m_rf_s3m, npft)
    call read_pft_parameter(ncid, 'm_rf_s4m', MicrobeDecompParamsInst%m_rf_s4m, npft)
    call read_pft_parameter(ncid, 'm_batm_f', MicrobeDecompParamsInst%m_batm_f, npft)
    call read_pft_parameter(ncid, 'm_fatm_f', MicrobeDecompParamsInst%m_fatm_f, npft)
    call read_pft_parameter(ncid, 'm_bdom_f', MicrobeDecompParamsInst%m_bdom_f, npft)
    call read_pft_parameter(ncid, 'm_fdom_f', MicrobeDecompParamsInst%m_fdom_f, npft)
    call read_pft_parameter(ncid, 'm_bs1_f', MicrobeDecompParamsInst%m_bs1_f, npft)
    call read_pft_parameter(ncid, 'm_bs2_f', MicrobeDecompParamsInst%m_bs2_f, npft)
    call read_pft_parameter(ncid, 'm_bs3_f', MicrobeDecompParamsInst%m_bs3_f, npft)
    call read_pft_parameter(ncid, 'm_fs1_f', MicrobeDecompParamsInst%m_fs1_f, npft)
    call read_pft_parameter(ncid, 'm_fs2_f', MicrobeDecompParamsInst%m_fs2_f, npft)
    call read_pft_parameter(ncid, 'm_fs3_f', MicrobeDecompParamsInst%m_fs3_f, npft)
    call read_pft_parameter(ncid, 'm_domb_f', MicrobeDecompParamsInst%m_domb_f, npft)
    call read_pft_parameter(ncid, 'm_domf_f', MicrobeDecompParamsInst%m_domf_f, npft)
    call read_pft_parameter(ncid, 'm_doms1_f', MicrobeDecompParamsInst%m_doms1_f, npft)
    call read_pft_parameter(ncid, 'm_doms2_f', MicrobeDecompParamsInst%m_doms2_f, npft)
    call read_pft_parameter(ncid, 'm_doms3_f', MicrobeDecompParamsInst%m_doms3_f, npft)
    call read_pft_parameter(ncid, 'cn_bacteria', MicrobeDecompParamsInst%cn_bacteria, npft)
    call read_pft_parameter(ncid, 'cn_fungi', MicrobeDecompParamsInst%cn_fungi, npft)

    call read_scalar_parameter(ncid, 'cn_dom', MicrobeDecompParamsInst%cn_dom)
    call read_scalar_parameter(ncid, 'bacteria_pool_cn', MicrobeDecompParamsInst%bacteria_pool_cn)
    call read_scalar_parameter(ncid, 'fungi_pool_cn', MicrobeDecompParamsInst%fungi_pool_cn)
    call read_scalar_parameter(ncid, 'cp_bacteria', MicrobeDecompParamsInst%cp_bacteria)
    call read_scalar_parameter(ncid, 'cp_fungi', MicrobeDecompParamsInst%cp_fungi)
    call read_scalar_parameter(ncid, 'cp_dom', MicrobeDecompParamsInst%cp_dom)
    call read_scalar_parameter(ncid, 'CUEmax', MicrobeDecompParamsInst%cue_max)
    call read_scalar_parameter(ncid, 'microbe_cue_cn_target', MicrobeDecompParamsInst%cue_cn_target)
    call read_scalar_parameter(ncid, 'microbe_allocation_cn_exponent', &
         MicrobeDecompParamsInst%allocation_cn_exponent)
    call read_scalar_parameter(ncid, 'bacteria_initial_c', MicrobeDecompParamsInst%bacteria_initial_c)
    call read_scalar_parameter(ncid, 'fungi_initial_c', MicrobeDecompParamsInst%fungi_initial_c)
    call read_scalar_parameter(ncid, 'dom_initial_c', MicrobeDecompParamsInst%dom_initial_c)
    call read_scalar_parameter(ncid, 'dom_som_diffusion_multiplier', &
         MicrobeDecompParamsInst%dom_som_diffusion_multiplier)
    call read_scalar_parameter(ncid, 'microbe_som2_q10', MicrobeDecompParamsInst%microbe_som2_q10)
    call read_scalar_parameter(ncid, 'microbe_som3_q10', MicrobeDecompParamsInst%microbe_som3_q10)
    call read_scalar_parameter(ncid, 'microbe_som4_q10', MicrobeDecompParamsInst%microbe_som4_q10)
    call read_scalar_parameter(ncid, 'microbe_dom_q10', MicrobeDecompParamsInst%microbe_dom_q10)

    call read_scalar_parameter(ncid, 'l1dom_f', MicrobeDecompParamsInst%l1dom_f)
    call read_scalar_parameter(ncid, 'l2dom_f', MicrobeDecompParamsInst%l2dom_f)
    call read_scalar_parameter(ncid, 'l3dom_f', MicrobeDecompParamsInst%l3dom_f)
    call read_scalar_parameter(ncid, 's1dom_f', MicrobeDecompParamsInst%s1dom_f)
    call read_scalar_parameter(ncid, 's2dom_f', MicrobeDecompParamsInst%s2dom_f)
    call read_scalar_parameter(ncid, 's3dom_f', MicrobeDecompParamsInst%s3dom_f)
    call read_scalar_parameter(ncid, 's4dom_f', MicrobeDecompParamsInst%s4dom_f)
    call read_scalar_parameter(ncid, 'l1s1_f', MicrobeDecompParamsInst%l1s1_f)
    call read_scalar_parameter(ncid, 'l2s2_f', MicrobeDecompParamsInst%l2s2_f)
    call read_scalar_parameter(ncid, 'l3s3_f', MicrobeDecompParamsInst%l3s3_f)
    call read_scalar_parameter(ncid, 's1s2_f', MicrobeDecompParamsInst%s1s2_f)
    call read_scalar_parameter(ncid, 's2s3_f', MicrobeDecompParamsInst%s2s3_f)
    call read_scalar_parameter(ncid, 's3s4_f', MicrobeDecompParamsInst%s3s4_f)

    call validate_parameters()

  contains

    subroutine read_pft_parameter(file, name, values, count)
      type(file_desc_t), intent(inout) :: file
      character(len=*), intent(in) :: name
      real(r8), allocatable, intent(out) :: values(:)
      integer, intent(in) :: count

      logical :: readv

      allocate(values(0:count-1))
      call ncd_io(trim(name), values, 'read', file, readvar=readv, posNOTonfile=.true.)
      if (.not. readv) then
         call endrun(msg=' ERROR: microbial decomposition parameter is missing: '//trim(name)//&
              errMsg(__FILE__, __LINE__))
      end if
    end subroutine read_pft_parameter

    subroutine read_scalar_parameter(file, name, value)
      type(file_desc_t), intent(inout) :: file
      character(len=*), intent(in) :: name
      real(r8), intent(out) :: value

      logical :: readv

      call ncd_io(trim(name), value, 'read', file, readvar=readv)
      if (.not. readv) then
         call endrun(msg=' ERROR: microbial decomposition parameter is missing: '//trim(name)//&
              errMsg(__FILE__, __LINE__))
      end if
    end subroutine read_scalar_parameter

  end subroutine readMicrobeDecompParams

  !-----------------------------------------------------------------------
  subroutine validate_parameters()
    integer :: pft
    real(r8) :: residual

    call require_probability_array('k_dom', MicrobeDecompParamsInst%k_dom, .false.)
    call require_probability_array('k_bacteria', MicrobeDecompParamsInst%k_bacteria, .false.)
    call require_probability_array('k_fungi', MicrobeDecompParamsInst%k_fungi, .false.)
    call require_fraction_array('m_rf_s1m', MicrobeDecompParamsInst%m_rf_s1m)
    call require_fraction_array('m_rf_s2m', MicrobeDecompParamsInst%m_rf_s2m)
    call require_fraction_array('m_rf_s3m', MicrobeDecompParamsInst%m_rf_s3m)
    call require_fraction_array('m_rf_s4m', MicrobeDecompParamsInst%m_rf_s4m)
    call require_fraction_array('m_batm_f', MicrobeDecompParamsInst%m_batm_f)
    call require_fraction_array('m_fatm_f', MicrobeDecompParamsInst%m_fatm_f)
    call require_fraction_array('m_bdom_f', MicrobeDecompParamsInst%m_bdom_f)
    call require_fraction_array('m_fdom_f', MicrobeDecompParamsInst%m_fdom_f)
    call require_fraction_array('m_bs1_f', MicrobeDecompParamsInst%m_bs1_f)
    call require_fraction_array('m_bs2_f', MicrobeDecompParamsInst%m_bs2_f)
    call require_fraction_array('m_bs3_f', MicrobeDecompParamsInst%m_bs3_f)
    call require_fraction_array('m_fs1_f', MicrobeDecompParamsInst%m_fs1_f)
    call require_fraction_array('m_fs2_f', MicrobeDecompParamsInst%m_fs2_f)
    call require_fraction_array('m_fs3_f', MicrobeDecompParamsInst%m_fs3_f)
    call require_fraction_array('m_domb_f', MicrobeDecompParamsInst%m_domb_f)
    call require_fraction_array('m_domf_f', MicrobeDecompParamsInst%m_domf_f)
    call require_fraction_array('m_doms1_f', MicrobeDecompParamsInst%m_doms1_f)
    call require_fraction_array('m_doms2_f', MicrobeDecompParamsInst%m_doms2_f)
    call require_fraction_array('m_doms3_f', MicrobeDecompParamsInst%m_doms3_f)
    call require_positive_array('cn_bacteria', MicrobeDecompParamsInst%cn_bacteria)
    call require_positive_array('cn_fungi', MicrobeDecompParamsInst%cn_fungi)

    call require_positive('cn_dom', MicrobeDecompParamsInst%cn_dom)
    call require_positive('bacteria_pool_cn', MicrobeDecompParamsInst%bacteria_pool_cn)
    call require_positive('fungi_pool_cn', MicrobeDecompParamsInst%fungi_pool_cn)
    call require_positive('cp_bacteria', MicrobeDecompParamsInst%cp_bacteria)
    call require_positive('cp_fungi', MicrobeDecompParamsInst%cp_fungi)
    call require_positive('cp_dom', MicrobeDecompParamsInst%cp_dom)
    call require_fraction('CUEmax', MicrobeDecompParamsInst%cue_max)
    call require_positive('microbe_cue_cn_target', MicrobeDecompParamsInst%cue_cn_target)
    call require_positive('microbe_allocation_cn_exponent', &
         MicrobeDecompParamsInst%allocation_cn_exponent)
    call require_nonnegative('bacteria_initial_c', MicrobeDecompParamsInst%bacteria_initial_c)
    call require_nonnegative('fungi_initial_c', MicrobeDecompParamsInst%fungi_initial_c)
    call require_nonnegative('dom_initial_c', MicrobeDecompParamsInst%dom_initial_c)
    call require_positive('dom_som_diffusion_multiplier', &
         MicrobeDecompParamsInst%dom_som_diffusion_multiplier)
    call require_positive('microbe_som2_q10', MicrobeDecompParamsInst%microbe_som2_q10)
    call require_positive('microbe_som3_q10', MicrobeDecompParamsInst%microbe_som3_q10)
    call require_positive('microbe_som4_q10', MicrobeDecompParamsInst%microbe_som4_q10)
    call require_positive('microbe_dom_q10', MicrobeDecompParamsInst%microbe_dom_q10)

    call require_fraction('l1dom_f', MicrobeDecompParamsInst%l1dom_f)
    call require_fraction('l2dom_f', MicrobeDecompParamsInst%l2dom_f)
    call require_fraction('l3dom_f', MicrobeDecompParamsInst%l3dom_f)
    call require_fraction('s1dom_f', MicrobeDecompParamsInst%s1dom_f)
    call require_fraction('s2dom_f', MicrobeDecompParamsInst%s2dom_f)
    call require_fraction('s3dom_f', MicrobeDecompParamsInst%s3dom_f)
    call require_fraction('s4dom_f', MicrobeDecompParamsInst%s4dom_f)
    call require_fraction('l1s1_f', MicrobeDecompParamsInst%l1s1_f)
    call require_fraction('l2s2_f', MicrobeDecompParamsInst%l2s2_f)
    call require_fraction('l3s3_f', MicrobeDecompParamsInst%l3s3_f)
    call require_fraction('s1s2_f', MicrobeDecompParamsInst%s1s2_f)
    call require_fraction('s2s3_f', MicrobeDecompParamsInst%s2s3_f)
    call require_fraction('s3s4_f', MicrobeDecompParamsInst%s3s4_f)

    call require_nonnegative_residual('litter1', 1._r8 - MicrobeDecompParamsInst%l1dom_f - &
         MicrobeDecompParamsInst%l1s1_f)
    call require_nonnegative_residual('litter2', 1._r8 - MicrobeDecompParamsInst%l2dom_f - &
         MicrobeDecompParamsInst%l2s2_f)
    call require_nonnegative_residual('litter3', 1._r8 - MicrobeDecompParamsInst%l3dom_f - &
         MicrobeDecompParamsInst%l3s3_f)
    call require_nonnegative_residual('soil1', 1._r8 - MicrobeDecompParamsInst%s1dom_f - &
         MicrobeDecompParamsInst%s1s2_f)
    call require_nonnegative_residual('soil2', 1._r8 - MicrobeDecompParamsInst%s2dom_f - &
         MicrobeDecompParamsInst%s2s3_f)
    call require_nonnegative_residual('soil3', 1._r8 - MicrobeDecompParamsInst%s3dom_f - &
         MicrobeDecompParamsInst%s3s4_f)

    do pft = lbound(MicrobeDecompParamsInst%k_dom, 1), &
         ubound(MicrobeDecompParamsInst%k_dom, 1)
       residual = 1._r8 - MicrobeDecompParamsInst%m_batm_f(pft) - &
            MicrobeDecompParamsInst%m_bdom_f(pft) - MicrobeDecompParamsInst%m_bs1_f(pft) - &
            MicrobeDecompParamsInst%m_bs2_f(pft) - MicrobeDecompParamsInst%m_bs3_f(pft)
       call require_nonnegative_residual('bacteria', residual, pft)

       residual = 1._r8 - MicrobeDecompParamsInst%m_fatm_f(pft) - &
            MicrobeDecompParamsInst%m_fdom_f(pft) - MicrobeDecompParamsInst%m_fs1_f(pft) - &
            MicrobeDecompParamsInst%m_fs2_f(pft) - MicrobeDecompParamsInst%m_fs3_f(pft)
       call require_nonnegative_residual('fungi', residual, pft)

       residual = 1._r8 - MicrobeDecompParamsInst%m_domb_f(pft) - &
            MicrobeDecompParamsInst%m_domf_f(pft) - MicrobeDecompParamsInst%m_doms1_f(pft) - &
            MicrobeDecompParamsInst%m_doms2_f(pft) - MicrobeDecompParamsInst%m_doms3_f(pft)
       call require_nonnegative_residual('DOM', residual, pft)
    end do

  end subroutine validate_parameters

  !-----------------------------------------------------------------------
  subroutine initMicrobeDecompCascade(bounds, cnstate_vars, cn_s1, cn_s2, cn_s3, cn_s4, &
       cp_s1, cp_s2, cp_s3, cp_s4, cwd_fcel, cwd_flig, spinup_vector)
    type(bounds_type), intent(in) :: bounds
    type(cnstate_type), intent(inout) :: cnstate_vars
    real(r8), intent(in) :: cn_s1, cn_s2, cn_s3, cn_s4
    real(r8), intent(in) :: cp_s1, cp_s2, cp_s3, cp_s4
    real(r8), intent(in) :: cwd_fcel, cwd_flig
    real(r8), intent(in) :: spinup_vector(:)

    integer :: c
    integer :: pft
    real(r8) :: fractions(microbe_decomp_ntransitions)
    real(r8) :: respirations(microbe_decomp_ntransitions)

    call initialize_pool_metadata(cn_s1, cn_s2, cn_s3, cn_s4, cp_s1, cp_s2, cp_s3, cp_s4, &
         spinup_vector)
    call initialize_transition_metadata()

    do c = bounds%begc, bounds%endc
       pft = dominant_column_pft(c)
       call build_column_pathways(pft, cwd_fcel, cwd_flig, fractions, respirations)
       call validateMicrobeDecompDonorClosure(fractions)
       cnstate_vars%pathfrac_decomp_cascade_col(c,1:nlevdecomp,:) = &
            spread(fractions, dim=1, ncopies=nlevdecomp)
       cnstate_vars%rf_decomp_cascade_col(c,1:nlevdecomp,:) = &
            spread(respirations, dim=1, ncopies=nlevdecomp)
    end do

  end subroutine initMicrobeDecompCascade

  !-----------------------------------------------------------------------
  subroutine initialize_pool_metadata(cn_s1, cn_s2, cn_s3, cn_s4, cp_s1, cp_s2, cp_s3, cp_s4, &
       spinup_vector)
    real(r8), intent(in) :: cn_s1, cn_s2, cn_s3, cn_s4
    real(r8), intent(in) :: cp_s1, cp_s2, cp_s3, cp_s4
    real(r8), intent(in) :: spinup_vector(:)

    integer, parameter :: i_soil1 = 5
    integer, parameter :: i_soil2 = 6
    integer, parameter :: i_soil3 = 7
    integer, parameter :: i_soil4 = 8
    integer, parameter :: i_atm = 0

    call set_pool(i_atm, 'atmosphere', 'atmosphere', 'atmosphere', '', &
         .false., .false., .false., .false., 0._r8, 0._r8, 0._r8)
    call set_pool(i_met_lit, 'litr1', 'LITR1', 'litter 1', 'L1', &
         .true., .true., .false., .false., 90._r8, 900._r8, 0._r8)
    call set_pool(i_cel_lit, 'litr2', 'LITR2', 'litter 2', 'L2', &
         .true., .true., .false., .false., 90._r8, 900._r8, 0._r8)
    call set_pool(i_lig_lit, 'litr3', 'LITR3', 'litter 3', 'L3', &
         .true., .true., .false., .false., 90._r8, 900._r8, 0._r8)
    call set_pool(i_cwd, 'cwd', 'CWD', 'coarse woody debris', 'CWD', &
         .true., .true., .true., .false., 500._r8, 5000._r8, 0._r8)
    call set_pool(i_soil1, 'soil1', 'SOIL1', 'soil 1', 'S1', &
         .false., .true., .false., .true., cn_s1, cp_s1, 0._r8)
    call set_pool(i_soil2, 'soil2', 'SOIL2', 'soil 2', 'S2', &
         .false., .true., .false., .true., cn_s2, cp_s2, 0._r8)
    call set_pool(i_soil3, 'soil3', 'SOIL3', 'soil 3', 'S3', &
         .false., .true., .false., .true., cn_s3, cp_s3, 0._r8)
    call set_pool(i_soil4, 'soil4', 'SOIL4', 'soil 4', 'S4', &
         .false., .true., .false., .true., cn_s4, cp_s4, 10._r8)
    call set_pool(i_bacteria, 'bacteria', 'BACTERIA', 'bacteria', 'B', &
         .false., .false., .false., .true., MicrobeDecompParamsInst%bacteria_pool_cn, &
         MicrobeDecompParamsInst%cp_bacteria, MicrobeDecompParamsInst%bacteria_initial_c)
    call set_pool(i_fungi, 'fungi', 'FUNGI', 'fungi', 'F', &
         .false., .false., .false., .true., MicrobeDecompParamsInst%fungi_pool_cn, &
         MicrobeDecompParamsInst%cp_fungi, MicrobeDecompParamsInst%fungi_initial_c)
    call set_pool(i_dom, 'dom', 'DOM', 'dissolved organic C', 'DOM', &
         .false., .false., .false., .true., MicrobeDecompParamsInst%cn_dom, &
         MicrobeDecompParamsInst%cp_dom, MicrobeDecompParamsInst%dom_initial_c)

    decomp_cascade_con%is_metabolic(i_met_lit) = .true.
    decomp_cascade_con%is_cellulose(i_cel_lit) = .true.
    decomp_cascade_con%is_lignin(i_lig_lit) = .true.
    decomp_cascade_con%is_litter(i_atm) = .true.
    decomp_cascade_con%is_dissolved(i_dom) = .true.
    decomp_cascade_con%is_microbial_biomass(i_bacteria) = .true.
    decomp_cascade_con%is_microbial_biomass(i_fungi) = .true.

    decomp_cascade_con%spinup_factor(i_atm) = 1._r8
    decomp_cascade_con%spinup_factor(i_met_lit) = spinup_vector(1)
    decomp_cascade_con%spinup_factor(i_cel_lit) = spinup_vector(2)
    decomp_cascade_con%spinup_factor(i_lig_lit) = spinup_vector(3)
    decomp_cascade_con%spinup_factor(i_cwd) = spinup_vector(8)
    decomp_cascade_con%spinup_factor(i_soil1) = spinup_vector(4)
    decomp_cascade_con%spinup_factor(i_soil2) = spinup_vector(5)
    decomp_cascade_con%spinup_factor(i_soil3) = spinup_vector(6)
    decomp_cascade_con%spinup_factor(i_soil4) = spinup_vector(7)
    decomp_cascade_con%spinup_factor(i_bacteria) = 1._r8
    decomp_cascade_con%spinup_factor(i_fungi) = 1._r8
    decomp_cascade_con%spinup_factor(i_dom) = 1._r8

  contains

    subroutine set_pool(index, restart_name, history_name, long_name, short_name, &
         floating_cn, floating_cp, cwd_pool, soil_pool, cn_ratio, cp_ratio, stock)
      integer, intent(in) :: index
      character(len=*), intent(in) :: restart_name, history_name, long_name, short_name
      logical, intent(in) :: floating_cn, floating_cp, cwd_pool, soil_pool
      real(r8), intent(in) :: cn_ratio, cp_ratio, stock

      decomp_cascade_con%floating_cn_ratio_decomp_pools(index) = floating_cn
      decomp_cascade_con%floating_cp_ratio_decomp_pools(index) = floating_cp
      decomp_cascade_con%decomp_pool_name_restart(index) = restart_name
      decomp_cascade_con%decomp_pool_name_history(index) = history_name
      decomp_cascade_con%decomp_pool_name_long(index) = long_name
      decomp_cascade_con%decomp_pool_name_short(index) = short_name
      decomp_cascade_con%is_litter(index) = index >= i_met_lit .and. index <= i_lig_lit
      decomp_cascade_con%is_soil(index) = soil_pool
      decomp_cascade_con%is_cwd(index) = cwd_pool
      decomp_cascade_con%initial_cn_ratio(index) = cn_ratio
      decomp_cascade_con%initial_cp_ratio(index) = cp_ratio
      decomp_cascade_con%initial_stock(index) = stock
      decomp_cascade_con%is_metabolic(index) = .false.
      decomp_cascade_con%is_cellulose(index) = .false.
      decomp_cascade_con%is_lignin(index) = .false.
      decomp_cascade_con%is_dissolved(index) = .false.
      decomp_cascade_con%is_microbial_biomass(index) = .false.
    end subroutine set_pool

  end subroutine initialize_pool_metadata

  !-----------------------------------------------------------------------
  subroutine initialize_transition_metadata()
    integer, parameter :: i_soil1 = 5
    integer, parameter :: i_soil2 = 6
    integer, parameter :: i_soil3 = 7
    integer, parameter :: i_soil4 = 8
    integer, parameter :: i_atm = 0

    call set_transition( 1, 'CWDL2', i_cwd,      i_cel_lit)
    call set_transition( 2, 'CWDL3', i_cwd,      i_lig_lit)
    call set_transition( 3, 'L1B',   i_met_lit,  i_bacteria)
    call set_transition( 4, 'L1F',   i_met_lit,  i_fungi)
    call set_transition( 5, 'L1S1',  i_met_lit,  i_soil1)
    call set_transition( 6, 'L2B',   i_cel_lit,  i_bacteria)
    call set_transition( 7, 'L2F',   i_cel_lit,  i_fungi)
    call set_transition( 8, 'L2S2',  i_cel_lit,  i_soil2)
    call set_transition( 9, 'L3B',   i_lig_lit,  i_bacteria)
    call set_transition(10, 'L3F',   i_lig_lit,  i_fungi)
    call set_transition(11, 'L3S3',  i_lig_lit,  i_soil3)
    call set_transition(12, 'S1B',   i_soil1,    i_bacteria)
    call set_transition(13, 'S1F',   i_soil1,    i_fungi)
    call set_transition(14, 'S1S2',  i_soil1,    i_soil2)
    call set_transition(15, 'S2B',   i_soil2,    i_bacteria)
    call set_transition(16, 'S2F',   i_soil2,    i_fungi)
    call set_transition(17, 'S2S3',  i_soil2,    i_soil3)
    call set_transition(18, 'S3B',   i_soil3,    i_bacteria)
    call set_transition(19, 'S3F',   i_soil3,    i_fungi)
    call set_transition(20, 'S3S4',  i_soil3,    i_soil4)
    call set_transition(21, 'S4B',   i_soil4,    i_bacteria)
    call set_transition(22, 'S4F',   i_soil4,    i_fungi)
    call set_transition(23, 'BS1',   i_bacteria, i_soil1)
    call set_transition(24, 'FS1',   i_fungi,    i_soil1)
    call set_transition(25, 'BS2',   i_bacteria, i_soil2)
    call set_transition(26, 'FS2',   i_fungi,    i_soil2)
    call set_transition(27, 'BS3',   i_bacteria, i_soil3)
    call set_transition(28, 'FS3',   i_fungi,    i_soil3)
    call set_transition(29, 'BS4',   i_bacteria, i_soil4)
    call set_transition(30, 'FS4',   i_fungi,    i_soil4)
    call set_transition(31, 'DOMS1', i_dom,      i_soil1)
    call set_transition(32, 'DOMS2', i_dom,      i_soil2)
    call set_transition(33, 'DOMS3', i_dom,      i_soil3)
    call set_transition(34, 'DOMS4', i_dom,      i_soil4)
    call set_transition(35, 'BDOM',  i_bacteria, i_dom)
    call set_transition(36, 'FDOM',  i_fungi,    i_dom)
    call set_transition(37, 'DOMB',  i_dom,      i_bacteria)
    call set_transition(38, 'DOMF',  i_dom,      i_fungi)
    call set_transition(39, 'BATM',  i_bacteria, i_atm)
    call set_transition(40, 'FATM',  i_fungi,    i_atm)
    call set_transition(41, 'L1DOM', i_met_lit,  i_dom)
    call set_transition(42, 'L2DOM', i_cel_lit,  i_dom)
    call set_transition(43, 'L3DOM', i_lig_lit,  i_dom)
    call set_transition(44, 'S1DOM', i_soil1,    i_dom)
    call set_transition(45, 'S2DOM', i_soil2,    i_dom)
    call set_transition(46, 'S3DOM', i_soil3,    i_dom)
    call set_transition(47, 'S4DOM', i_soil4,    i_dom)

  contains

    subroutine set_transition(index, name, donor, receiver)
      integer, intent(in) :: index, donor, receiver
      character(len=*), intent(in) :: name

      decomp_cascade_con%cascade_step_name(index) = name
      decomp_cascade_con%cascade_donor_pool(index) = donor
      decomp_cascade_con%cascade_receiver_pool(index) = receiver
    end subroutine set_transition

  end subroutine initialize_transition_metadata

  !-----------------------------------------------------------------------
  subroutine build_column_pathways(pft, cwd_fcel, cwd_flig, fractions, respirations)
    integer, intent(in) :: pft
    real(r8), intent(in) :: cwd_fcel, cwd_flig
    real(r8), intent(out) :: fractions(microbe_decomp_ntransitions)
    real(r8), intent(out) :: respirations(microbe_decomp_ntransitions)

    real(r8) :: bacteria_fraction
    real(r8) :: fungi_fraction
    real(r8) :: cue
    real(r8) :: residual

    fractions = 0._r8
    respirations = 0._r8

    fractions(1) = cwd_fcel
    fractions(2) = cwd_flig

    call substrate_paths(90._r8, MicrobeDecompParamsInst%l1dom_f, &
         MicrobeDecompParamsInst%l1s1_f, pft, cue_from_cn(90._r8), &
         fractions(3), fractions(4), fractions(5), fractions(41), &
         respirations(3), respirations(4))
    call substrate_paths(90._r8, MicrobeDecompParamsInst%l2dom_f, &
         MicrobeDecompParamsInst%l2s2_f, pft, cue_from_cn(90._r8), &
         fractions(6), fractions(7), fractions(8), fractions(42), &
         respirations(6), respirations(7))
    call substrate_paths(90._r8, MicrobeDecompParamsInst%l3dom_f, &
         MicrobeDecompParamsInst%l3s3_f, pft, cue_from_cn(90._r8), &
         fractions(9), fractions(10), fractions(11), fractions(43), &
         respirations(9), respirations(10))

    call substrate_paths(decomp_cascade_con%initial_cn_ratio(5), MicrobeDecompParamsInst%s1dom_f, &
         MicrobeDecompParamsInst%s1s2_f, pft, MicrobeDecompParamsInst%m_rf_s1m(pft), &
         fractions(12), fractions(13), fractions(14), fractions(44), &
         respirations(12), respirations(13))
    call substrate_paths(decomp_cascade_con%initial_cn_ratio(6), MicrobeDecompParamsInst%s2dom_f, &
         MicrobeDecompParamsInst%s2s3_f, pft, MicrobeDecompParamsInst%m_rf_s2m(pft), &
         fractions(15), fractions(16), fractions(17), fractions(45), &
         respirations(15), respirations(16))
    call substrate_paths(decomp_cascade_con%initial_cn_ratio(7), MicrobeDecompParamsInst%s3dom_f, &
         MicrobeDecompParamsInst%s3s4_f, pft, MicrobeDecompParamsInst%m_rf_s3m(pft), &
         fractions(18), fractions(19), fractions(20), fractions(46), &
         respirations(18), respirations(19))

    call microbial_allocation(decomp_cascade_con%initial_cn_ratio(8), pft, &
         bacteria_fraction, fungi_fraction)
    residual = 1._r8 - MicrobeDecompParamsInst%s4dom_f
    fractions(21) = bacteria_fraction * residual
    fractions(22) = fungi_fraction * residual
    fractions(47) = MicrobeDecompParamsInst%s4dom_f
    cue = MicrobeDecompParamsInst%m_rf_s4m(pft)
    respirations(21:22) = 1._r8 - cue

    fractions(23) = MicrobeDecompParamsInst%m_bs1_f(pft)
    fractions(25) = MicrobeDecompParamsInst%m_bs2_f(pft)
    fractions(27) = MicrobeDecompParamsInst%m_bs3_f(pft)
    fractions(29) = max(0._r8, 1._r8 - MicrobeDecompParamsInst%m_batm_f(pft) - &
         MicrobeDecompParamsInst%m_bdom_f(pft) - MicrobeDecompParamsInst%m_bs1_f(pft) - &
         MicrobeDecompParamsInst%m_bs2_f(pft) - MicrobeDecompParamsInst%m_bs3_f(pft))
    fractions(35) = MicrobeDecompParamsInst%m_bdom_f(pft)
    fractions(39) = MicrobeDecompParamsInst%m_batm_f(pft)
    respirations(39) = 1._r8

    fractions(24) = MicrobeDecompParamsInst%m_fs1_f(pft)
    fractions(26) = MicrobeDecompParamsInst%m_fs2_f(pft)
    fractions(28) = MicrobeDecompParamsInst%m_fs3_f(pft)
    fractions(30) = max(0._r8, 1._r8 - MicrobeDecompParamsInst%m_fatm_f(pft) - &
         MicrobeDecompParamsInst%m_fdom_f(pft) - MicrobeDecompParamsInst%m_fs1_f(pft) - &
         MicrobeDecompParamsInst%m_fs2_f(pft) - MicrobeDecompParamsInst%m_fs3_f(pft))
    fractions(36) = MicrobeDecompParamsInst%m_fdom_f(pft)
    fractions(40) = MicrobeDecompParamsInst%m_fatm_f(pft)
    respirations(40) = 1._r8

    fractions(31) = MicrobeDecompParamsInst%m_doms1_f(pft)
    fractions(32) = MicrobeDecompParamsInst%m_doms2_f(pft)
    fractions(33) = MicrobeDecompParamsInst%m_doms3_f(pft)
    fractions(34) = max(0._r8, 1._r8 - MicrobeDecompParamsInst%m_domb_f(pft) - &
         MicrobeDecompParamsInst%m_domf_f(pft) - MicrobeDecompParamsInst%m_doms1_f(pft) - &
         MicrobeDecompParamsInst%m_doms2_f(pft) - MicrobeDecompParamsInst%m_doms3_f(pft))
    fractions(37) = MicrobeDecompParamsInst%m_domb_f(pft)
    fractions(38) = MicrobeDecompParamsInst%m_domf_f(pft)

  contains

    subroutine substrate_paths(substrate_cn, dom_fraction, direct_fraction, pft_index, cue_value, &
         bacteria_path, fungi_path, direct_path, dom_path, bacteria_respiration, fungi_respiration)
      real(r8), intent(in) :: substrate_cn, dom_fraction, direct_fraction, cue_value
      integer, intent(in) :: pft_index
      real(r8), intent(out) :: bacteria_path, fungi_path, direct_path, dom_path
      real(r8), intent(out) :: bacteria_respiration, fungi_respiration

      real(r8) :: available
      real(r8) :: bacteria_share
      real(r8) :: fungi_share

      call microbial_allocation(substrate_cn, pft_index, bacteria_share, fungi_share)
      available = 1._r8 - dom_fraction - direct_fraction
      bacteria_path = available * bacteria_share
      fungi_path = available * fungi_share
      direct_path = direct_fraction
      dom_path = dom_fraction
      bacteria_respiration = 1._r8 - cue_value
      fungi_respiration = 1._r8 - cue_value
    end subroutine substrate_paths

  end subroutine build_column_pathways

  !-----------------------------------------------------------------------
  subroutine microbial_allocation(substrate_cn, pft, bacteria_fraction, fungi_fraction)
    real(r8), intent(in) :: substrate_cn
    integer, intent(in) :: pft
    real(r8), intent(out) :: bacteria_fraction
    real(r8), intent(out) :: fungi_fraction

    real(r8) :: bacteria_weight
    real(r8) :: fungi_weight

    bacteria_weight = (MicrobeDecompParamsInst%cn_bacteria(pft) / substrate_cn) ** &
         MicrobeDecompParamsInst%allocation_cn_exponent
    fungi_weight = (MicrobeDecompParamsInst%cn_fungi(pft) / substrate_cn) ** &
         MicrobeDecompParamsInst%allocation_cn_exponent
    bacteria_fraction = bacteria_weight / (bacteria_weight + fungi_weight)
    fungi_fraction = 1._r8 - bacteria_fraction
  end subroutine microbial_allocation

  !-----------------------------------------------------------------------
  real(r8) function cue_from_cn(substrate_cn)
    real(r8), intent(in) :: substrate_cn

    cue_from_cn = MicrobeDecompParamsInst%cue_max * min(1._r8, &
         MicrobeDecompParamsInst%cue_cn_target / &
         (substrate_cn * MicrobeDecompParamsInst%cue_max))
  end function cue_from_cn

  !-----------------------------------------------------------------------
  subroutine validateMicrobeDecompDonorClosure(fractions)
    real(r8), intent(in) :: fractions(microbe_decomp_ntransitions)

    integer :: donor
    integer :: transition
    real(r8) :: total
    real(r8), parameter :: tolerance = 100._r8 * epsilon(1._r8)

    if (any(fractions < -tolerance) .or. any(fractions > 1._r8 + tolerance)) then
       call endrun(msg=' ERROR: microbial decomposition path fraction is outside [0,1].'//&
            errMsg(__FILE__, __LINE__))
    end if

    do donor = 1, microbe_decomp_npools
       total = 0._r8
       do transition = 1, microbe_decomp_ntransitions
          if (decomp_cascade_con%cascade_donor_pool(transition) == donor) then
             total = total + fractions(transition)
          end if
       end do
       if (abs(total - 1._r8) > tolerance) then
          call endrun(msg=' ERROR: microbial decomposition donor path fractions do not sum to one.'//&
               errMsg(__FILE__, __LINE__))
       end if
    end do
  end subroutine validateMicrobeDecompDonorClosure

  !-----------------------------------------------------------------------
  subroutine getMicrobeDecompTimestepRates(c, dtd, dom_rate, bacteria_rate, fungi_rate)
    integer, intent(in) :: c
    real(r8), intent(in) :: dtd
    real(r8), intent(out) :: dom_rate, bacteria_rate, fungi_rate

    integer :: pft

    pft = dominant_column_pft(c)
    dom_rate = 1._r8 - (1._r8 - MicrobeDecompParamsInst%k_dom(pft)) ** dtd
    bacteria_rate = 1._r8 - (1._r8 - MicrobeDecompParamsInst%k_bacteria(pft)) ** dtd
    fungi_rate = 1._r8 - (1._r8 - MicrobeDecompParamsInst%k_fungi(pft)) ** dtd
  end subroutine getMicrobeDecompTimestepRates

  !-----------------------------------------------------------------------
  integer function dominant_column_pft(c)
    integer, intent(in) :: c

    integer :: p
    real(r8) :: largest_weight

    dominant_column_pft = 0
    if (.not. (col_pp%is_soil(c) .or. col_pp%is_crop(c))) return

    largest_weight = -1._r8
    do p = col_pp%pfti(c), col_pp%pftf(c)
       if (veg_pp%wtcol(p) > largest_weight) then
          dominant_column_pft = veg_pp%itype(p)
          largest_weight = veg_pp%wtcol(p)
       end if
    end do

    if (dominant_column_pft < lbound(MicrobeDecompParamsInst%k_dom, 1) .or. &
         dominant_column_pft > ubound(MicrobeDecompParamsInst%k_dom, 1)) then
       call endrun(msg=' ERROR: column PFT index is outside the microbial parameter dimension.'//&
            errMsg(__FILE__, __LINE__))
    end if
  end function dominant_column_pft

  !-----------------------------------------------------------------------
  subroutine require_fraction(name, value)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: value

    if (.not. ieee_is_finite(value) .or. value < 0._r8 .or. value > 1._r8) then
       call endrun(msg=' ERROR: microbial parameter must be in [0,1]: '//trim(name)//&
            errMsg(__FILE__, __LINE__))
    end if
  end subroutine require_fraction

  subroutine require_fraction_array(name, values)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: values(:)

    if (any(.not. ieee_is_finite(values)) .or. &
         any(values < 0._r8) .or. any(values > 1._r8)) then
       call endrun(msg=' ERROR: microbial PFT parameter must be in [0,1]: '//trim(name)//&
            errMsg(__FILE__, __LINE__))
    end if
  end subroutine require_fraction_array

  subroutine require_probability_array(name, values, allow_one)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: values(:)
    logical, intent(in) :: allow_one

    if (any(.not. ieee_is_finite(values)) .or. &
         any(values < 0._r8) .or. any(values > 1._r8) .or. &
         (.not. allow_one .and. any(values == 1._r8))) then
       call endrun(msg=' ERROR: microbial daily turnover parameter must be in [0,1): '//trim(name)//&
            errMsg(__FILE__, __LINE__))
    end if
  end subroutine require_probability_array

  subroutine require_positive(name, value)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: value

    if (.not. ieee_is_finite(value) .or. value <= 0._r8) then
       call endrun(msg=' ERROR: microbial parameter must be positive: '//trim(name)//&
            errMsg(__FILE__, __LINE__))
    end if
  end subroutine require_positive

  subroutine require_positive_array(name, values)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: values(:)

    if (any(.not. ieee_is_finite(values)) .or. any(values <= 0._r8)) then
       call endrun(msg=' ERROR: microbial PFT parameter must be positive: '//trim(name)//&
            errMsg(__FILE__, __LINE__))
    end if
  end subroutine require_positive_array

  subroutine require_nonnegative(name, value)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: value

    if (.not. ieee_is_finite(value) .or. value < 0._r8) then
       call endrun(msg=' ERROR: microbial parameter must be nonnegative: '//trim(name)//&
            errMsg(__FILE__, __LINE__))
    end if
  end subroutine require_nonnegative

  subroutine require_nonnegative_residual(name, residual, pft)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: residual
    integer, optional, intent(in) :: pft

    character(len=32) :: pft_context
    real(r8), parameter :: tolerance = 100._r8 * epsilon(1._r8)

    if (.not. ieee_is_finite(residual) .or. residual < -tolerance) then
       if (present(pft)) then
          write(pft_context, '(a,i0)') ' at PFT ', pft
       else
          pft_context = ''
       end if
       call endrun(msg=' ERROR: microbial pathway fractions have a negative residual for '//&
            trim(name)//trim(pft_context)//errMsg(__FILE__, __LINE__))
    end if
  end subroutine require_nonnegative_residual

end module MicrobeDecompMod
