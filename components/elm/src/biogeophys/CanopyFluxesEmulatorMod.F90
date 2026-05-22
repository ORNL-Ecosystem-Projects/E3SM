module CanopyFluxesEmulatorMod

  ! Placeholder hook for a future canopy-flux emulator.
  !
  ! The V1 emulator interface is intentionally limited to downstream-required
  ! non-FATES outputs:
  ! - canopy water / irrigation state needed by hydrology
  ! - canopy latent and sensible fluxes needed by soil-energy coupling
  ! - soil-moisture stress and root weighting diagnostics
  ! - aerodynamic exchange terms consumed by chemistry / CH4 / map-to-atm
  ! - aggregate photosynthesis and stomatal outputs needed by BGC
  !
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_const_mod         , only : SHR_CONST_RGAS
  use decompMod             , only : bounds_type
  use atm2lndType           , only : atm2lnd_type
  use CanopyStateType       , only : canopystate_type, perchroot, perchroot_alt
  use CNStateType           , only : cnstate_type
  use EnergyFluxType        , only : energyflux_type
  use FrictionvelocityType  , only : frictionvel_type
  use SoilStateType         , only : soilstate_type
  use SoilHydrologyType     , only : soilhydrology_type
  use SolarAbsorbedType     , only : solarabs_type
  use SurfaceAlbedoType     , only : surfalb_type
  use CH4Mod                , only : ch4_type
  use PhotosynthesisType    , only : photosyns_type
  use VegetationType        , only : veg_pp
  use LandunitType          , only : lun_pp
  use landunit_varcon       , only : istsoil
  use WaterfluxType         , only : waterflux_vars
  use VegetationPropertiesType, only : veg_vp
  use VegetationDataType    , only : veg_es, veg_ws, veg_wf, veg_ef
  use ColumnType            , only : col_pp
  use ColumnDataType        , only : col_ws, col_es, col_ef, col_wf
  use CanopyFluxesMod       , only : CanopyFluxes
  use PhotosynthesisMod     , only : params_inst
  use SoilMoistStressMod    , only : calc_effective_soilporosity, calc_volumetric_h2oliq
  use SoilMoistStressMod    , only : calc_root_moist_stress, set_perchroot_opt
  use elm_varpar            , only : numpft, nlevgrnd
   use elm_varcon            , only : denice, denh2o, sb, cpair, hvap, secspday
  use spmdMod               , only : masterproc, mpicom
   use elm_time_manager      , only : get_curr_date, get_curr_calday, get_days_per_year
   use elm_time_manager      , only : get_nstep, get_step_size
   use elm_varctl            , only : caseid
  use shr_mpi_mod           , only : shr_mpi_gathScatVInit, shr_mpi_gatherv
  use netcdf                , only : nf90_clobber, nf90_noerr, nf90_double, nf90_int, nf90_global
  use netcdf                , only : nf90_create, nf90_open, nf90_def_dim, nf90_def_var
  use netcdf                , only : nf90_put_var, nf90_put_att, nf90_enddef, nf90_inq_dimid
  use netcdf                , only : nf90_inquire_dimension, nf90_inq_varid, nf90_close
  use netcdf                , only : nf90_strerror, nf90_write, nf90_unlimited, nf90_nowrite
  use netcdf                , only : nf90_get_var, nf90_get_att, nf90_inquire_variable

  implicit none
  private

   logical, public :: use_canopyflux_emulator = .false.
  logical, public :: write_canopyflux_training_data = .true.
  logical, public :: randomize_canopyflux_training_traits = .true.
  logical, public :: debug_canopyflux_emulator = .false.
  integer, public :: debug_canopyflux_patch = 2
  integer, public :: debug_canopyflux_compare_nstep = 8
  character(len=256), public :: canopyflux_model_dir = '/code/E3SM-Peatlands/components/elm/tools/canopyflux_emulator/canopyflux_models/'

  integer, parameter, public :: canopyflux_emulator_profile_layers = 10
  integer, parameter, public :: canopyflux_emulator_rootfr_layers = 10
  integer, parameter, public :: canopyflux_emulator_base_features = 28
  integer, parameter, public :: canopyflux_emulator_post_profile_features = 6
  integer, parameter, public :: canopyflux_emulator_num_features = canopyflux_emulator_base_features + &
       4 * canopyflux_emulator_profile_layers + canopyflux_emulator_rootfr_layers + &
       canopyflux_emulator_post_profile_features
  integer, parameter, public :: canopyflux_emulator_num_targets = 18
  integer, parameter, public :: canopyflux_emulator_num_debug = 23
  integer, parameter, public :: canopyflux_emulator_num_params = 31
  real(r8), allocatable, save :: captured_features(:,:)
  real(r8), allocatable, save :: captured_targets(:,:)
  real(r8), allocatable, save :: captured_debug(:,:)
  logical , allocatable, save :: captured_natveg_mask(:)
  integer, save :: captured_begp = 0
  integer, save :: captured_endp = -1
  logical, save :: capture_ready = .false.
  logical, save :: reported_required_pfts = .false.
  logical, save :: reported_emulator_outputs = .false.
  logical, save :: canopyflux_training_file_initialized = .false.
  integer, save :: info3330_reset_nstep = -1
  real(r8), save :: canopyflux_native_time_total = 0._r8
  real(r8), save :: canopyflux_emulator_time_total = 0._r8
  real(r8), save :: canopyflux_emulator_inference_time_total = 0._r8
  real(r8), save :: canopyflux_emulator_apply_time_total = 0._r8
  real(r8), save :: canopyflux_emulator_rebuild_time_total = 0._r8
  integer , save :: canopyflux_native_call_count = 0
  integer , save :: canopyflux_emulator_call_count = 0
  integer , save :: canopyflux_native_patch_count_total = 0
  integer , save :: canopyflux_emulator_patch_count_total = 0
  integer , save :: canopyflux_emulator_group_count_total = 0
  integer , save :: canopyflux_emulator_group_patch_count_total = 0
  integer , save :: canopyflux_emulator_max_group_size = 0

  type :: canopyflux_dense_layer_type
     real(r8), allocatable :: weight(:,:)
     real(r8), allocatable :: bias(:)
  end type canopyflux_dense_layer_type

  type :: canopyflux_model_type
     logical :: is_loaded = .false.
     integer :: patch_index = -1
     integer, allocatable :: input_map(:)
     integer, allocatable :: target_map(:)
     integer, allocatable :: nonnegative_target(:)
     real(r8), allocatable :: x_mean(:)
     real(r8), allocatable :: x_std(:)
     real(r8), allocatable :: y_mean(:)
     real(r8), allocatable :: y_std(:)
     type(canopyflux_dense_layer_type), allocatable :: layers(:)
  end type canopyflux_model_type

  type(canopyflux_model_type), allocatable, save :: canopyflux_models(:)
  integer, parameter :: canopyflux_timing_repeat_count = 1000
  real(r8), public :: canopyflux_day_threshold = 1._r8
  real(r8), parameter :: canopyflux_training_flnr_min = 0.03_r8
  real(r8), parameter :: canopyflux_training_flnr_max = 0.2_r8
  real(r8), parameter :: canopyflux_training_mbbopt_min = 4._r8
  real(r8), parameter :: canopyflux_training_mbbopt_max = 13._r8
   real(r8), parameter :: canopyflux_training_pco2_min = 18._r8
   real(r8), parameter :: canopyflux_training_pco2_max = 80._r8
   real(r8), public :: canopyflux_training_start_delay_years = 10._r8
   logical, save :: canopyflux_training_start_initialized = .false.
   integer, save :: canopyflux_training_start_year = -huge(1)
   real(r8), save :: canopyflux_training_start_calday = 0._r8
   integer, save :: canopyflux_training_start_sec = 0

  character(len=*), parameter :: canopyflux_feature_names = &
       'forc_t,forc_q,forc_pbot,forc_lwrad,forc_rain,forc_snow,forc_swrad,' // &
       'forc_solad_vis,forc_solad_nir,forc_solai_vis,forc_solai_nir,' // &
       'wind_speed,forc_rho,forc_pco2,frac_veg_nosno_patch,elai_patch,esai_patch,frac_sno,snow_depth,' // &
       'zwt_col,h2ocan,t10,dayl,max_dayl,forc_hgt_u_patch,z0mg_col,frac_h2osfc,coszen_col,' // &
       'smp_l_lev1,smp_l_lev2,smp_l_lev3,smp_l_lev4,smp_l_lev5,smp_l_lev6,smp_l_lev7,smp_l_lev8,' // &
       'smp_l_lev9,smp_l_lev10,' // &
       'eff_porosity_lev1,eff_porosity_lev2,eff_porosity_lev3,eff_porosity_lev4,' // &
       'eff_porosity_lev5,eff_porosity_lev6,eff_porosity_lev7,eff_porosity_lev8,' // &
       'eff_porosity_lev9,eff_porosity_lev10,' // &
       'h2osoi_liqvol_lev1,h2osoi_liqvol_lev2,h2osoi_liqvol_lev3,h2osoi_liqvol_lev4,' // &
       'h2osoi_liqvol_lev5,h2osoi_liqvol_lev6,h2osoi_liqvol_lev7,h2osoi_liqvol_lev8,' // &
       'h2osoi_liqvol_lev9,h2osoi_liqvol_lev10,' // &
       'rootfr_lev1,rootfr_lev2,rootfr_lev3,rootfr_lev4,rootfr_lev5,rootfr_lev6,rootfr_lev7,rootfr_lev8,' // &
       'rootfr_lev9,rootfr_lev10,' // &
       't_grnd,t_h2osfc,ugust,' // &
       't_soisno_lev1,t_soisno_lev2,t_soisno_lev3,t_soisno_lev4,t_soisno_lev5,t_soisno_lev6,t_soisno_lev7,t_soisno_lev8,' // &
       't_soisno_lev9,t_soisno_lev10,' // &
       'flnr,mbbopt,sabg_snow'

  character(len=*), parameter :: canopyflux_target_names = &
       'fpsn_patch,eflx_sh_tot,qflx_tran_veg,qflx_evap_canopy,' // &
       'qflx_ev_soil,qflx_ev_snow,qflx_ev_h2osfc,qflx_prec_grnd,' // &
       't_veg,t_ref2m,hs_soil,hs_top_snow,hs_h2osfc,dhsdT,' // &
       'laisun_patch,laisha_patch,psnsun_patch,psnsha_patch'
  real(r8), parameter :: canopyflux_emulator_elai_native_fallback = 1.0e-8_r8
  real(r8), parameter :: canopyflux_emulator_elai_fpsn_cutoff = 1.0e-2_r8
  real(r8), parameter :: canopyflux_emulator_frac_veg_nosno_native_fallback = tiny(1._r8)

  character(len=*), parameter :: canopyflux_param_names = &
       'pft_type,c3psn,nonvascular,dleaf,smpso,smpsc,slatop,leafcn,flnr,fnitr,' // &
       'br_mr_pft,q10_mr_pft,i_vc,s_vc,act25,kcha,koha,cpha,vcmaxha,jmaxha,' // &
       'tpuha,lmrha,vcmaxhd,jmaxhd,tpuhd,lmrhd,lmrse,qe,theta_cj,bbbopt,' // &
       'mbbopt'

  type, public :: canopyflux_emulator_input_view_type
     real(r8), pointer :: forc_lwrad(:) => null()
     real(r8), pointer :: forc_rain(:) => null()
     real(r8), pointer :: forc_snow(:) => null()
     real(r8), pointer :: forc_swrad(:) => null()
     real(r8), pointer :: forc_solad(:,:) => null()
     real(r8), pointer :: forc_solai(:,:) => null()
     real(r8), pointer :: forc_q(:) => null()
     real(r8), pointer :: forc_pbot(:) => null()
     real(r8), pointer :: forc_rho(:) => null()
     real(r8), pointer :: forc_t(:) => null()
     real(r8), pointer :: forc_u(:) => null()
     real(r8), pointer :: forc_v(:) => null()
     real(r8), pointer :: wsresp(:) => null()
     real(r8), pointer :: tau_est(:) => null()
     real(r8), pointer :: ugust(:) => null()
     real(r8), pointer :: forc_pco2(:) => null()
     real(r8), pointer :: sabv_patch(:) => null()
     real(r8), pointer :: sabg_snow_patch(:) => null()
     integer , pointer :: frac_veg_nosno_patch(:) => null()
     real(r8), pointer :: elai_patch(:) => null()
     real(r8), pointer :: esai_patch(:) => null()
     real(r8), pointer :: laisun_patch(:) => null()
     real(r8), pointer :: laisha_patch(:) => null()
     real(r8), pointer :: displa_patch(:) => null()
     real(r8), pointer :: htop_patch(:) => null()
     real(r8), pointer :: fwet(:) => null()
     real(r8), pointer :: fdry(:) => null()
     real(r8), pointer :: frac_sno(:) => null()
     real(r8), pointer :: snow_depth(:) => null()
     real(r8), pointer :: h2osfc(:) => null()
     real(r8), pointer :: frac_h2osfc_act(:) => null()
     real(r8), pointer :: h2osoi_ice(:,:) => null()
     real(r8), pointer :: h2osoi_liq(:,:) => null()
     real(r8), pointer :: h2osoi_vol(:,:) => null()
     real(r8), pointer :: h2osoi_liqvol(:,:) => null()
     real(r8), pointer :: zwt_col(:) => null()
     real(r8), pointer :: watsat_col(:,:) => null()
     real(r8), pointer :: smp_l_col(:,:) => null()
     real(r8), pointer :: eff_porosity_col(:,:) => null()
     real(r8), pointer :: rootfr_patch(:,:) => null()
     real(r8), pointer :: soilbeta_col(:) => null()
     real(r8), pointer :: h2ocan(:) => null()
     real(r8), pointer :: qflx_tran_veg(:) => null()
     real(r8), pointer :: btran_patch(:) => null()
     real(r8), pointer :: h2o_moss_wc(:) => null()
     real(r8), pointer :: t10(:) => null()
     real(r8), pointer :: t_veg(:) => null()
     real(r8), pointer :: dayl(:) => null()
     real(r8), pointer :: max_dayl(:) => null()
     real(r8), pointer :: forc_hgt_u_patch(:) => null()
     real(r8), pointer :: z0mg_col(:) => null()
     real(r8), pointer :: frac_h2osfc(:) => null()
     real(r8), pointer :: coszen_col(:) => null()
     real(r8), pointer :: t_h2osfc(:) => null()
     real(r8), pointer :: t_soisno(:,:) => null()
     real(r8), pointer :: t_grnd(:) => null()
     real(r8), pointer :: thm(:) => null()
     real(r8), pointer :: emv(:) => null()
     real(r8), pointer :: emg(:) => null()
     real(r8), pointer :: qg_snow(:) => null()
     real(r8), pointer :: qg_soil(:) => null()
     real(r8), pointer :: qg_h2osfc(:) => null()
     real(r8), pointer :: qg(:) => null()
     real(r8), pointer :: dqgdT(:) => null()
  contains
     procedure, public :: bind => bind_canopyflux_emulator_input_view
  end type canopyflux_emulator_input_view_type

  type, public :: canopyflux_emulator_output_view_type
     real(r8), pointer :: ram1_patch(:) => null()
     real(r8), pointer :: rb1_patch(:) => null()
     real(r8), pointer :: h2ocan(:) => null()
     integer , pointer :: n_irrig_steps_left(:) => null()
     real(r8), pointer :: irrig_rate(:) => null()
     real(r8), pointer :: qflx_tran_veg(:) => null()
     real(r8), pointer :: qflx_evap_veg(:) => null()
     real(r8), pointer :: qflx_evap_soi(:) => null()
     real(r8), pointer :: qflx_ev_snow(:) => null()
     real(r8), pointer :: qflx_ev_soil(:) => null()
     real(r8), pointer :: qflx_ev_h2osfc(:) => null()
     real(r8), pointer :: qflx_prec_grnd(:) => null()
     real(r8), pointer :: qflx_h2osfc_surf_col(:) => null()
     real(r8), pointer :: qflx_gross_infl_soil_col(:) => null()
     real(r8), pointer :: qflx_adv_col(:,:) => null()
     real(r8), pointer :: qflx_rootsoi_col(:,:) => null()
     real(r8), pointer :: grnd_ch4_cond_patch(:) => null()
     real(r8), pointer :: btran_patch(:) => null()
     real(r8), pointer :: btran2_patch(:) => null()
     real(r8), pointer :: rresis_patch(:,:) => null()
     real(r8), pointer :: canopy_cond_patch(:) => null()
     real(r8), pointer :: eff_porosity_col(:,:) => null()
     real(r8), pointer :: rootr_patch(:,:) => null()
     real(r8), pointer :: h2osoi_liqvol(:,:) => null()
     real(r8), pointer :: cgrnds(:) => null()
     real(r8), pointer :: cgrndl(:) => null()
     real(r8), pointer :: dlrad(:) => null()
     real(r8), pointer :: ulrad(:) => null()
     real(r8), pointer :: cgrnd(:) => null()
     real(r8), pointer :: eflx_sh_snow(:) => null()
     real(r8), pointer :: eflx_sh_h2osfc(:) => null()
     real(r8), pointer :: eflx_sh_soil(:) => null()
     real(r8), pointer :: eflx_sh_veg(:) => null()
     real(r8), pointer :: eflx_sh_grnd(:) => null()
     real(r8), pointer :: eflx_sh_tot(:) => null()
     real(r8), pointer :: eflx_hs_soil_col(:) => null()
     real(r8), pointer :: eflx_hs_top_snow_col(:) => null()
     real(r8), pointer :: eflx_hs_h2osfc_col(:) => null()
     real(r8), pointer :: eflx_dhsdT_col(:) => null()
     real(r8), pointer :: rssun_patch(:) => null()
     real(r8), pointer :: rssha_patch(:) => null()
     real(r8), pointer :: lmrsun_patch(:) => null()
     real(r8), pointer :: lmrsha_patch(:) => null()
     real(r8), pointer :: laisun_patch(:) => null()
     real(r8), pointer :: laisha_patch(:) => null()
     real(r8), pointer :: fpsn_patch(:) => null()
     real(r8), pointer :: psnsun_patch(:) => null()
     real(r8), pointer :: psnsha_patch(:) => null()
     real(r8), pointer :: t_veg(:) => null()
     real(r8), pointer :: t_ref2m(:) => null()
   contains
     procedure, public :: bind => bind_canopyflux_emulator_output_view
  end type canopyflux_emulator_output_view_type

  real(r8), parameter :: canopyflux_min_conductance = 1.e-6_r8

  public :: CanopyFluxesEmulator
  public :: predict_canopyfluxes
  public :: bind_canopyflux_emulator_input_view
  public :: bind_canopyflux_emulator_output_view
  public :: assemble_canopyflux_emulator_features
  public :: assemble_canopyflux_emulator_targets
  public :: assemble_canopyflux_emulator_params
  public :: begin_canopyflux_training_capture
  public :: capture_canopyflux_training_debug
  public :: finalize_canopyflux_training_capture
  public :: flush_canopyflux_training_capture

contains

  logical pure function is_natveg_patch(p) result(is_natveg)

    integer, intent(in) :: p

    is_natveg = (lun_pp%itype(veg_pp%landunit(p)) == istsoil)

  end function is_natveg_patch

  integer pure function capture_slot(p) result(slot)

    integer, intent(in) :: p

    slot = p - captured_begp + 1

  end function capture_slot

  pure real(r8) function patch_wind_speed(u_patch, v_patch) result(speed)

    real(r8), intent(in) :: u_patch
    real(r8), intent(in) :: v_patch

    speed = sqrt(u_patch * u_patch + v_patch * v_patch)

  end function patch_wind_speed

  pure real(r8) function canopyflux_conductance_from_ram1(ram1) result(gh_canopy)

    real(r8), intent(in) :: ram1

    gh_canopy = 1._r8 / max(ram1, canopyflux_min_conductance)

  end function canopyflux_conductance_from_ram1

  pure real(r8) function canopyflux_ram1_from_conductance(gh_canopy) result(ram1)

    real(r8), intent(in) :: gh_canopy

    ram1 = 1._r8 / max(gh_canopy, canopyflux_min_conductance)

  end function canopyflux_ram1_from_conductance

  logical function canopyflux_training_capture_enabled() result(enabled)

    integer :: year, mon, day, sec
    real(r8) :: elapsed_years, days_per_year

    if (.not. write_canopyflux_training_data) then
       enabled = .false.
       return
    end if

    call get_curr_date(year, mon, day, sec)
    if (.not. canopyflux_training_start_initialized) then
       canopyflux_training_start_year = year
       canopyflux_training_start_calday = get_curr_calday()
       canopyflux_training_start_sec = sec
       canopyflux_training_start_initialized = .true.
    end if

    if (canopyflux_training_start_delay_years <= 0._r8) then
       enabled = .true.
       return
    end if

    days_per_year = real(get_days_per_year(), r8)
    elapsed_years = real(year - canopyflux_training_start_year, r8) + &
         (((get_curr_calday() - canopyflux_training_start_calday) * secspday) + &
         real(sec - canopyflux_training_start_sec, r8)) / (days_per_year * secspday)
    if (elapsed_years < 0._r8) elapsed_years = 0._r8
    enabled = (elapsed_years >= canopyflux_training_start_delay_years)

  end function canopyflux_training_capture_enabled

  subroutine bind_canopyflux_emulator_input_view(this, atm2lnd_vars, canopystate_vars, energyflux_vars, &
       frictionvel_vars, soilstate_vars, solarabs_vars, surfalb_vars, soilhydrology_vars)

    use TopounitDataType      , only : top_as, top_af
    use GridcellType          , only : grc_pp

    class(canopyflux_emulator_input_view_type), intent(inout) :: this
    type(atm2lnd_type)       , intent(inout), target :: atm2lnd_vars
    type(canopystate_type)   , intent(inout), target :: canopystate_vars
    type(energyflux_type)    , intent(inout), target :: energyflux_vars
    type(frictionvel_type)   , intent(inout), target :: frictionvel_vars
    type(soilstate_type)     , intent(inout), target :: soilstate_vars
    type(solarabs_type)      , intent(inout), target :: solarabs_vars
    type(surfalb_type)       , intent(inout), target :: surfalb_vars
    type(soilhydrology_type) , intent(in)   , target :: soilhydrology_vars

    this%forc_lwrad => top_af%lwrad
    this%forc_rain => top_af%rain
    this%forc_snow => top_af%snow
    this%forc_swrad => top_af%solar
    this%forc_solad => top_af%solad
    this%forc_solai => top_af%solai
    this%forc_q => top_as%qbot
    this%forc_pbot => top_as%pbot
    this%forc_rho => top_as%rhobot
    this%forc_t => top_as%tbot
    this%forc_u => top_as%ubot
    this%forc_v => top_as%vbot
    this%wsresp => top_as%wsresp
    this%tau_est => top_as%tau_est
    this%ugust => top_as%ugust
    this%forc_pco2 => top_as%pco2bot
    this%sabv_patch => solarabs_vars%sabv_patch
    this%sabg_snow_patch => solarabs_vars%sabg_snow_patch
    this%frac_veg_nosno_patch => canopystate_vars%frac_veg_nosno_patch
    this%elai_patch => canopystate_vars%elai_patch
    this%esai_patch => canopystate_vars%esai_patch
    this%laisun_patch => canopystate_vars%laisun_patch
    this%laisha_patch => canopystate_vars%laisha_patch
    this%displa_patch => canopystate_vars%displa_patch
    this%htop_patch => canopystate_vars%htop_patch
    this%fwet => veg_ws%fwet
    this%fdry => veg_ws%fdry
    this%frac_sno => col_ws%frac_sno_eff
    this%snow_depth => col_ws%snow_depth
    this%h2osfc => col_ws%h2osfc
    this%frac_h2osfc_act => col_ws%frac_h2osfc_act
    this%h2osoi_ice => col_ws%h2osoi_ice
    this%h2osoi_liq => col_ws%h2osoi_liq
    this%h2osoi_vol => col_ws%h2osoi_vol
    this%h2osoi_liqvol => col_ws%h2osoi_liqvol
    this%zwt_col => soilhydrology_vars%zwt_col
    this%watsat_col => soilstate_vars%watsat_col
    this%smp_l_col => soilstate_vars%smp_l_col
    this%eff_porosity_col => soilstate_vars%eff_porosity_col
    this%rootfr_patch => soilstate_vars%rootfr_patch
    this%soilbeta_col => soilstate_vars%soilbeta_col
    this%h2ocan => veg_ws%h2ocan
    this%qflx_tran_veg => veg_wf%qflx_tran_veg
    this%btran_patch => energyflux_vars%btran_patch
    this%t10 => veg_es%t_a10
    this%t_veg => veg_es%t_veg
    this%dayl => grc_pp%dayl
    this%max_dayl => grc_pp%max_dayl
    this%forc_hgt_u_patch => frictionvel_vars%forc_hgt_u_patch
    this%z0mg_col => frictionvel_vars%z0mg_col
    this%frac_h2osfc => col_ws%frac_h2osfc
    this%coszen_col => surfalb_vars%coszen_col
    this%t_h2osfc => col_es%t_h2osfc
    this%t_soisno => col_es%t_soisno
    this%t_grnd => col_es%t_grnd
    this%thm => veg_es%thm
    this%emv => veg_es%emv
    this%emg => col_es%emg
    this%qg_snow => col_ws%qg_snow
    this%qg_soil => col_ws%qg_soil
    this%qg_h2osfc => col_ws%qg_h2osfc
    this%qg => col_ws%qg
    this%dqgdT => col_ws%dqgdT

  end subroutine bind_canopyflux_emulator_input_view

  subroutine bind_canopyflux_emulator_output_view(this, canopystate_vars, energyflux_vars, &
       frictionvel_vars, soilstate_vars, ch4_vars, photosyns_vars)

    class(canopyflux_emulator_output_view_type), intent(inout) :: this
    type(canopystate_type)   , intent(inout), target :: canopystate_vars
    type(energyflux_type)    , intent(inout), target :: energyflux_vars
    type(frictionvel_type)   , intent(inout), target :: frictionvel_vars
    type(soilstate_type)     , intent(inout), target :: soilstate_vars
    type(ch4_type)           , intent(inout), target :: ch4_vars
    type(photosyns_type)     , intent(inout), target :: photosyns_vars

    this%ram1_patch => frictionvel_vars%ram1_patch
    this%rb1_patch => frictionvel_vars%rb1_patch
    this%h2ocan => veg_ws%h2ocan
    this%n_irrig_steps_left => veg_wf%n_irrig_steps_left
    this%irrig_rate => veg_wf%irrig_rate
    this%qflx_tran_veg => veg_wf%qflx_tran_veg
    this%qflx_evap_veg => veg_wf%qflx_evap_veg
    this%qflx_evap_soi => veg_wf%qflx_evap_soi
    this%qflx_ev_snow => veg_wf%qflx_ev_snow
    this%qflx_ev_soil => veg_wf%qflx_ev_soil
    this%qflx_ev_h2osfc => veg_wf%qflx_ev_h2osfc
    this%qflx_prec_grnd => veg_wf%qflx_prec_grnd
    this%qflx_h2osfc_surf_col => col_wf%qflx_h2osfc_surf
    this%qflx_gross_infl_soil_col => col_wf%qflx_gross_infl_soil
    this%qflx_adv_col => col_wf%qflx_adv
    this%qflx_rootsoi_col => col_wf%qflx_rootsoi
    this%grnd_ch4_cond_patch => ch4_vars%grnd_ch4_cond_patch
    this%btran_patch => energyflux_vars%btran_patch
    this%btran2_patch => energyflux_vars%btran2_patch
    this%rresis_patch => energyflux_vars%rresis_patch
    this%canopy_cond_patch => energyflux_vars%canopy_cond_patch
    this%eff_porosity_col => soilstate_vars%eff_porosity_col
    this%rootr_patch => soilstate_vars%rootr_patch
    this%h2osoi_liqvol => col_ws%h2osoi_liqvol
    this%cgrnds => veg_ef%cgrnds
    this%cgrndl => veg_ef%cgrndl
    this%dlrad => veg_ef%dlrad
    this%ulrad => veg_ef%ulrad
    this%cgrnd => veg_ef%cgrnd
    this%eflx_sh_snow => veg_ef%eflx_sh_snow
    this%eflx_sh_h2osfc => veg_ef%eflx_sh_h2osfc
    this%eflx_sh_soil => veg_ef%eflx_sh_soil
    this%eflx_sh_veg => veg_ef%eflx_sh_veg
    this%eflx_sh_grnd => veg_ef%eflx_sh_grnd
    this%eflx_sh_tot => veg_ef%eflx_sh_tot
    this%eflx_hs_soil_col => col_ef%eflx_hs_soil
    this%eflx_hs_top_snow_col => col_ef%eflx_hs_top_snow
    this%eflx_hs_h2osfc_col => col_ef%eflx_hs_h2osfc
    this%eflx_dhsdT_col => col_ef%eflx_dhsdT
    this%rssun_patch => photosyns_vars%rssun_patch
    this%rssha_patch => photosyns_vars%rssha_patch
    this%lmrsun_patch => photosyns_vars%lmrsun_patch
    this%lmrsha_patch => photosyns_vars%lmrsha_patch
    this%laisun_patch => canopystate_vars%laisun_patch
    this%laisha_patch => canopystate_vars%laisha_patch
    this%fpsn_patch => photosyns_vars%fpsn_patch
    this%psnsun_patch => photosyns_vars%psnsun_patch
    this%psnsha_patch => photosyns_vars%psnsha_patch
    this%t_veg => veg_es%t_veg
    this%t_ref2m => veg_es%t_ref2m

  end subroutine bind_canopyflux_emulator_output_view

  subroutine assemble_canopyflux_emulator_features(bounds, num_nolakeurbanp, filter_nolakeurbanp, &
     inputs, features, flnr_override, mbbopt_override, pco2_override)

    type(bounds_type)                         , intent(in)  :: bounds
    integer                                   , intent(in)  :: num_nolakeurbanp
    integer                                   , intent(in)  :: filter_nolakeurbanp(:)
    type(canopyflux_emulator_input_view_type) , intent(in)  :: inputs
    real(r8)                                  , intent(out) :: features(:,:)
    real(r8), intent(in), optional            :: flnr_override(:)
    real(r8), intent(in), optional            :: mbbopt_override(:)
   real(r8), intent(in), optional            :: pco2_override(:)

    integer :: fp, p, c, g, t, lev, idx, ivt

    if (size(features, 1) /= num_nolakeurbanp) then
       error stop 'assemble_canopyflux_emulator_features: unexpected sample dimension'
    end if

    if (size(features, 2) /= canopyflux_emulator_num_features) then
       error stop 'assemble_canopyflux_emulator_features: unexpected feature dimension'
    end if

    do fp = 1, num_nolakeurbanp
       p = filter_nolakeurbanp(fp)
       c = veg_pp%column(p)
       g = veg_pp%gridcell(p)
       t = veg_pp%topounit(p)
       ivt = veg_pp%itype(p)

       ! Feature ordering is part of the emulator contract.
       features(fp,  1) = inputs%forc_t(t)                ! near-surface air temperature [K]
       features(fp,  2) = inputs%forc_q(t)                ! near-surface specific humidity [kg kg-1]
       features(fp,  3) = inputs%forc_pbot(t)             ! surface pressure [Pa]
       features(fp,  4) = inputs%forc_lwrad(t)            ! downward longwave radiation [W m-2]
       features(fp,  5) = inputs%forc_rain(t)             ! rainfall rate [mm s-1]
       features(fp,  6) = inputs%forc_snow(t)             ! snowfall rate [mm s-1]
       features(fp,  7) = inputs%forc_swrad(t)            ! total downward shortwave radiation [W m-2]
       features(fp,  8) = inputs%forc_solad(t, 1)         ! direct visible shortwave radiation [W m-2]
       features(fp,  9) = inputs%forc_solad(t, 2)         ! direct NIR shortwave radiation [W m-2]
       features(fp, 10) = inputs%forc_solai(t, 1)         ! diffuse visible shortwave radiation [W m-2]
       features(fp, 11) = inputs%forc_solai(t, 2)         ! diffuse NIR shortwave radiation [W m-2]
       features(fp, 12) = patch_wind_speed(inputs%forc_u(t), inputs%forc_v(t)) ! horizontal wind speed [m s-1]
       features(fp, 13) = inputs%forc_rho(t)              ! air density [kg m-3]
       if (present(pco2_override)) then
          features(fp, 14) = pco2_override(fp)            ! canopy CO2 partial pressure [Pa]
       else
          features(fp, 14) = inputs%forc_pco2(t)          ! canopy CO2 partial pressure [Pa]
       end if
       features(fp, 15) = inputs%frac_veg_nosno_patch(p)  ! exposed vegetation fraction [-]
       features(fp, 16) = inputs%elai_patch(p)            ! snow-buried effective leaf area index [m2 m-2]
       features(fp, 17) = inputs%esai_patch(p)            ! snow-buried effective stem area index [m2 m-2]
       features(fp, 18) = inputs%frac_sno(c)              ! ground snow cover fraction [-]
       features(fp, 19) = inputs%snow_depth(c)            ! snow depth [m]
       features(fp, 20) = inputs%zwt_col(c)               ! water table depth [m]
       features(fp, 21) = inputs%h2ocan(p)                ! canopy water storage [mm]
       features(fp, 22) = inputs%t10(p)                   ! 10-day mean 2 m air temperature for acclimation [K]
       features(fp, 23) = inputs%dayl(g)                  ! current daylength [s]
       features(fp, 24) = inputs%max_dayl(g)              ! annual maximum daylength [s]
       features(fp, 25) = inputs%forc_hgt_u_patch(p)      ! forcing wind reference height [m]
       features(fp, 26) = inputs%z0mg_col(c)              ! ground roughness length for momentum [m]
       features(fp, 27) = inputs%frac_h2osfc(c)           ! surface water fraction [-]
       features(fp, 28) = inputs%coszen_col(c)            ! cosine solar zenith angle [-]

       ! Top-N soil matric potential profile [mm].
       idx = 29
       do lev = 1, canopyflux_emulator_profile_layers
          features(fp, idx + lev - 1) = inputs%smp_l_col(c, lev)
       end do

       ! Top-N effective porosity profile [m3 m-3].
       idx = idx + canopyflux_emulator_profile_layers
       do lev = 1, canopyflux_emulator_profile_layers
          features(fp, idx + lev - 1) = inputs%eff_porosity_col(c, lev)
       end do

       ! Top-N volumetric liquid water content profile [m3 m-3].
       idx = idx + canopyflux_emulator_profile_layers
       do lev = 1, canopyflux_emulator_profile_layers
          features(fp, idx + lev - 1) = inputs%h2osoi_liqvol(c, lev)
       end do

       ! Top-N root fraction profile for this patch [-].
       idx = idx + canopyflux_emulator_profile_layers
       do lev = 1, canopyflux_emulator_rootfr_layers
          features(fp, idx + lev - 1) = inputs%rootfr_patch(p, lev)
       end do

       idx = idx + canopyflux_emulator_rootfr_layers
       features(fp, idx    ) = inputs%t_grnd(c)               ! ground surface temperature [K]
       features(fp, idx + 1) = inputs%t_h2osfc(c)             ! surface-water temperature [K]
       features(fp, idx + 2) = inputs%ugust(t)                ! atmospheric gustiness [m s-1]

       ! Top-N soil/snow temperature profile [K].
       idx = idx + 3
       do lev = 1, canopyflux_emulator_profile_layers
          features(fp, idx + lev - 1) = inputs%t_soisno(c, lev)
       end do

       idx = idx + canopyflux_emulator_profile_layers
       if (present(flnr_override)) then
          features(fp, idx    ) = flnr_override(fp)
       else
          features(fp, idx    ) = veg_vp%flnr(ivt)
       end if
       if (present(mbbopt_override)) then
          features(fp, idx + 1) = mbbopt_override(fp)
       else
          features(fp, idx + 1) = veg_vp%mbbopt(ivt)
       end if
       features(fp, idx + 2) = inputs%sabg_snow_patch(p)  ! solar radiation absorbed by snow [W m-2]
    end do

  end subroutine assemble_canopyflux_emulator_features

  subroutine sample_canopyflux_training_traits(flnr_sample, mbbopt_sample)

    real(r8), intent(out) :: flnr_sample
    real(r8), intent(out) :: mbbopt_sample

    real(r8) :: draw

    call random_number(draw)
    flnr_sample = canopyflux_training_flnr_min + &
         (canopyflux_training_flnr_max - canopyflux_training_flnr_min) * draw
    call random_number(draw)
    mbbopt_sample = canopyflux_training_mbbopt_min + &
         (canopyflux_training_mbbopt_max - canopyflux_training_mbbopt_min) * draw

  end subroutine sample_canopyflux_training_traits

   subroutine sample_canopyflux_training_pco2(pco2_sample)

      real(r8), intent(out) :: pco2_sample

      real(r8) :: draw

      call random_number(draw)
      pco2_sample = canopyflux_training_pco2_min + &
             (canopyflux_training_pco2_max - canopyflux_training_pco2_min) * draw

   end subroutine sample_canopyflux_training_pco2

  subroutine assemble_canopyflux_emulator_targets(num_nolakeurbanp, filter_nolakeurbanp, inputs, outputs, targets)

    integer                                    , intent(in)  :: num_nolakeurbanp
    integer                                    , intent(in)  :: filter_nolakeurbanp(:)
    type(canopyflux_emulator_input_view_type)  , intent(in)  :: inputs
    type(canopyflux_emulator_output_view_type) , intent(in)  :: outputs
    real(r8)                                   , intent(out) :: targets(:,:)

    integer :: fp, p, c

    if (size(targets, 1) /= num_nolakeurbanp) then
       error stop 'assemble_canopyflux_emulator_targets: unexpected sample dimension'
    end if

    if (size(targets, 2) /= canopyflux_emulator_num_targets) then
       error stop 'assemble_canopyflux_emulator_targets: unexpected target dimension'
    end if

    do fp = 1, num_nolakeurbanp
       p = filter_nolakeurbanp(fp)
       c = veg_pp%column(p)
       ! Target ordering is part of the emulator contract.
       targets(fp,  1) = outputs%fpsn_patch(p)            ! total photosynthesis [umol CO2 m-2 s-1]
       targets(fp,  2) = outputs%eflx_sh_tot(p)           ! total sensible heat flux [W m-2]
       targets(fp,  3) = outputs%qflx_tran_veg(p)         ! transpiration flux [mm s-1]
       targets(fp,  4) = max(0._r8, outputs%qflx_evap_veg(p) - outputs%qflx_tran_veg(p)) ! canopy-water evaporation [mm s-1]
       targets(fp,  5) = outputs%qflx_ev_soil(p)          ! soil evaporation flux [mm s-1]
       targets(fp,  6) = outputs%qflx_ev_snow(p)          ! snow evaporation/sublimation flux [mm s-1]
       targets(fp,  7) = outputs%qflx_ev_h2osfc(p)        ! surface-water evaporation flux [mm s-1]
       targets(fp,  8) = outputs%qflx_prec_grnd(p)        ! water onto ground including canopy runoff [mm s-1]
       targets(fp,  9) = outputs%t_veg(p)                 ! vegetation temperature [K]
       targets(fp, 10) = outputs%t_ref2m(p)               ! 2 m reference air temperature [K]
       targets(fp, 11) = outputs%eflx_hs_soil_col(c)      ! implicit soil-surface heat flux [W m-2]
       targets(fp, 12) = outputs%eflx_hs_top_snow_col(c)  ! implicit snow-surface heat flux [W m-2]
       targets(fp, 13) = outputs%eflx_hs_h2osfc_col(c)    ! implicit standing-water heat flux [W m-2]
       targets(fp, 14) = outputs%eflx_dhsdT_col(c)        ! deriv. of energy flux wrt surface temp [W m-2 K-1]
       targets(fp, 15) = inputs%laisun_patch(p)           ! sunlit projected leaf area index [m2 m-2]
       targets(fp, 16) = inputs%laisha_patch(p)           ! shaded projected leaf area index [m2 m-2]
       targets(fp, 17) = outputs%psnsun_patch(p)          ! sunlit leaf photosynthesis [umol CO2 m-2 s-1]
       targets(fp, 18) = outputs%psnsha_patch(p)          ! shaded leaf photosynthesis [umol CO2 m-2 s-1]
    end do

  end subroutine assemble_canopyflux_emulator_targets

  subroutine assemble_canopyflux_emulator_params(num_nolakeurbanp, filter_nolakeurbanp, params)

    integer  , intent(in)  :: num_nolakeurbanp
    integer  , intent(in)  :: filter_nolakeurbanp(:)
    real(r8) , intent(out) :: params(:,:)

    integer :: fp, p, ivt

    if (size(params, 1) /= num_nolakeurbanp) then
       error stop 'assemble_canopyflux_emulator_params: unexpected sample dimension'
    end if

    if (size(params, 2) /= canopyflux_emulator_num_params) then
       error stop 'assemble_canopyflux_emulator_params: unexpected parameter dimension'
    end if

    do fp = 1, num_nolakeurbanp
       p = filter_nolakeurbanp(fp)
       ivt = veg_pp%itype(p)

       params(fp,  1) = real(ivt, r8)
       params(fp,  2) = veg_vp%c3psn(ivt)
       params(fp,  3) = veg_vp%nonvascular(ivt)
       params(fp,  4) = veg_vp%dleaf(ivt)
       params(fp,  5) = veg_vp%smpso(ivt)
       params(fp,  6) = veg_vp%smpsc(ivt)
       params(fp,  7) = veg_vp%slatop(ivt)
       params(fp,  8) = veg_vp%leafcn(ivt)
       params(fp,  9) = veg_vp%flnr(ivt)
       params(fp, 10) = veg_vp%fnitr(ivt)
       params(fp, 11) = 0._r8  ! br_mr_pft not available in this configuration
       params(fp, 12) = 0._r8  ! q10_mr_pft not available in this configuration
       params(fp, 13) = veg_vp%i_vc(ivt)
       params(fp, 14) = veg_vp%s_vc(ivt)
       params(fp, 15) = veg_vp%act25(ivt)
       params(fp, 16) = veg_vp%kcha(ivt)
       params(fp, 17) = veg_vp%koha(ivt)
       params(fp, 18) = veg_vp%cpha(ivt)
       params(fp, 19) = veg_vp%vcmaxha(ivt)
       params(fp, 20) = veg_vp%jmaxha(ivt)
       params(fp, 21) = veg_vp%tpuha(ivt)
       params(fp, 22) = veg_vp%lmrha(ivt)
       params(fp, 23) = veg_vp%vcmaxhd(ivt)
       params(fp, 24) = veg_vp%jmaxhd(ivt)
       params(fp, 25) = veg_vp%tpuhd(ivt)
       params(fp, 26) = veg_vp%lmrhd(ivt)
       params(fp, 27) = veg_vp%lmrse(ivt)
       params(fp, 28) = veg_vp%qe(ivt)
       params(fp, 29) = veg_vp%theta_cj(ivt)
       params(fp, 30) = veg_vp%bbbopt(ivt)
       params(fp, 31) = veg_vp%mbbopt(ivt)
    end do

  end subroutine assemble_canopyflux_emulator_params

  subroutine begin_canopyflux_training_capture(bounds_proc)

    type(bounds_type), intent(in) :: bounds_proc
    integer :: num_local_patches

   if (.not. canopyflux_training_capture_enabled()) return

    captured_begp = bounds_proc%begp
    captured_endp = bounds_proc%endp
    num_local_patches = captured_endp - captured_begp + 1

    if (.not. allocated(captured_features)) then
       allocate(captured_features(num_local_patches, canopyflux_emulator_num_features))
       allocate(captured_targets(num_local_patches, canopyflux_emulator_num_targets))
       allocate(captured_debug(num_local_patches, canopyflux_emulator_num_debug))
       allocate(captured_natveg_mask(num_local_patches))
    else if (size(captured_features, 1) /= num_local_patches) then
       deallocate(captured_features, captured_targets, captured_debug, captured_natveg_mask)
       allocate(captured_features(num_local_patches, canopyflux_emulator_num_features))
       allocate(captured_targets(num_local_patches, canopyflux_emulator_num_targets))
       allocate(captured_debug(num_local_patches, canopyflux_emulator_num_debug))
       allocate(captured_natveg_mask(num_local_patches))
    end if

    captured_features(:,:) = 0._r8
    captured_targets(:,:) = 0._r8
    captured_debug(:,:) = 0._r8
    captured_natveg_mask(:) = .false.
    capture_ready = .true.

  end subroutine begin_canopyflux_training_capture

  subroutine capture_canopyflux_training_features(num_nolakeurbanp, filter_nolakeurbanp, features, btran)

    integer , intent(in) :: num_nolakeurbanp
    integer , intent(in) :: filter_nolakeurbanp(:)
    real(r8), intent(in) :: features(:,:)
    real(r8), intent(in) :: btran(:)

    integer :: fp, p, slot

    if (.not. capture_ready) return

    do fp = 1, num_nolakeurbanp
       p = filter_nolakeurbanp(fp)
       if (.not. is_natveg_patch(p)) cycle
       slot = capture_slot(p)
       captured_features(slot,:) = features(fp,:)
       captured_debug(slot,1) = btran(fp)
       captured_natveg_mask(slot) = .true.
    end do

  end subroutine capture_canopyflux_training_features

  subroutine capture_canopyflux_training_debug(num_nolakeurbanp, filter_nolakeurbanp)

    integer , intent(in) :: num_nolakeurbanp
    integer , intent(in) :: filter_nolakeurbanp(:)

    integer :: fp, p, c, slot

    if (.not. capture_ready) return

   do fp = 1, num_nolakeurbanp
       p = filter_nolakeurbanp(fp)
       if (.not. is_natveg_patch(p)) cycle
       c = veg_pp%column(p)
       slot = capture_slot(p)
       captured_debug(slot,2) = col_wf%qflx_infl(c)
       captured_debug(slot,3) = col_ws%frac_h2osfc_act(c)
       captured_debug(slot,4) = col_ws%h2osfc(c)
       captured_debug(slot,5) = col_wf%qflx_h2osfc_surf(c)
       captured_debug(slot,6) = col_wf%qflx_gross_infl_soil(c)
       captured_debug(slot,7) = col_ws%frac_h2osfc(c)
       captured_debug(slot,8) = col_wf%qflx_adv(c,1)
       captured_debug(slot,9) = col_wf%qflx_infl(c) - col_wf%qflx_rootsoi(c,1) - col_wf%qflx_adv(c,1)
       captured_debug(slot,10) = col_wf%qflx_drain(c)
       captured_debug(slot,11) = col_wf%qflx_rsub_sat(c)
       captured_debug(slot,12) = col_wf%qflx_rootsoi(c,1)
       captured_debug(slot,13) = col_wf%qflx_rootsoi(c,2)
       captured_debug(slot,14) = col_wf%qflx_rootsoi(c,3)
       captured_debug(slot,15) = col_wf%qflx_rootsoi(c,4)
       captured_debug(slot,16) = col_wf%qflx_rootsoi(c,5)
       captured_debug(slot,17) = col_wf%qflx_rootsoi(c,6)
       captured_debug(slot,18) = col_wf%qflx_rootsoi(c,7)
       captured_debug(slot,19) = col_wf%qflx_rootsoi(c,8)
       captured_debug(slot,20) = col_wf%qflx_rootsoi(c,9)
       captured_debug(slot,21) = col_wf%qflx_rootsoi(c,10)
       captured_debug(slot,22) = col_ws%h2osno(c)
       captured_debug(slot,23) = col_wf%qflx_snow_grnd(c)
    end do

  end subroutine capture_canopyflux_training_debug

  subroutine capture_canopyflux_training_targets(num_nolakeurbanp, filter_nolakeurbanp, targets)

    integer , intent(in) :: num_nolakeurbanp
    integer , intent(in) :: filter_nolakeurbanp(:)
    real(r8), intent(in) :: targets(:,:)

    integer :: fp, p, slot

    if (.not. capture_ready) return

    do fp = 1, num_nolakeurbanp
      p = filter_nolakeurbanp(fp)
      if (.not. is_natveg_patch(p)) cycle
      slot = capture_slot(p)
      captured_targets(slot,:) = targets(fp,:)
    end do

  end subroutine capture_canopyflux_training_targets

  subroutine finalize_canopyflux_training_capture(bounds, num_nolakeurbanp, filter_nolakeurbanp, &
       atm2lnd_vars, canopystate_vars, energyflux_vars, frictionvel_vars, soilstate_vars, solarabs_vars, surfalb_vars, &
       ch4_vars, photosyns_vars, soilhydrology_vars)

    type(bounds_type)         , intent(in)    :: bounds
    integer                   , intent(in)    :: num_nolakeurbanp
    integer                   , intent(in)    :: filter_nolakeurbanp(:)
    type(atm2lnd_type)        , intent(inout) :: atm2lnd_vars
    type(canopystate_type)    , intent(inout) :: canopystate_vars
    type(energyflux_type)     , intent(inout) :: energyflux_vars
    type(frictionvel_type)    , intent(inout) :: frictionvel_vars
    type(soilstate_type)      , intent(inout) :: soilstate_vars
    type(solarabs_type)       , intent(inout) :: solarabs_vars
    type(surfalb_type)        , intent(inout) :: surfalb_vars
    type(ch4_type)            , intent(inout) :: ch4_vars
    type(photosyns_type)      , intent(inout) :: photosyns_vars
    type(soilhydrology_type)  , intent(in)    :: soilhydrology_vars
    type(canopyflux_emulator_input_view_type)  :: inputs
    type(canopyflux_emulator_output_view_type) :: outputs
    real(r8) :: targets(num_nolakeurbanp, canopyflux_emulator_num_targets)

   if (.not. canopyflux_training_capture_enabled()) return
    if (.not. capture_ready) return
    if (num_nolakeurbanp <= 0) return

    call inputs%bind(atm2lnd_vars, canopystate_vars, energyflux_vars, frictionvel_vars, &
         soilstate_vars, solarabs_vars, surfalb_vars, soilhydrology_vars)
    call outputs%bind(canopystate_vars, energyflux_vars, frictionvel_vars, &
         soilstate_vars, ch4_vars, photosyns_vars)
    call assemble_canopyflux_emulator_targets(num_nolakeurbanp, filter_nolakeurbanp, inputs, outputs, targets)
    call capture_canopyflux_training_targets(num_nolakeurbanp, filter_nolakeurbanp, targets)

  end subroutine finalize_canopyflux_training_capture

  subroutine flush_canopyflux_training_capture()

    real(r8), allocatable :: local_patch_ids(:)
    real(r8), allocatable :: local_cell_ids(:)
    real(r8), allocatable :: local_features_flat(:)
    real(r8), allocatable :: local_targets_flat(:)
    real(r8), allocatable :: local_debug_flat(:)
    real(r8), allocatable :: local_params_flat(:)
    real(r8), pointer :: global_patch_ids(:) => null()
    real(r8), pointer :: global_cell_ids(:) => null()
    real(r8), pointer :: global_features_flat(:) => null()
    real(r8), pointer :: global_targets_flat(:) => null()
    real(r8), pointer :: global_debug_flat(:) => null()
    real(r8), pointer :: global_params_flat(:) => null()
    real(r8), allocatable :: global_features(:,:)
    real(r8), allocatable :: global_targets(:,:)
    real(r8), allocatable :: global_debug(:,:)
    real(r8), allocatable :: global_params(:,:)
    integer :: local_natveg_count
    integer :: global_natveg_count
    integer :: slot, sample, p
    integer :: year, mon, day, sec, nstep
    integer, pointer :: patch_counts(:) => null()
    integer, pointer :: patch_displs(:) => null()
    integer, pointer :: cell_counts(:) => null()
    integer, pointer :: cell_displs(:) => null()
    integer, pointer :: feature_counts(:) => null()
    integer, pointer :: feature_displs(:) => null()
    integer, pointer :: target_counts(:) => null()
    integer, pointer :: target_displs(:) => null()
    integer, pointer :: debug_counts(:) => null()
    integer, pointer :: debug_displs(:) => null()
    integer, pointer :: param_counts(:) => null()
    integer, pointer :: param_displs(:) => null()
    real(r8), allocatable :: local_params(:,:)

   if (.not. canopyflux_training_capture_enabled()) return
    if (.not. capture_ready) return

    local_natveg_count = count(captured_natveg_mask)
    allocate(local_params(size(captured_features,1), canopyflux_emulator_num_params))
    local_params(:,:) = 0._r8
    allocate(local_patch_ids(local_natveg_count))
    allocate(local_cell_ids(local_natveg_count))
    allocate(local_features_flat(local_natveg_count * canopyflux_emulator_num_features))
    allocate(local_targets_flat(local_natveg_count * canopyflux_emulator_num_targets))
    allocate(local_debug_flat(local_natveg_count * canopyflux_emulator_num_debug))
    allocate(local_params_flat(local_natveg_count * canopyflux_emulator_num_params))

    call assemble_canopyflux_emulator_params(size(captured_features,1), [(captured_begp + slot - 1, slot=1,size(captured_features,1))], local_params)

    sample = 0
    do slot = 1, size(captured_natveg_mask)
       if (.not. captured_natveg_mask(slot)) cycle
       sample = sample + 1
       p = captured_begp + slot - 1
       local_patch_ids(sample) = real(veg_pp%itype(p), r8)
       local_cell_ids(sample)  = real(p, r8)
       local_features_flat((sample-1)*canopyflux_emulator_num_features + 1 : sample*canopyflux_emulator_num_features) = &
            captured_features(slot,:)
       local_targets_flat((sample-1)*canopyflux_emulator_num_targets + 1 : sample*canopyflux_emulator_num_targets) = &
            captured_targets(slot,:)
       local_debug_flat((sample-1)*canopyflux_emulator_num_debug + 1 : sample*canopyflux_emulator_num_debug) = &
            captured_debug(slot,:)
       local_params_flat((sample-1)*canopyflux_emulator_num_params + 1 : sample*canopyflux_emulator_num_params) = &
            local_params(slot,:)
    end do

    call shr_mpi_gathScatvInit(mpicom, 0, local_patch_ids, global_patch_ids, patch_counts, patch_displs)
    call shr_mpi_gatherv(local_patch_ids, size(local_patch_ids), global_patch_ids, patch_counts, patch_displs, 0, mpicom)

    call shr_mpi_gathScatvInit(mpicom, 0, local_cell_ids, global_cell_ids, cell_counts, cell_displs)
    call shr_mpi_gatherv(local_cell_ids, size(local_cell_ids), global_cell_ids, cell_counts, cell_displs, 0, mpicom)

    call shr_mpi_gathScatvInit(mpicom, 0, local_features_flat, global_features_flat, feature_counts, feature_displs)
    call shr_mpi_gatherv(local_features_flat, size(local_features_flat), global_features_flat, feature_counts, feature_displs, 0, mpicom)

    call shr_mpi_gathScatvInit(mpicom, 0, local_targets_flat, global_targets_flat, target_counts, target_displs)
    call shr_mpi_gatherv(local_targets_flat, size(local_targets_flat), global_targets_flat, target_counts, target_displs, 0, mpicom)

    call shr_mpi_gathScatvInit(mpicom, 0, local_debug_flat, global_debug_flat, debug_counts, debug_displs)
    call shr_mpi_gatherv(local_debug_flat, size(local_debug_flat), global_debug_flat, debug_counts, debug_displs, 0, mpicom)

    call shr_mpi_gathScatvInit(mpicom, 0, local_params_flat, global_params_flat, param_counts, param_displs)
    call shr_mpi_gatherv(local_params_flat, size(local_params_flat), global_params_flat, param_counts, param_displs, 0, mpicom)

    if (masterproc) then
       global_natveg_count = size(global_patch_ids)
       if (global_natveg_count == 0) then
          capture_ready = .false.
          deallocate(global_patch_ids, global_features_flat, global_targets_flat, global_debug_flat, global_params_flat)
          if (associated(patch_counts)) deallocate(patch_counts, patch_displs)
          if (associated(feature_counts)) deallocate(feature_counts, feature_displs)
          if (associated(target_counts)) deallocate(target_counts, target_displs)
          if (associated(debug_counts)) deallocate(debug_counts, debug_displs)
          if (associated(param_counts)) deallocate(param_counts, param_displs)
          deallocate(local_patch_ids, local_features_flat, local_targets_flat, local_debug_flat, local_params_flat, local_params)
          return
       end if
       allocate(global_features(global_natveg_count, canopyflux_emulator_num_features))
       allocate(global_targets(global_natveg_count, canopyflux_emulator_num_targets))
       allocate(global_debug(global_natveg_count, canopyflux_emulator_num_debug))
       allocate(global_params(global_natveg_count, canopyflux_emulator_num_params))
       do sample = 1, global_natveg_count
          global_features(sample,:) = global_features_flat((sample-1)*canopyflux_emulator_num_features + 1 : &
               sample*canopyflux_emulator_num_features)
          global_targets(sample,:) = global_targets_flat((sample-1)*canopyflux_emulator_num_targets + 1 : &
               sample*canopyflux_emulator_num_targets)
          global_debug(sample,:) = global_debug_flat((sample-1)*canopyflux_emulator_num_debug + 1 : &
               sample*canopyflux_emulator_num_debug)
          global_params(sample,:) = global_params_flat((sample-1)*canopyflux_emulator_num_params + 1 : &
               sample*canopyflux_emulator_num_params)
       end do
       call get_curr_date(year, mon, day, sec)
       nstep = get_nstep()
       call write_canopyflux_training_netcdf(global_patch_ids, global_cell_ids, global_features, global_targets, global_debug, global_params, year, mon, day, sec, nstep)
       deallocate(global_features, global_targets, global_debug, global_params)
    end if

    if (associated(patch_counts)) deallocate(patch_counts, patch_displs)
    if (associated(cell_counts)) deallocate(cell_counts, cell_displs)
    if (associated(feature_counts)) deallocate(feature_counts, feature_displs)
    if (associated(target_counts)) deallocate(target_counts, target_displs)
    if (associated(debug_counts)) deallocate(debug_counts, debug_displs)
    if (associated(param_counts)) deallocate(param_counts, param_displs)
    if (associated(global_patch_ids)) deallocate(global_patch_ids)
    if (associated(global_cell_ids)) deallocate(global_cell_ids)
    if (associated(global_features_flat)) deallocate(global_features_flat)
    if (associated(global_targets_flat)) deallocate(global_targets_flat)
    if (associated(global_debug_flat)) deallocate(global_debug_flat)
    if (associated(global_params_flat)) deallocate(global_params_flat)
    deallocate(local_patch_ids, local_cell_ids, local_features_flat, local_targets_flat, local_debug_flat, local_params_flat, local_params)

    capture_ready = .false.

  end subroutine flush_canopyflux_training_capture

  subroutine write_canopyflux_training_netcdf(patch_ids_r8, cell_ids_r8, features, targets, debug, params, year, mon, day, sec, nstep)

    real(r8), intent(in) :: patch_ids_r8(:)
    real(r8), intent(in) :: cell_ids_r8(:)
    real(r8), intent(in) :: features(:,:)
    real(r8), intent(in) :: targets(:,:)
    real(r8), intent(in) :: debug(:,:)
    real(r8), intent(in) :: params(:,:)
    integer , intent(in) :: year, mon, day, sec, nstep

    integer :: ncid, dim_sample, dim_feature, dim_target, dim_debug, dim_param
    integer :: var_patch, var_cell, var_nstep, var_year, var_month, var_day, var_sec
    integer :: var_features, var_targets, var_debug, var_params
    integer :: sample_start
    integer :: sample_count
    integer :: file_feature_count, file_target_count, file_debug_count, file_param_count
    integer :: status
    integer, allocatable :: patch_ids(:)
    integer, allocatable :: cell_ids(:)
    integer, allocatable :: nstep_vec(:), year_vec(:), month_vec(:), day_vec(:), sec_vec(:)
    real(r8), allocatable :: features_out(:,:), targets_out(:,:), debug_out(:,:), params_out(:,:)
    character(len=256) :: file_name
    logical :: file_exists

    allocate(patch_ids(size(patch_ids_r8)))
    patch_ids(:) = int(patch_ids_r8(:))
    allocate(cell_ids(size(cell_ids_r8)))
    cell_ids(:) = int(cell_ids_r8(:))
    sample_count = size(patch_ids)
    allocate(nstep_vec(sample_count), year_vec(sample_count), month_vec(sample_count), day_vec(sample_count), sec_vec(sample_count))
    allocate(features_out(canopyflux_emulator_num_features, sample_count))
    allocate(targets_out(canopyflux_emulator_num_targets, sample_count))
    allocate(debug_out(canopyflux_emulator_num_debug, sample_count))
    allocate(params_out(canopyflux_emulator_num_params, sample_count))
    file_debug_count = -1

    nstep_vec(:) = nstep
    year_vec(:) = year
    month_vec(:) = mon
    day_vec(:) = day
    sec_vec(:) = sec
    features_out(:,:) = transpose(features)
    targets_out(:,:) = transpose(targets)
    debug_out(:,:) = transpose(debug)
    params_out(:,:) = transpose(params)

    write(file_name,'(a,".canopyflux_training.nc")') trim(caseid)
    inquire(file=trim(file_name), exist=file_exists)
    if (.not. canopyflux_training_file_initialized) file_exists = .false.

    if (file_exists) then
       call netcdf_check(nf90_open(trim(file_name), nf90_write, ncid), 'open '//trim(file_name))
       call netcdf_check(nf90_inq_dimid(ncid, 'feature', dim_feature), 'inq_dimid feature')
       call netcdf_check(nf90_inquire_dimension(ncid, dim_feature, len=file_feature_count), 'inquire_dimension feature')
       call netcdf_check(nf90_inq_dimid(ncid, 'target', dim_target), 'inq_dimid target')
       call netcdf_check(nf90_inquire_dimension(ncid, dim_target, len=file_target_count), 'inquire_dimension target')
       status = nf90_inq_dimid(ncid, 'debug', dim_debug)
       if (status == nf90_noerr) then
          call netcdf_check(nf90_inquire_dimension(ncid, dim_debug, len=file_debug_count), 'inquire_dimension debug')
       else
          file_exists = .false.
       end if
       ! Require cell_index variable; recreate the file if it is absent (old format).
       status = nf90_inq_varid(ncid, 'cell_index', var_cell)
       if (status /= nf90_noerr) file_exists = .false.
       call netcdf_check(nf90_inq_dimid(ncid, 'param', dim_param), 'inq_dimid param')
       call netcdf_check(nf90_inquire_dimension(ncid, dim_param, len=file_param_count), 'inquire_dimension param')
       call netcdf_check(nf90_close(ncid), 'close')

       if (file_feature_count /= canopyflux_emulator_num_features .or. &
            file_target_count /= canopyflux_emulator_num_targets .or. &
            (file_exists .and. file_debug_count /= canopyflux_emulator_num_debug) .or. &
            file_param_count /= canopyflux_emulator_num_params) then
          file_exists = .false.
       end if
    end if

    if (.not. file_exists) then
       call netcdf_check(nf90_create(trim(file_name), nf90_clobber, ncid), 'create '//trim(file_name))
       call netcdf_check(nf90_put_att(ncid, nf90_global, 'title', 'ELM canopy-flux training data'), 'put_att title')
       call netcdf_check(nf90_put_att(ncid, nf90_global, 'case', trim(caseid)), 'put_att case')
       call netcdf_check(nf90_put_att(ncid, nf90_global, 'patch_filter', 'natural vegetation only (istsoil)'), 'put_att patch_filter')
       call netcdf_check(nf90_put_att(ncid, nf90_global, 'patch_index_meaning', 'veg_pp%itype (PFT index), not raw patch id'), 'put_att patch_index_meaning')
       call netcdf_check(nf90_put_att(ncid, nf90_global, 'feature_names', canopyflux_feature_names), 'put_att feature_names')
       call netcdf_check(nf90_put_att(ncid, nf90_global, 'target_names', canopyflux_target_names), 'put_att target_names')
       call netcdf_check(nf90_put_att(ncid, nf90_global, 'debug_names', &
            'btran,qflx_infl_col,frac_h2osfc_act,h2osfc,qflx_h2osfc_surf,' // &
            'qflx_gross_infl_soil,frac_h2osfc,q12_col,dW1dt_col,qflx_drain_col,' // &
            'qflx_rsub_sat_col,qflx_rootsoi,qflx_rootsoi_lev2,qflx_rootsoi_lev3,' // &
            'qflx_rootsoi_lev4,qflx_rootsoi_lev5,qflx_rootsoi_lev6,qflx_rootsoi_lev7,' // &
            'qflx_rootsoi_lev8,qflx_rootsoi_lev9,qflx_rootsoi_lev10,h2osno,' // &
            'qflx_snow_grnd'), 'put_att debug_names')
       call netcdf_check(nf90_put_att(ncid, nf90_global, 'param_names', canopyflux_param_names), 'put_att param_names')

       call netcdf_check(nf90_def_dim(ncid, 'sample', nf90_unlimited, dim_sample), 'def_dim sample')
       call netcdf_check(nf90_def_dim(ncid, 'feature', canopyflux_emulator_num_features, dim_feature), 'def_dim feature')
       call netcdf_check(nf90_def_dim(ncid, 'target', canopyflux_emulator_num_targets, dim_target), 'def_dim target')
       call netcdf_check(nf90_def_dim(ncid, 'debug', canopyflux_emulator_num_debug, dim_debug), 'def_dim debug')
       call netcdf_check(nf90_def_dim(ncid, 'param', canopyflux_emulator_num_params, dim_param), 'def_dim param')

       call netcdf_check(nf90_def_var(ncid, 'patch_index', nf90_int, (/dim_sample/), var_patch), 'def_var patch_index')
       call netcdf_check(nf90_def_var(ncid, 'cell_index',  nf90_int, (/dim_sample/), var_cell),  'def_var cell_index')
       call netcdf_check(nf90_def_var(ncid, 'nstep', nf90_int, (/dim_sample/), var_nstep), 'def_var nstep')
       call netcdf_check(nf90_def_var(ncid, 'year', nf90_int, (/dim_sample/), var_year), 'def_var year')
       call netcdf_check(nf90_def_var(ncid, 'month', nf90_int, (/dim_sample/), var_month), 'def_var month')
       call netcdf_check(nf90_def_var(ncid, 'day', nf90_int, (/dim_sample/), var_day), 'def_var day')
       call netcdf_check(nf90_def_var(ncid, 'sec', nf90_int, (/dim_sample/), var_sec), 'def_var sec')
       call netcdf_check(nf90_def_var(ncid, 'features', nf90_double, (/dim_feature, dim_sample/), var_features), 'def_var features')
       call netcdf_check(nf90_def_var(ncid, 'targets', nf90_double, (/dim_target, dim_sample/), var_targets), 'def_var targets')
       call netcdf_check(nf90_def_var(ncid, 'debug', nf90_double, (/dim_debug, dim_sample/), var_debug), 'def_var debug')
       call netcdf_check(nf90_def_var(ncid, 'parameters', nf90_double, (/dim_param, dim_sample/), var_params), 'def_var parameters')
       call netcdf_check(nf90_enddef(ncid), 'enddef')
       sample_start = 1
    else
       call netcdf_check(nf90_open(trim(file_name), nf90_write, ncid), 'open '//trim(file_name))
       call netcdf_check(nf90_inq_dimid(ncid, 'sample', dim_sample), 'inq_dimid sample')
       call netcdf_check(nf90_inquire_dimension(ncid, dim_sample, len=sample_start), 'inquire_dimension sample')
       sample_start = sample_start + 1
       call netcdf_check(nf90_inq_varid(ncid, 'patch_index', var_patch), 'inq_varid patch_index')
       call netcdf_check(nf90_inq_varid(ncid, 'cell_index',  var_cell),  'inq_varid cell_index')
       call netcdf_check(nf90_inq_varid(ncid, 'nstep', var_nstep), 'inq_varid nstep')
       call netcdf_check(nf90_inq_varid(ncid, 'year', var_year), 'inq_varid year')
       call netcdf_check(nf90_inq_varid(ncid, 'month', var_month), 'inq_varid month')
       call netcdf_check(nf90_inq_varid(ncid, 'day', var_day), 'inq_varid day')
       call netcdf_check(nf90_inq_varid(ncid, 'sec', var_sec), 'inq_varid sec')
       call netcdf_check(nf90_inq_varid(ncid, 'features', var_features), 'inq_varid features')
       call netcdf_check(nf90_inq_varid(ncid, 'targets', var_targets), 'inq_varid targets')
       call netcdf_check(nf90_inq_varid(ncid, 'debug', var_debug), 'inq_varid debug')
       call netcdf_check(nf90_inq_varid(ncid, 'parameters', var_params), 'inq_varid parameters')
    end if

    call netcdf_check(nf90_put_var(ncid, var_patch, patch_ids, start=(/sample_start/), count=(/sample_count/)), 'put_var patch_index')
    call netcdf_check(nf90_put_var(ncid, var_cell,  cell_ids,  start=(/sample_start/), count=(/sample_count/)), 'put_var cell_index')
    call netcdf_check(nf90_put_var(ncid, var_nstep, nstep_vec, start=(/sample_start/), count=(/sample_count/)), 'put_var nstep')
    call netcdf_check(nf90_put_var(ncid, var_year, year_vec, start=(/sample_start/), count=(/sample_count/)), 'put_var year')
    call netcdf_check(nf90_put_var(ncid, var_month, month_vec, start=(/sample_start/), count=(/sample_count/)), 'put_var month')
    call netcdf_check(nf90_put_var(ncid, var_day, day_vec, start=(/sample_start/), count=(/sample_count/)), 'put_var day')
    call netcdf_check(nf90_put_var(ncid, var_sec, sec_vec, start=(/sample_start/), count=(/sample_count/)), 'put_var sec')
    call netcdf_check(nf90_put_var(ncid, var_features, features_out, start=(/1, sample_start/), &
         count=(/canopyflux_emulator_num_features, sample_count/)), 'put_var features')
    call netcdf_check(nf90_put_var(ncid, var_targets, targets_out, start=(/1, sample_start/), &
         count=(/canopyflux_emulator_num_targets, sample_count/)), 'put_var targets')
    call netcdf_check(nf90_put_var(ncid, var_debug, debug_out, start=(/1, sample_start/), &
         count=(/canopyflux_emulator_num_debug, sample_count/)), 'put_var debug')
    call netcdf_check(nf90_put_var(ncid, var_params, params_out, start=(/1, sample_start/), &
         count=(/canopyflux_emulator_num_params, sample_count/)), 'put_var parameters')
    call netcdf_check(nf90_close(ncid), 'close')
    canopyflux_training_file_initialized = .true.

    deallocate(patch_ids, cell_ids, nstep_vec, year_vec, month_vec, day_vec, sec_vec, features_out, targets_out, debug_out, params_out)

  end subroutine write_canopyflux_training_netcdf

  subroutine netcdf_check(status, action)

    integer         , intent(in) :: status
    character(len=*), intent(in) :: action

    if (status /= nf90_noerr) then
       error stop 'CanopyFluxesEmulatorMod netcdf error during '//trim(action)//': '//trim(nf90_strerror(status))
    end if

  end subroutine netcdf_check

  elemental pure real(r8) function gelu_activation(x) result(y)

    real(r8), intent(in) :: x

    y = 0.5_r8 * x * (1._r8 + tanh(0.7978845608028654_r8 * (x + 0.044715_r8 * x**3)))

  end function gelu_activation

  subroutine split_csv_names(csv_text, names)

    character(len=*), intent(in) :: csv_text
    character(len=64), allocatable, intent(out) :: names(:)

    integer :: token_count
    integer :: start_pos
    integer :: end_pos
    integer :: idx
    character(len=:), allocatable :: token

    token_count = 1
    do idx = 1, len_trim(csv_text)
      if (csv_text(idx:idx) == ',') token_count = token_count + 1
    end do

    allocate(names(token_count))
    names(:) = ''

    start_pos = 1
    idx = 0
    do
       end_pos = index(csv_text(start_pos:), ',')
       idx = idx + 1
       if (end_pos == 0) then
          token = adjustl(trim(csv_text(start_pos:)))
          names(idx) = token
          exit
       else
          token = adjustl(trim(csv_text(start_pos:start_pos + end_pos - 2)))
          names(idx) = token
          start_pos = start_pos + end_pos
       end if
    end do

  end subroutine split_csv_names

  integer pure function find_name_index(name, name_list) result(idx)

    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: name_list(:)
    integer :: item

    idx = 0
    do item = 1, size(name_list)
       if (trim(name) == trim(name_list(item))) then
          idx = item
          return
       end if
    end do

  end function find_name_index

  logical pure function is_log1p_target(canonical_idx) result(use_log1p)

    integer, intent(in) :: canonical_idx

    select case (canonical_idx)
    case (9, 23)
       use_log1p = .true.
    case default
       use_log1p = .false.
    end select

  end function is_log1p_target

  subroutine clear_canopyflux_model(model)

    type(canopyflux_model_type), intent(inout) :: model

    if (allocated(model%input_map)) deallocate(model%input_map)
    if (allocated(model%target_map)) deallocate(model%target_map)
    if (allocated(model%nonnegative_target)) deallocate(model%nonnegative_target)
    if (allocated(model%x_mean)) deallocate(model%x_mean)
    if (allocated(model%x_std)) deallocate(model%x_std)
    if (allocated(model%y_mean)) deallocate(model%y_mean)
    if (allocated(model%y_std)) deallocate(model%y_std)
    if (allocated(model%layers)) deallocate(model%layers)
    model%is_loaded = .false.
    model%patch_index = -1

  end subroutine clear_canopyflux_model

  subroutine canopyflux_model_path_for_pft(ivt, file_name)

    integer         , intent(in)  :: ivt
    character(len=*), intent(out) :: file_name

    write(file_name, '(a,"/pft_",i0,".nc")') trim(canopyflux_model_dir), ivt

  end subroutine canopyflux_model_path_for_pft

  subroutine log_required_canopyflux_pfts(filter_nolakeurbanp)

    integer, intent(in) :: filter_nolakeurbanp(:)

    logical :: seen_pft(numpft)
    integer :: fp, p, ivt

    if (reported_required_pfts) return
    if (.not. masterproc) return

    seen_pft(:) = .false.
    do fp = 1, size(filter_nolakeurbanp)
       p = filter_nolakeurbanp(fp)
       if (p < 1 .or. p > size(veg_pp%itype)) cycle
       ivt = veg_pp%itype(p)
       if (ivt >= 1 .and. ivt <= numpft) seen_pft(ivt) = .true.
    end do

    write(*,*) 'CanopyFluxesEmulator runtime PFT indices on master task:'
    do ivt = 1, numpft
       if (seen_pft(ivt)) write(*,*) '  ivt =', ivt
    end do
    reported_required_pfts = .true.

  end subroutine log_required_canopyflux_pfts

  subroutine ensure_canopyflux_model_loaded(ivt)

    integer, intent(in) :: ivt
    character(len=512) :: file_name
    logical :: file_exists

    if (ivt < 1 .or. ivt > numpft) then
       error stop 'CanopyFluxesEmulatorMod: invalid pft index for emulator model'
    end if

    if (.not. allocated(canopyflux_models)) then
       allocate(canopyflux_models(numpft))
    end if

    if (canopyflux_models(ivt)%is_loaded) return
    if (len_trim(canopyflux_model_dir) == 0) then
       error stop 'CanopyFluxesEmulatorMod: canopyflux_model_dir is empty'
    end if

    call canopyflux_model_path_for_pft(ivt, file_name)
    inquire(file=trim(file_name), exist=file_exists)
    if (.not. file_exists) then
       error stop 'CanopyFluxesEmulatorMod: no emulator model file found for requested PFT'
    end if
    call load_canopyflux_model(trim(file_name), canopyflux_models(ivt))

  end subroutine ensure_canopyflux_model_loaded

  subroutine load_canopyflux_model(file_name, model)

    character(len=*), intent(in) :: file_name
    type(canopyflux_model_type), intent(inout) :: model

    integer :: ncid
    integer :: varid
    integer :: input_dimid, output_dimid
    integer :: input_count, output_count
    integer :: bias_dimids(1)
    integer :: layer_idx
    integer :: status
    integer :: out_count, in_count
    integer :: prev_width
    real(r8), allocatable :: weight_tmp(:,:)
    integer :: canonical_idx
    character(len=8192) :: input_names_text
    character(len=8192) :: target_names_text
    character(len=64), allocatable :: canonical_input_names(:)
    character(len=64), allocatable :: canonical_target_names(:)
    character(len=64), allocatable :: file_input_names(:)
    character(len=64), allocatable :: file_target_names(:)
    character(len=32) :: patch_index_text
    character(len=32) :: weight_name
    character(len=32) :: bias_name

    call clear_canopyflux_model(model)

    call split_csv_names(canopyflux_feature_names, canonical_input_names)
    call split_csv_names(canopyflux_target_names, canonical_target_names)

    call netcdf_check(nf90_open(trim(file_name), nf90_nowrite, ncid), 'open '//trim(file_name))
    call netcdf_check(nf90_inq_dimid(ncid, 'input', input_dimid), 'inq_dimid input')
    call netcdf_check(nf90_inquire_dimension(ncid, input_dimid, len=input_count), 'inquire_dimension input')
    call netcdf_check(nf90_inq_dimid(ncid, 'output', output_dimid), 'inq_dimid output')
    call netcdf_check(nf90_inquire_dimension(ncid, output_dimid, len=output_count), 'inquire_dimension output')

    input_names_text = ''
    target_names_text = ''
    patch_index_text = ''
    call netcdf_check(nf90_get_att(ncid, nf90_global, 'input_names', input_names_text), 'get_att input_names')
    call netcdf_check(nf90_get_att(ncid, nf90_global, 'target_names', target_names_text), 'get_att target_names')
    status = nf90_get_att(ncid, nf90_global, 'patch_index', model%patch_index)
    if (status /= nf90_noerr) then
       status = nf90_get_att(ncid, nf90_global, 'patch_index', patch_index_text)
       if (status == nf90_noerr) read(patch_index_text, *) model%patch_index
    end if

    call split_csv_names(trim(input_names_text), file_input_names)
    call split_csv_names(trim(target_names_text), file_target_names)

    if (size(file_input_names) /= input_count) then
       error stop 'CanopyFluxesEmulatorMod: input_names length does not match model input dimension'
    end if
    if (size(file_target_names) /= output_count) then
       error stop 'CanopyFluxesEmulatorMod: target_names length does not match model output dimension'
    end if
    if (input_count /= size(canonical_input_names)) then
       error stop 'CanopyFluxesEmulatorMod: model input count does not match current emulator contract'
    end if
    allocate(model%input_map(input_count), model%x_mean(input_count), model%x_std(input_count))
    allocate(model%target_map(output_count), model%y_mean(output_count), model%y_std(output_count))
    allocate(model%nonnegative_target(output_count))

    do canonical_idx = 1, input_count
       model%input_map(canonical_idx) = find_name_index(file_input_names(canonical_idx), canonical_input_names)
       if (model%input_map(canonical_idx) == 0) then
          error stop 'CanopyFluxesEmulatorMod: unknown input feature in model file'
       end if
    end do

    do canonical_idx = 1, output_count
       model%target_map(canonical_idx) = find_name_index(file_target_names(canonical_idx), canonical_target_names)
       if (model%target_map(canonical_idx) == 0) then
          error stop 'CanopyFluxesEmulatorMod: unknown target field in model file'
       end if
    end do

    if (output_count /= 12 .and. output_count /= 14 .and. output_count /= 18) then
        error stop 'CanopyFluxesEmulatorMod: model output count does not match supported emulator contracts'
    end if

    call netcdf_check(nf90_inq_varid(ncid, 'x_mean', varid), 'inq_varid x_mean')
    call netcdf_check(nf90_get_var(ncid, varid, model%x_mean), 'get_var x_mean')
    call netcdf_check(nf90_inq_varid(ncid, 'x_std', varid), 'inq_varid x_std')
    call netcdf_check(nf90_get_var(ncid, varid, model%x_std), 'get_var x_std')
    call netcdf_check(nf90_inq_varid(ncid, 'y_mean', varid), 'inq_varid y_mean')
    call netcdf_check(nf90_get_var(ncid, varid, model%y_mean), 'get_var y_mean')
    call netcdf_check(nf90_inq_varid(ncid, 'y_std', varid), 'inq_varid y_std')
    call netcdf_check(nf90_get_var(ncid, varid, model%y_std), 'get_var y_std')
    call netcdf_check(nf90_inq_varid(ncid, 'nonnegative_target', varid), 'inq_varid nonnegative_target')
    call netcdf_check(nf90_get_var(ncid, varid, model%nonnegative_target), 'get_var nonnegative_target')

    layer_idx = 0
    do
       write(weight_name, '("weight_", i0)') layer_idx
       status = nf90_inq_varid(ncid, trim(weight_name), varid)
       if (status /= nf90_noerr) exit
       layer_idx = layer_idx + 1
    end do

    if (layer_idx == 0) then
       error stop 'CanopyFluxesEmulatorMod: no dense layers found in model file'
    end if

    allocate(model%layers(layer_idx))
    do layer_idx = 1, size(model%layers)
       write(weight_name, '("weight_", i0)') layer_idx - 1
       write(bias_name  , '("bias_", i0)') layer_idx - 1
       call netcdf_check(nf90_inq_varid(ncid, trim(bias_name), varid), 'inq_varid '//trim(bias_name))
       call netcdf_check(nf90_inquire_variable(ncid, varid, dimids=bias_dimids), 'inquire_variable '//trim(bias_name))
       call netcdf_check(nf90_inquire_dimension(ncid, bias_dimids(1), len=out_count), 'inquire_dimension '//trim(bias_name))
       allocate(model%layers(layer_idx)%bias(out_count))
       call netcdf_check(nf90_get_var(ncid, varid, model%layers(layer_idx)%bias), 'get_var '//trim(bias_name))

       if (layer_idx == 1) then
          prev_width = size(model%x_mean)
       else
          prev_width = size(model%layers(layer_idx - 1)%bias)
       end if
       in_count = prev_width

       call netcdf_check(nf90_inq_varid(ncid, trim(weight_name), varid), 'inq_varid '//trim(weight_name))
       allocate(model%layers(layer_idx)%weight(out_count, in_count))
       allocate(weight_tmp(in_count, out_count))
       call netcdf_check(nf90_get_var(ncid, varid, weight_tmp), 'get_var '//trim(weight_name))
       model%layers(layer_idx)%weight(:,:) = transpose(weight_tmp)
       deallocate(weight_tmp)
    end do

    call netcdf_check(nf90_close(ncid), 'close '//trim(file_name))
    model%is_loaded = .true.

  end subroutine load_canopyflux_model

  subroutine run_canopyflux_model(model, features, predictions)

    type(canopyflux_model_type), intent(in) :: model
    real(r8), intent(in) :: features(:,:)
    real(r8), intent(out) :: predictions(:,:)

    integer :: layer_idx
    integer :: output_idx
    integer :: input_idx
    integer :: sample_count
    real(r8), allocatable :: layer_input(:,:)
    real(r8), allocatable :: layer_output(:,:)

    if (.not. model%is_loaded) then
       error stop 'CanopyFluxesEmulatorMod: model must be loaded before inference'
    end if

    if (size(predictions, 1) /= size(features, 1)) then
       error stop 'CanopyFluxesEmulatorMod: prediction sample dimension mismatch'
    end if
    if (size(predictions, 2) /= size(model%target_map)) then
       error stop 'CanopyFluxesEmulatorMod: prediction target dimension mismatch'
    end if

    sample_count = size(features, 1)
    allocate(layer_input(sample_count, size(model%input_map)))
    do input_idx = 1, size(model%input_map)
       layer_input(:, input_idx) = (features(:, model%input_map(input_idx)) - model%x_mean(input_idx)) / &
            model%x_std(input_idx)
    end do

    do layer_idx = 1, size(model%layers)
       allocate(layer_output(sample_count, size(model%layers(layer_idx)%bias)))
       layer_output(:,:) = matmul(layer_input, transpose(model%layers(layer_idx)%weight))
       do output_idx = 1, size(model%layers(layer_idx)%bias)
          layer_output(:, output_idx) = layer_output(:, output_idx) + model%layers(layer_idx)%bias(output_idx)
       end do
       if (layer_idx < size(model%layers)) then
          layer_output(:,:) = gelu_activation(layer_output)
       end if
       call move_alloc(layer_output, layer_input)
    end do

    predictions(:,:) = layer_input(:,:)
    do output_idx = 1, size(model%y_std)
       predictions(:, output_idx) = predictions(:, output_idx) * model%y_std(output_idx) + model%y_mean(output_idx)
       if (is_log1p_target(model%target_map(output_idx))) then
          predictions(:, output_idx) = exp(predictions(:, output_idx)) - 1._r8
       end if
       if (model%nonnegative_target(output_idx) /= 0) then
          predictions(:, output_idx) = max(predictions(:, output_idx), 0._r8)
       end if
    end do

    if (allocated(layer_input)) deallocate(layer_input)

  end subroutine run_canopyflux_model

  subroutine apply_canopyflux_predictions(filter_nolakeurbanp, inputs, outputs, model, predictions)

    integer, intent(in) :: filter_nolakeurbanp(:)
    type(canopyflux_emulator_input_view_type), intent(in) :: inputs
    type(canopyflux_emulator_output_view_type), intent(inout) :: outputs
    type(canopyflux_model_type), intent(in) :: model
    real(r8), intent(in) :: predictions(:,:)

    integer :: fp, p, c, t, target_idx, canonical_idx
    real(r8) :: dtime
    real(r8) :: h2ocan_old
    real(r8) :: h2ocan_pred
    real(r8) :: qflx_tran_pred
    real(r8) :: qflx_evap_canopy_pred
    real(r8) :: qflx_evap_canopy_max
    real(r8) :: weighted_component_sum
    real(r8) :: snow_frac
    real(r8) :: h2osfc_frac
    real(r8) :: soil_frac
    real(r8) :: t_veg_delta_pred
    real(r8) :: tlbef
    real(r8) :: dt_veg
    real(r8) :: lw_grnd
    real(r8) :: cf_bare
    real(r8) :: canopy_latent_flux
    real(r8) :: canopy_air_lw
    real(r8) :: canopy_self_lw
    real(r8) :: canopy_ground_lw

    dtime = real(get_step_size(), r8)
    if (dtime <= 0._r8) then
       error stop 'CanopyFluxesEmulatorMod: nonpositive timestep in apply_canopyflux_predictions'
    end if

    do fp = 1, size(filter_nolakeurbanp)
       p = filter_nolakeurbanp(fp)
       c = veg_pp%column(p)
       t = veg_pp%topounit(p)
       h2ocan_old = outputs%h2ocan(p)
       h2ocan_pred = h2ocan_old
       qflx_tran_pred = outputs%qflx_tran_veg(p)
       qflx_evap_canopy_pred = max(0._r8, outputs%qflx_evap_veg(p) - outputs%qflx_tran_veg(p))
       t_veg_delta_pred = 0._r8
       do target_idx = 1, size(model%target_map)
          canonical_idx = model%target_map(target_idx)
          select case (canonical_idx)
          case (1)
             h2ocan_pred = predictions(fp, target_idx)
          case (2)
             outputs%n_irrig_steps_left(p) = max(0, nint(predictions(fp, target_idx)))
          case (3)
             outputs%irrig_rate(p) = predictions(fp, target_idx)
          case (4)
             qflx_tran_pred = predictions(fp, target_idx)
          case (5)
             qflx_evap_canopy_pred = predictions(fp, target_idx)
          case (6)
             outputs%qflx_ev_snow(p) = predictions(fp, target_idx)
          case (7)
             outputs%qflx_ev_soil(p) = predictions(fp, target_idx)
          case (8)
             outputs%qflx_ev_h2osfc(p) = predictions(fp, target_idx)
          case (9)
             outputs%ram1_patch(p) = canopyflux_ram1_from_conductance(predictions(fp, target_idx))
          case (10)
             outputs%rb1_patch(p) = predictions(fp, target_idx)
          case (11)
             outputs%grnd_ch4_cond_patch(p) = predictions(fp, target_idx)
          case (12)
             outputs%canopy_cond_patch(p) = predictions(fp, target_idx)
          case (13)
             outputs%cgrnds(p) = predictions(fp, target_idx)
          case (14)
             outputs%cgrndl(p) = predictions(fp, target_idx)
          case (15)
             outputs%laisun_patch(p) = max(0._r8, predictions(fp, target_idx))
          case (16)
             outputs%laisha_patch(p) = max(0._r8, predictions(fp, target_idx))
          case (17)
             outputs%psnsun_patch(p) = predictions(fp, target_idx)
          case (18)
             outputs%psnsha_patch(p) = predictions(fp, target_idx)
          case (19)
             outputs%rssun_patch(p) = predictions(fp, target_idx)
          case (20)
             outputs%rssha_patch(p) = predictions(fp, target_idx)
          case (21)
             outputs%lmrsun_patch(p) = predictions(fp, target_idx)
          case (22)
             outputs%lmrsha_patch(p) = predictions(fp, target_idx)
          case (23)
             outputs%psnsun_patch(p) = predictions(fp, target_idx)
          case (24)
             outputs%psnsha_patch(p) = predictions(fp, target_idx)
          case (25)
             t_veg_delta_pred = predictions(fp, target_idx)
          case (26)
             outputs%t_ref2m(p) = predictions(fp, target_idx)
          case default
             error stop 'CanopyFluxesEmulatorMod: unsupported target index in model file'
          end select
       end do

       call debug_print_canopyflux_psn_state('raw', p, inputs, outputs)

       ! Ground-flux partitions are support-limited by snow, bare-soil, and
       ! surface-water area. Rebuild total ground evaporation directly from the
       ! physically active component fluxes.
       snow_frac = max(0._r8, min(1._r8, inputs%frac_sno(c)))
       h2osfc_frac = max(0._r8, min(1._r8, inputs%frac_h2osfc(c)))
       soil_frac = max(0._r8, 1._r8 - snow_frac - h2osfc_frac)

       if (snow_frac <= 0._r8) outputs%qflx_ev_snow(p) = 0._r8
       if (soil_frac <= 0._r8) outputs%qflx_ev_soil(p) = 0._r8
       if (h2osfc_frac <= 0._r8) outputs%qflx_ev_h2osfc(p) = 0._r8

       weighted_component_sum = snow_frac * outputs%qflx_ev_snow(p) + &
            soil_frac * outputs%qflx_ev_soil(p) + &
            h2osfc_frac * outputs%qflx_ev_h2osfc(p)
       outputs%qflx_evap_soi(p) = weighted_component_sum

       ! Mirror native canopy-water guards: transpiration cannot be negative,
       ! canopy-water evaporation cannot exceed available canopy water,
       ! and total vegetated evaporation is transpiration plus canopy-water evaporation.
       outputs%qflx_tran_veg(p) = max(0._r8, qflx_tran_pred)
       if (inputs%elai_patch(p) <= canopyflux_emulator_elai_native_fallback .or. &
            inputs%frac_veg_nosno_patch(p) <= canopyflux_emulator_frac_veg_nosno_native_fallback) then
          outputs%qflx_tran_veg(p) = 0._r8
       end if
       qflx_evap_canopy_max = max(0._r8, h2ocan_old) / dtime
       qflx_evap_canopy_pred = max(0._r8, qflx_evap_canopy_pred)
       outputs%qflx_evap_veg(p) = outputs%qflx_tran_veg(p) + min(qflx_evap_canopy_pred, qflx_evap_canopy_max)
       cf_bare = inputs%forc_pbot(t) / (SHR_CONST_RGAS * 0.001_r8 * inputs%thm(c)) * 1.e06_r8
       outputs%rssun_patch(p) = max(outputs%rssun_patch(p), 1._r8 / 1.e15_r8 * cf_bare)
       outputs%rssha_patch(p) = max(outputs%rssha_patch(p), 1._r8 / 1.e15_r8 * cf_bare)
       if (inputs%sabv_patch(p) <= tiny(1._r8)) then
          outputs%psnsun_patch(p) = 0._r8
          outputs%psnsha_patch(p) = 0._r8
       end if
       if (inputs%elai_patch(p) <= canopyflux_emulator_elai_fpsn_cutoff .or. &
            inputs%frac_veg_nosno_patch(p) <= canopyflux_emulator_frac_veg_nosno_native_fallback) then
          outputs%psnsun_patch(p) = 0._r8
          outputs%psnsha_patch(p) = 0._r8
       end if
       outputs%fpsn_patch(p) = outputs%psnsun_patch(p) * outputs%laisun_patch(p) + &
            outputs%psnsha_patch(p) * outputs%laisha_patch(p)
       outputs%h2ocan(p) = max(0._r8, max(0._r8, h2ocan_old) + &
            (outputs%qflx_tran_veg(p) - outputs%qflx_evap_veg(p)) * dtime)
       outputs%cgrnd(p) = outputs%cgrnds(p) + outputs%cgrndl(p) * col_ef%htvp(c)
       outputs%t_veg(p) = inputs%forc_t(t) + t_veg_delta_pred

       ! Mirror native longwave diagnostics from the updated canopy temperature
       ! instead of trusting them as free ML targets.
       tlbef = inputs%t_veg(p)
       dt_veg = outputs%t_veg(p) - tlbef
       lw_grnd = snow_frac * inputs%t_soisno(c, col_pp%snl(c) + 1)**4 + &
            soil_frac * inputs%t_soisno(c, 1)**4 + &
            h2osfc_frac * inputs%t_h2osfc(c)**4
       outputs%dlrad(p) = (1._r8 - inputs%emv(p)) * inputs%emg(c) * inputs%forc_lwrad(t) + &
            inputs%emv(p) * inputs%emg(c) * sb * tlbef**3 * (tlbef + 4._r8 * dt_veg)
       outputs%ulrad(p) = ((1._r8 - inputs%emg(c)) * (1._r8 - inputs%emv(p))**2 * inputs%forc_lwrad(t)) + &
            inputs%emv(p) * (1._r8 + (1._r8 - inputs%emg(c)) * (1._r8 - inputs%emv(p))) * &
            sb * tlbef**3 * (tlbef + 4._r8 * dt_veg) + &
            inputs%emg(c) * (1._r8 - inputs%emv(p)) * sb * lw_grnd

       canopy_air_lw = inputs%emv(p) * (1._r8 + (1._r8 - inputs%emv(p)) * (1._r8 - inputs%emg(c))) * &
            inputs%forc_lwrad(t)
       canopy_self_lw = - (2._r8 - inputs%emv(p) * (1._r8 - inputs%emg(c))) * inputs%emv(p) * sb * &
            outputs%t_veg(p)**4
       canopy_ground_lw = inputs%emv(p) * inputs%emg(c) * sb * lw_grnd
       canopy_latent_flux = hvap * outputs%qflx_evap_veg(p)
       outputs%eflx_sh_veg(p) = inputs%sabv_patch(p) + canopy_air_lw + canopy_self_lw + &
            canopy_ground_lw - canopy_latent_flux

       call debug_print_canopyflux_longwave_state('emulator', p, inputs, outputs)
       call debug_print_canopyflux_psn_state('postproc', p, inputs, outputs)
    end do

  end subroutine apply_canopyflux_predictions

  subroutine rebuild_canopyflux_derived_state(bounds, num_nolakeurbanp, filter_nolakeurbanp, &
       inputs, outputs, canopystate_vars, energyflux_vars, soilstate_vars)

    type(bounds_type)                            , intent(in)    :: bounds
    integer                                      , intent(in)    :: num_nolakeurbanp
    integer                                      , intent(in)    :: filter_nolakeurbanp(:)
    type(canopyflux_emulator_input_view_type)    , intent(in)    :: inputs
    type(canopyflux_emulator_output_view_type)   , intent(inout) :: outputs
    type(canopystate_type)                       , intent(in)    :: canopystate_vars
    type(energyflux_type)                        , intent(inout) :: energyflux_vars
    type(soilstate_type)                         , intent(inout) :: soilstate_vars

    integer, allocatable :: filterc(:)
    integer, allocatable :: jtop(:)
    integer :: fp, p, c, num_filterc

    do fp = 1, num_nolakeurbanp
       p = filter_nolakeurbanp(fp)
       outputs%btran_patch(p) = 0._r8
       outputs%btran2_patch(p) = 0._r8
       outputs%rootr_patch(p,:) = 0._r8
       outputs%rresis_patch(p,:) = 0._r8
    end do

    allocate(filterc(num_nolakeurbanp))
    num_filterc = 0
    do fp = 1, num_nolakeurbanp
       p = filter_nolakeurbanp(fp)
       if (p < 1 .or. p > size(veg_pp%column)) cycle
       c = veg_pp%column(p)
       if (num_filterc == 0 .or. .not. any(filterc(1:num_filterc) == c)) then
          num_filterc = num_filterc + 1
          filterc(num_filterc) = c
       end if
    end do

    if (num_filterc > 0) then
       call calc_effective_soilporosity(bounds,                          &
            ubj = nlevgrnd,                                              &
            numf = num_filterc,                                          &
            filter = filterc(1:num_filterc),                             &
            watsat = inputs%watsat_col(bounds%begc:bounds%endc, 1:nlevgrnd), &
            h2osoi_ice = inputs%h2osoi_ice(bounds%begc:bounds%endc, 1:nlevgrnd), &
            denice = denice,                                             &
            eff_por = outputs%eff_porosity_col(bounds%begc:bounds%endc, 1:nlevgrnd))

       allocate(jtop(bounds%begc:bounds%endc))
       jtop(bounds%begc:bounds%endc) = 1
       call calc_volumetric_h2oliq(bounds,                                    &
            jtop = jtop(bounds%begc:bounds%endc),                             &
            lbj = 1,                                                          &
            ubj = nlevgrnd,                                                   &
            numf = num_filterc,                                               &
            filter = filterc(1:num_filterc),                                  &
            eff_porosity = outputs%eff_porosity_col(bounds%begc:bounds%endc, 1:nlevgrnd), &
            h2osoi_liq = inputs%h2osoi_liq(bounds%begc:bounds%endc, 1:nlevgrnd), &
            denh2o = denh2o,                                                  &
            vol_liq = outputs%h2osoi_liqvol(bounds%begc:bounds%endc, 1:nlevgrnd))
       deallocate(jtop)
    end if

    call set_perchroot_opt(perchroot, perchroot_alt)
    call calc_root_moist_stress(bounds,     &
         nlevgrnd = nlevgrnd,               &
         fn = num_nolakeurbanp,             &
         filterp = filter_nolakeurbanp,     &
         canopystate_vars = canopystate_vars, &
         energyflux_vars = energyflux_vars, &
         soilstate_vars = soilstate_vars)

    deallocate(filterc)

  end subroutine rebuild_canopyflux_derived_state

  subroutine debug_print_canopyflux_predictions(filter_nolakeurbanp, model, predictions)

    integer, intent(in) :: filter_nolakeurbanp(:)
    type(canopyflux_model_type), intent(in) :: model
    real(r8), intent(in) :: predictions(:,:)

    integer :: fp, p, target_idx, canonical_idx

    if (.not. debug_canopyflux_emulator) return
    if (reported_emulator_outputs) return
    if (.not. masterproc) return

    do fp = 1, size(filter_nolakeurbanp)
       p = filter_nolakeurbanp(fp)
       if (p < 1 .or. p > size(veg_pp%itype)) cycle
       write(*,*) 'CanopyFluxesEmulator predicted outputs for first active patch:'
       write(*,*) '  patch =', p, ' ivt =', veg_pp%itype(p)
       do target_idx = 1, size(model%target_map)
          canonical_idx = model%target_map(target_idx)
          if (canonical_idx <= 0) cycle
          write(*,'(a,i0,a,es24.16)') '  target(', canonical_idx, ') = ', predictions(fp, target_idx)
       end do
       reported_emulator_outputs = .true.
       exit
    end do

  end subroutine debug_print_canopyflux_predictions

  subroutine canopyflux_info3330_path(file_name)

    character(len=*), intent(out) :: file_name

    character(len=512) :: base_dir
    integer :: slash_pos

    base_dir = trim(canopyflux_model_dir)
    if (len_trim(base_dir) > 0) then
       if (base_dir(len_trim(base_dir):len_trim(base_dir)) == '/') then
          base_dir = base_dir(:len_trim(base_dir) - 1)
       end if
    end if
    slash_pos = scan(trim(base_dir), '/', back=.true.)
    if (slash_pos <= 0) then
       file_name = 'info3330.txt'
    else
       file_name = trim(base_dir(:slash_pos - 1)) // '/info3330.txt'
    end if

  end subroutine canopyflux_info3330_path

  subroutine reset_info3330_file_if_needed()

    integer :: unitno
    integer :: nstep
    character(len=512) :: file_name

    if (.not. masterproc) return
    nstep = get_nstep()
    if (nstep /= 3330) return
    if (info3330_reset_nstep == nstep) return

    call canopyflux_info3330_path(file_name)
    open(newunit=unitno, file=trim(file_name), status='replace', action='write')
    close(unitno)
    info3330_reset_nstep = nstep

  end subroutine reset_info3330_file_if_needed

  subroutine debug_dump_canopyflux_features(filter_nolakeurbanp, features)

    integer , intent(in) :: filter_nolakeurbanp(:)
    real(r8), intent(in) :: features(:,:)

    character(len=64), allocatable :: feature_names(:)
    character(len=512) :: file_name
    integer :: fp, p, feature_idx
    integer :: nstep
    integer :: unitno

    nstep = get_nstep()
    if (nstep /= 3330) return
    if (.not. masterproc) return

    call reset_info3330_file_if_needed()
    call canopyflux_info3330_path(file_name)
    call split_csv_names(canopyflux_feature_names, feature_names)

    open(newunit=unitno, file=trim(file_name), status='old', action='write', position='append')
    do fp = 1, size(filter_nolakeurbanp)
       p = filter_nolakeurbanp(fp)
       if (p /= 2 .and. p /= 8) cycle
       write(unitno,*) 'CanopyFluxesEmulator feature dump: nstep=', nstep, ' patch=', p, ' ivt=', veg_pp%itype(p)
       do feature_idx = 1, size(feature_names)
          write(unitno,'(a,a,a,es24.16)') '  ', trim(feature_names(feature_idx)), ' = ', features(fp, feature_idx)
       end do
    end do
    close(unitno)

    deallocate(feature_names)

  end subroutine debug_dump_canopyflux_features

  real(r8) function canopyflux_output_value(inputs, outputs, p, canonical_idx) result(value)

    type(canopyflux_emulator_input_view_type), intent(in) :: inputs
    type(canopyflux_emulator_output_view_type), intent(in) :: outputs
    integer, intent(in) :: p
    integer, intent(in) :: canonical_idx
    integer :: t

    t = veg_pp%topounit(p)

    select case (canonical_idx)
    case (1)
       value = outputs%h2ocan(p)
    case (2)
       value = real(outputs%n_irrig_steps_left(p), r8)
    case (3)
       value = outputs%irrig_rate(p)
    case (4)
       value = outputs%qflx_tran_veg(p)
    case (5)
       value = max(0._r8, outputs%qflx_evap_veg(p) - outputs%qflx_tran_veg(p))
    case (6)
       value = outputs%qflx_ev_snow(p)
    case (7)
       value = outputs%qflx_ev_soil(p)
    case (8)
       value = outputs%qflx_ev_h2osfc(p)
    case (9)
       value = canopyflux_conductance_from_ram1(outputs%ram1_patch(p))
    case (10)
       value = outputs%rb1_patch(p)
    case (11)
       value = outputs%grnd_ch4_cond_patch(p)
    case (12)
       value = outputs%canopy_cond_patch(p)
    case (13)
       value = outputs%cgrnds(p)
    case (14)
       value = outputs%cgrndl(p)
    case (15)
       value = outputs%laisun_patch(p)
    case (16)
       value = outputs%laisha_patch(p)
    case (17)
       value = outputs%psnsun_patch(p)
    case (18)
       value = outputs%psnsha_patch(p)
    case (19)
       value = outputs%rssun_patch(p)
    case (20)
       value = outputs%rssha_patch(p)
    case (21)
       value = outputs%lmrsun_patch(p)
    case (22)
       value = outputs%lmrsha_patch(p)
    case (23)
       value = outputs%psnsun_patch(p)
    case (24)
       value = outputs%psnsha_patch(p)
    case (25)
       value = outputs%t_veg(p) - inputs%forc_t(t)
    case (26)
       value = outputs%t_ref2m(p)
    case default
       value = 0._r8
    end select

  end function canopyflux_output_value

  subroutine debug_dump_canopyflux_targets(filter_nolakeurbanp, inputs, outputs, model, predictions)

    integer, intent(in) :: filter_nolakeurbanp(:)
    type(canopyflux_emulator_input_view_type), intent(in) :: inputs
    type(canopyflux_emulator_output_view_type), intent(in) :: outputs
    type(canopyflux_model_type), intent(in) :: model
    real(r8), intent(in) :: predictions(:,:)

    character(len=64), allocatable :: target_names(:)
    character(len=512) :: file_name
    integer :: fp, p, target_idx, canonical_idx
    integer :: unitno

    if (get_nstep() /= 3330) return
    if (.not. masterproc) return

    call reset_info3330_file_if_needed()
    call canopyflux_info3330_path(file_name)
    call split_csv_names(canopyflux_target_names, target_names)

    open(newunit=unitno, file=trim(file_name), status='old', action='write', position='append')
    do fp = 1, size(filter_nolakeurbanp)
       p = filter_nolakeurbanp(fp)
       if (p /= 2 .and. p /= 8) cycle
       write(unitno,*) 'CanopyFluxesEmulator target dump: nstep=', get_nstep(), ' patch=', p, ' ivt=', veg_pp%itype(p)
       do target_idx = 1, size(model%target_map)
          canonical_idx = model%target_map(target_idx)
          write(unitno,'(a,a,a,es24.16,a,es24.16)') '  ', trim(target_names(canonical_idx)), &
               ' raw = ', predictions(fp, target_idx), ' final = ', canopyflux_output_value(inputs, outputs, p, canonical_idx)
       end do
       write(unitno,'(a,es24.16)') '  qflx_evap_soi_final = ', outputs%qflx_evap_soi(p)
       write(unitno,'(a,es24.16)') '  fpsn_final = ', outputs%fpsn_patch(p)
    end do
    close(unitno)

    deallocate(target_names)

  end subroutine debug_dump_canopyflux_targets

  subroutine debug_print_canopyflux_psn_state(stage, p, inputs, outputs)

    character(len=*), intent(in) :: stage
    integer, intent(in) :: p
    type(canopyflux_emulator_input_view_type), intent(in) :: inputs
    type(canopyflux_emulator_output_view_type), intent(in) :: outputs

    if (.not. debug_canopyflux_emulator) return
    if (debug_canopyflux_patch > 0 .and. p /= debug_canopyflux_patch) return

    write(*,*) 'CanopyFluxesEmulator ', trim(stage), ': nstep=', get_nstep(), &
         ' patch=', p, ' psnsun=', outputs%psnsun_patch(p), &
         ' psnsha=', outputs%psnsha_patch(p), ' fpsn=', outputs%fpsn_patch(p), &
         ' sabv=', inputs%sabv_patch(p), ' elai=', inputs%elai_patch(p), &
         ' laisun=', inputs%laisun_patch(p), ' laisha=', inputs%laisha_patch(p)

  end subroutine debug_print_canopyflux_psn_state

  subroutine debug_print_canopyflux_longwave_state(stage, p, inputs, outputs)

    character(len=*), intent(in) :: stage
    integer, intent(in) :: p
    type(canopyflux_emulator_input_view_type), intent(in) :: inputs
    type(canopyflux_emulator_output_view_type), intent(in) :: outputs
    integer :: c, t

    if (.not. debug_canopyflux_emulator) return
    if (get_nstep() /= debug_canopyflux_compare_nstep) return
    if (debug_canopyflux_patch > 0 .and. p /= debug_canopyflux_patch) return
    if (.not. masterproc) return

    c = veg_pp%column(p)
    t = veg_pp%topounit(p)
    write(*,*) 'CanopyFluxes ', trim(stage), ' longwave: nstep=', get_nstep(), &
         ' patch=', p, ' ivt=', veg_pp%itype(p), ' t_veg=', outputs%t_veg(p), &
         ' dlrad=', outputs%dlrad(p), ' ulrad=', outputs%ulrad(p), &
         ' forc_lwrad=', inputs%forc_lwrad(t), ' emv=', inputs%emv(p), &
         ' emg=', inputs%emg(c), ' frac_sno=', inputs%frac_sno(c), &
         ' frac_h2osfc=', inputs%frac_h2osfc(c), ' eflx_sh_tot=', outputs%eflx_sh_veg(p) + outputs%eflx_sh_grnd(p), &
         ' eflx_sh_veg=', outputs%eflx_sh_veg(p), ' eflx_sh_grnd=', outputs%eflx_sh_grnd(p), &
         ' eflx_lh_tot=', veg_ef%eflx_lh_tot(p), ' eflx_soil_grnd=', veg_ef%eflx_soil_grnd(p), &
         ' cgrnds=', outputs%cgrnds(p), ' cgrnd=', outputs%cgrnd(p)

  end subroutine debug_print_canopyflux_longwave_state

  subroutine log_canopyflux_timing(label, elapsed, patch_count)

    character(len=*), intent(in) :: label
    real(r8), intent(in) :: elapsed
    integer, intent(in) :: patch_count
    real(r8) :: avg_time

    if (.not. masterproc) return

    if (patch_count > 0) then
       avg_time = elapsed / real(patch_count, r8)
    else
       avg_time = 0._r8
    end if

    write(*,'(a,1x,a,1x,a,i0,1x,a,es12.5,1x,a,es12.5)') &
         'CanopyFlux timing:', trim(label), 'patches=', patch_count, 'elapsed_s=', elapsed, 'elapsed_s_per_patch=', avg_time

  end subroutine log_canopyflux_timing

  subroutine log_canopyflux_timing_summary()

    real(r8) :: native_per_call
    real(r8) :: emulator_per_call
    real(r8) :: native_per_patch
    real(r8) :: emulator_per_patch
    real(r8) :: emulator_inference_fraction
    real(r8) :: emulator_apply_fraction
    real(r8) :: emulator_rebuild_fraction
    real(r8) :: emulator_mean_group_size

    if (.not. masterproc) return

    native_per_call = 0._r8
    emulator_per_call = 0._r8
    native_per_patch = 0._r8
    emulator_per_patch = 0._r8
    emulator_inference_fraction = 0._r8
    emulator_apply_fraction = 0._r8
    emulator_rebuild_fraction = 0._r8
    emulator_mean_group_size = 0._r8

    if (canopyflux_native_call_count > 0) then
       native_per_call = canopyflux_native_time_total / real(canopyflux_native_call_count, r8)
    end if
    if (canopyflux_emulator_call_count > 0) then
       emulator_per_call = canopyflux_emulator_time_total / real(canopyflux_emulator_call_count, r8)
    end if
    if (canopyflux_native_patch_count_total > 0) then
       native_per_patch = canopyflux_native_time_total / real(canopyflux_native_patch_count_total, r8)
    end if
    if (canopyflux_emulator_patch_count_total > 0) then
       emulator_per_patch = canopyflux_emulator_time_total / real(canopyflux_emulator_patch_count_total, r8)
    end if
    if (canopyflux_emulator_time_total > 0._r8) then
       emulator_inference_fraction = canopyflux_emulator_inference_time_total / canopyflux_emulator_time_total
       emulator_apply_fraction = canopyflux_emulator_apply_time_total / canopyflux_emulator_time_total
       emulator_rebuild_fraction = canopyflux_emulator_rebuild_time_total / canopyflux_emulator_time_total
    end if
    if (canopyflux_emulator_group_count_total > 0) then
       emulator_mean_group_size = real(canopyflux_emulator_group_patch_count_total, r8) / &
            real(canopyflux_emulator_group_count_total, r8)
    end if

    write(*,'(a,1x,a,i0,1x,a,es12.5,1x,a,es12.5,1x,a,es12.5)') &
         'CanopyFlux timing summary:', 'native_calls=', canopyflux_native_call_count, &
         'native_total_s=', canopyflux_native_time_total, 'native_s_per_call=', native_per_call, &
         'native_s_per_patch=', native_per_patch
    write(*,'(a,1x,a,i0,1x,a,es12.5,1x,a,es12.5,1x,a,es12.5)') &
         'CanopyFlux timing summary:', 'emulator_calls=', canopyflux_emulator_call_count, &
         'emulator_total_s=', canopyflux_emulator_time_total, 'emulator_s_per_call=', emulator_per_call, &
         'emulator_s_per_patch=', emulator_per_patch
    write(*,'(a,1x,a,es12.5,1x,a,es12.5,1x,a,es12.5)') &
         'CanopyFlux timing summary:', 'emulator_inference_s=', canopyflux_emulator_inference_time_total, &
         'emulator_apply_s=', canopyflux_emulator_apply_time_total, &
         'emulator_rebuild_s=', canopyflux_emulator_rebuild_time_total
    write(*,'(a,1x,a,f8.4,1x,a,f8.4,1x,a,f8.4)') &
         'CanopyFlux timing summary:', 'emulator_inference_frac=', emulator_inference_fraction, &
         'emulator_apply_frac=', emulator_apply_fraction, &
         'emulator_rebuild_frac=', emulator_rebuild_fraction
    write(*,'(a,1x,a,i0,1x,a,f8.3,1x,a,i0)') &
         'CanopyFlux timing summary:', 'emulator_groups=', canopyflux_emulator_group_count_total, &
         'emulator_mean_group_size=', emulator_mean_group_size, &
         'emulator_max_group_size=', canopyflux_emulator_max_group_size

  end subroutine log_canopyflux_timing_summary

  subroutine debug_print_canopyflux_late_state(filter_nolakeurbanp, inputs, outputs)

    integer, intent(in) :: filter_nolakeurbanp(:)
    type(canopyflux_emulator_input_view_type), intent(in) :: inputs
    type(canopyflux_emulator_output_view_type), intent(in) :: outputs

    integer :: fp, p

    do fp = 1, size(filter_nolakeurbanp)
       p = filter_nolakeurbanp(fp)
       call debug_print_canopyflux_psn_state('late', p, inputs, outputs)
    end do

  end subroutine debug_print_canopyflux_late_state

  subroutine cache_emulator_photosyn_outputs(filter_nolakeurbanp, outputs, psnsun_saved, psnsha_saved, fpsn_saved)

    integer, intent(in) :: filter_nolakeurbanp(:)
    type(canopyflux_emulator_output_view_type), intent(in) :: outputs
    real(r8), intent(out) :: psnsun_saved(:)
    real(r8), intent(out) :: psnsha_saved(:)
    real(r8), intent(out) :: fpsn_saved(:)

    integer :: fp, p

    do fp = 1, size(filter_nolakeurbanp)
       p = filter_nolakeurbanp(fp)
       psnsun_saved(fp) = outputs%psnsun_patch(p)
       psnsha_saved(fp) = outputs%psnsha_patch(p)
       fpsn_saved(fp) = outputs%fpsn_patch(p)
    end do

  end subroutine cache_emulator_photosyn_outputs

  subroutine restore_emulator_photosyn_outputs(filter_nolakeurbanp, outputs, psnsun_saved, psnsha_saved, fpsn_saved)

    integer, intent(in) :: filter_nolakeurbanp(:)
    type(canopyflux_emulator_output_view_type), intent(inout) :: outputs
    real(r8), intent(in) :: psnsun_saved(:)
    real(r8), intent(in) :: psnsha_saved(:)
    real(r8), intent(in) :: fpsn_saved(:)

    integer :: fp, p

    do fp = 1, size(filter_nolakeurbanp)
       p = filter_nolakeurbanp(fp)
       outputs%psnsun_patch(p) = psnsun_saved(fp)
       outputs%psnsha_patch(p) = psnsha_saved(fp)
       outputs%fpsn_patch(p) = fpsn_saved(fp)
    end do

  end subroutine restore_emulator_photosyn_outputs

  subroutine predict_canopyfluxes(bounds, num_nolakeurbanp, filter_nolakeurbanp, inputs, outputs)

    type(bounds_type)                            , intent(in)    :: bounds
    integer                                      , intent(in)    :: num_nolakeurbanp
    integer                                      , intent(in)    :: filter_nolakeurbanp(:)
    type(canopyflux_emulator_input_view_type)    , intent(in)    :: inputs
    type(canopyflux_emulator_output_view_type)   , intent(inout) :: outputs
    real(r8) :: features(num_nolakeurbanp, canopyflux_emulator_num_features)
    real(r8), allocatable :: predictions(:,:)
    real(r8), allocatable :: group_features(:,:)
    real(r8), allocatable :: batched_group_features(:,:)
    real(r8), allocatable :: batched_predictions(:,:)
    integer, allocatable :: group_filter(:)
    logical :: grouped(num_nolakeurbanp)
    integer :: fp, other_fp, p, c, ivt
    integer :: group_count, group_fp, other_ivt
    integer :: batch_count, repeat_idx, sample_offset
    real(r8) :: t_start, t_end, elapsed
    real(r8) :: inference_elapsed, apply_elapsed

    call assemble_canopyflux_emulator_features(bounds, num_nolakeurbanp, filter_nolakeurbanp, &
         inputs, features)
    call debug_dump_canopyflux_features(filter_nolakeurbanp, features)
    call log_required_canopyflux_pfts(filter_nolakeurbanp)

    grouped(:) = .false.
    do fp = 1, num_nolakeurbanp
       p = filter_nolakeurbanp(fp)
       if (p < 1 .or. p > size(veg_pp%itype)) cycle
       c = veg_pp%column(p)
    end do

    do fp = 1, num_nolakeurbanp
       if (grouped(fp)) cycle
       p = filter_nolakeurbanp(fp)
       if (p < 1 .or. p > size(veg_pp%itype)) then
          grouped(fp) = .true.
          cycle
       end if
       ivt = veg_pp%itype(p)
       call ensure_canopyflux_model_loaded(ivt)

       group_count = 0
       do other_fp = fp, num_nolakeurbanp
          if (grouped(other_fp)) cycle
          p = filter_nolakeurbanp(other_fp)
          if (p < 1 .or. p > size(veg_pp%itype)) cycle
          other_ivt = veg_pp%itype(p)
          if (other_ivt == ivt) then
             group_count = group_count + 1
          end if
       end do

       canopyflux_emulator_group_count_total = canopyflux_emulator_group_count_total + 1
       canopyflux_emulator_group_patch_count_total = canopyflux_emulator_group_patch_count_total + group_count
       canopyflux_emulator_max_group_size = max(canopyflux_emulator_max_group_size, group_count)

       allocate(group_filter(group_count))
       allocate(group_features(group_count, canopyflux_emulator_num_features))
       allocate(predictions(group_count, size(canopyflux_models(ivt)%target_map)))
       batch_count = group_count * canopyflux_timing_repeat_count
       allocate(batched_group_features(batch_count, canopyflux_emulator_num_features))
       allocate(batched_predictions(batch_count, size(canopyflux_models(ivt)%target_map)))

       group_fp = 0
       do other_fp = fp, num_nolakeurbanp
          if (grouped(other_fp)) cycle
          p = filter_nolakeurbanp(other_fp)
          if (p < 1 .or. p > size(veg_pp%itype)) cycle
          other_ivt = veg_pp%itype(p)
          if (other_ivt == ivt) then
             group_fp = group_fp + 1
             group_filter(group_fp) = p
             group_features(group_fp,:) = features(other_fp,:)
             grouped(other_fp) = .true.
          end if
       end do

       do repeat_idx = 1, canopyflux_timing_repeat_count
          sample_offset = (repeat_idx - 1) * group_count
          batched_group_features(sample_offset + 1:sample_offset + group_count,:) = group_features(:,:)
       end do
       call cpu_time(t_start)
       call run_canopyflux_model(canopyflux_models(ivt), batched_group_features, batched_predictions)
       call cpu_time(t_end)
       inference_elapsed = t_end - t_start
       canopyflux_emulator_inference_time_total = canopyflux_emulator_inference_time_total + inference_elapsed
       predictions(:,:) = batched_predictions(1:group_count,:)
       call debug_print_canopyflux_predictions(group_filter, canopyflux_models(ivt), predictions)
       call cpu_time(t_start)
       call apply_canopyflux_predictions(group_filter, inputs, outputs, canopyflux_models(ivt), predictions)
       call cpu_time(t_end)
       apply_elapsed = t_end - t_start
       elapsed = inference_elapsed + apply_elapsed
       canopyflux_emulator_apply_time_total = canopyflux_emulator_apply_time_total + apply_elapsed
       canopyflux_emulator_time_total = canopyflux_emulator_time_total + elapsed
       canopyflux_emulator_call_count = canopyflux_emulator_call_count + 1
       canopyflux_emulator_patch_count_total = canopyflux_emulator_patch_count_total + batch_count
       call log_canopyflux_timing('emulator', elapsed, batch_count)
       call debug_dump_canopyflux_targets(group_filter, inputs, outputs, canopyflux_models(ivt), predictions)
       deallocate(group_filter, group_features, batched_group_features, batched_predictions, predictions)
    end do

  end subroutine predict_canopyfluxes

  subroutine CanopyFluxesEmulator(bounds, num_nolakeurbanp, filter_nolakeurbanp, &
       atm2lnd_vars, canopystate_vars, cnstate_vars, energyflux_vars, &
       frictionvel_vars, soilstate_vars, solarabs_vars, surfalb_vars, &
       ch4_vars, photosyns_vars, soilhydrology_vars)

    type(bounds_type)         , intent(in)    :: bounds
    integer                   , intent(in)    :: num_nolakeurbanp
    integer                   , intent(in)    :: filter_nolakeurbanp(:)
    type(atm2lnd_type)        , intent(inout) :: atm2lnd_vars
    type(canopystate_type)    , intent(inout) :: canopystate_vars
    type(cnstate_type)        , intent(inout) :: cnstate_vars
    type(energyflux_type)     , intent(inout) :: energyflux_vars
    type(frictionvel_type)    , intent(inout) :: frictionvel_vars
    type(soilstate_type)      , intent(inout) :: soilstate_vars
    type(solarabs_type)       , intent(inout) :: solarabs_vars
    type(surfalb_type)        , intent(inout) :: surfalb_vars
    type(ch4_type)            , intent(inout) :: ch4_vars
    type(photosyns_type)      , intent(inout) :: photosyns_vars
    type(soilhydrology_type)  , intent(in)    :: soilhydrology_vars
    type(canopyflux_emulator_input_view_type) :: inputs
    real(r8) :: features(num_nolakeurbanp, canopyflux_emulator_num_features)
    real(r8), allocatable :: training_flnr(:)
    real(r8), allocatable :: training_mbbopt(:)
   real(r8), allocatable :: training_pco2(:)
    real(r8), allocatable :: flnr_saved(:)
    real(r8), allocatable :: mbbopt_saved(:)
   real(r8), allocatable :: forc_pco2_saved(:)
    real(r8) :: training_flnr_value, training_mbbopt_value
   real(r8) :: training_pco2_value
   logical :: training_capture_active

    call inputs%bind(atm2lnd_vars, canopystate_vars, energyflux_vars, frictionvel_vars, &
         soilstate_vars, solarabs_vars, surfalb_vars, soilhydrology_vars)
    call set_perchroot_opt(perchroot, perchroot_alt)
    ! The native canopy-flux path zeros btran/btran2 before accumulating root stress.
    ! The emulator capture path needs to do the same, otherwise the training debug
    ! capture can inherit spval from the EnergyFluxType allocation state.
    energyflux_vars%btran_patch(:) = 0._r8
    energyflux_vars%btran2_patch(:) = 0._r8
    call calc_root_moist_stress(bounds,     &
         nlevgrnd = nlevgrnd,               &
         fn = num_nolakeurbanp,             &
         filterp = filter_nolakeurbanp,     &
         canopystate_vars = canopystate_vars, &
         energyflux_vars = energyflux_vars, &
         soilstate_vars = soilstate_vars)
    training_capture_active = canopyflux_training_capture_enabled()
    if (training_capture_active) then
       allocate(training_pco2(num_nolakeurbanp))
       call sample_canopyflux_training_pco2(training_pco2_value)
       training_pco2(:) = training_pco2_value
       if (randomize_canopyflux_training_traits) then
          allocate(training_flnr(num_nolakeurbanp))
          allocate(training_mbbopt(num_nolakeurbanp))
          call sample_canopyflux_training_traits(training_flnr_value, training_mbbopt_value)
          training_flnr(:) = training_flnr_value
          training_mbbopt(:) = training_mbbopt_value
          call assemble_canopyflux_emulator_features(bounds, num_nolakeurbanp, filter_nolakeurbanp, &
               inputs, features, flnr_override=training_flnr, mbbopt_override=training_mbbopt, pco2_override=training_pco2)
       else
          call assemble_canopyflux_emulator_features(bounds, num_nolakeurbanp, filter_nolakeurbanp, &
               inputs, features, pco2_override=training_pco2)
       end if
       call capture_canopyflux_training_features(num_nolakeurbanp, filter_nolakeurbanp, features, energyflux_vars%btran_patch)
       allocate(forc_pco2_saved(lbound(inputs%forc_pco2, 1):ubound(inputs%forc_pco2, 1)))
       forc_pco2_saved(:) = inputs%forc_pco2(:)
       inputs%forc_pco2(:) = training_pco2_value
       if (randomize_canopyflux_training_traits) then
          allocate(flnr_saved(lbound(veg_vp%flnr, 1):ubound(veg_vp%flnr, 1)))
          allocate(mbbopt_saved(lbound(veg_vp%mbbopt, 1):ubound(veg_vp%mbbopt, 1)))
          flnr_saved(:) = veg_vp%flnr(:)
          mbbopt_saved(:) = veg_vp%mbbopt(:)
          veg_vp%flnr(:) = training_flnr_value
          veg_vp%mbbopt(:) = training_mbbopt_value
       end if
    end if

    call CanopyFluxes(bounds, num_nolakeurbanp, filter_nolakeurbanp, &
         atm2lnd_vars, canopystate_vars, cnstate_vars, energyflux_vars, &
         frictionvel_vars, soilstate_vars, solarabs_vars, surfalb_vars, &
         ch4_vars, photosyns_vars)

    if (training_capture_active) then
       inputs%forc_pco2(:) = forc_pco2_saved(:)
       deallocate(forc_pco2_saved, training_pco2)
       if (randomize_canopyflux_training_traits) then
          veg_vp%flnr(:) = flnr_saved(:)
          veg_vp%mbbopt(:) = mbbopt_saved(:)
          deallocate(flnr_saved, mbbopt_saved, training_flnr, training_mbbopt)
       end if
    end if

  end subroutine CanopyFluxesEmulator

end module CanopyFluxesEmulatorMod
