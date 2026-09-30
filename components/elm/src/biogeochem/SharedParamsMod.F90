module SharedParamsMod

  !-----------------------------------------------------------------------
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  implicit none
  save

  ! ParamsShareInst.  PGI wants the type decl. public but the instance
  ! is indeed protected.  A generic private statement at the start of the module
  ! overrides the protected functionality with PGI

  type, public :: ParamsShareType

      real(r8), pointer :: Q10_mr                => null() ! temperature dependence for maintenance respiraton
      real(r8), pointer :: Q10_hr                => null() ! temperature dependence for heterotrophic respiration
      real(r8), pointer :: minpsi                => null() ! minimum soil water potential for heterotrophic resp
      real(r8), pointer :: cwd_fcel              => null() ! cellulose fraction of coarse woody debris
      real(r8), pointer :: cwd_flig              => null() ! lignin fraction of coarse woody debris
      real(r8), pointer :: froz_q10              => null() ! separate q10 for frozen soil respiration rates
      real(r8), pointer :: decomp_depth_efolding => null() ! e-folding depth for reduction in decomposition (m)
      real(r8), pointer :: mino2lim              => null() ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate
      real(r8), pointer :: organic_max           => null() ! organic matter content (kg/m3) where soil is assumed to act like peat

  end type ParamsShareType

  type(ParamsShareType),protected :: ParamsShareInst

  !$acc declare create(ParamsShareInst)
  logical, public :: anoxia_wtsat = .false.
  integer, public :: nlev_soildecomp_standard = 5
  real(r8), public :: moss_capillary_max_demand_fraction = 0.10_r8
  real(r8), public :: moss_capillary_connectivity_timescale_days = 30._r8
  real(r8), public :: moss_water_layer_thickness = 0.05_r8
  real(r8), public :: moss_water_saturated_suction = 10._r8
  real(r8), public :: moss_water_clapp_hornberger_b = 3.5_r8
  real(r8), public :: moss_water_hydraulic_conductivity_sat = 1.5_r8
  real(r8), public :: moss_water_conductivity_exponent = 5._r8
  real(r8), public :: moss_initial_carbon = 10._r8
  real(r8), public :: moss_carbon_fraction_dry_mass = 0.5_r8
  real(r8), public :: moss_water_initial_content = 5._r8
  real(r8), public :: moss_water_content_min = 0.5_r8
  real(r8), public :: moss_water_content_max = 8._r8
  real(r8), public :: moss_water_drainage_threshold = 8.5_r8
  real(r8), public :: moss_water_substrate_rewetting_timescale_days = 1._r8
  real(r8), public :: peat_compaction_surface_density = 20._r8
  real(r8), public :: peat_compaction_deep_density = 90._r8
  real(r8), public :: peat_compaction_initial_acrotelm_depth = 0.25_r8
  real(r8), public :: peat_compaction_transition_width = 0.10_r8
  real(r8), public :: peat_compaction_zwt_running_mean_days = 30._r8
  real(r8), public :: peat_compaction_zwt_smoothing_years = 10._r8
  real(r8), public :: peat_compaction_growing_season_start_doy = 121._r8
  real(r8), public :: peat_compaction_growing_season_end_doy = 305._r8
  real(r8), public :: peat_compaction_timescale_years = 1._r8
  real(r8), public :: soil_ice_impedance_exponent = 6._r8
  real(r8), public :: peat_hydraulic_b_max = 12._r8
  real(r8), public :: peat_hydraulic_hksat_min = 1.e-4_r8
  real(r8), public :: peat_hydraulic_carbon_fraction = 0.56_r8
  real(r8), public :: peat_hydraulic_particle_density = 1260._r8
  real(r8), public :: peat_hydraulic_bulk_density_min = 35._r8
  real(r8), public :: peat_hydraulic_bulk_density_max = 210._r8
  !$acc declare create(anoxia_wtsat)
  !$acc declare create(nlev_soildecomp_standard)
  !$acc declare copyin(moss_capillary_max_demand_fraction)
  !$acc declare copyin(moss_capillary_connectivity_timescale_days)
  !$acc declare copyin(moss_water_layer_thickness)
  !$acc declare copyin(moss_water_saturated_suction)
  !$acc declare copyin(moss_water_clapp_hornberger_b)
  !$acc declare copyin(moss_water_hydraulic_conductivity_sat)
  !$acc declare copyin(moss_water_conductivity_exponent)
  !$acc declare copyin(moss_initial_carbon)
  !$acc declare copyin(moss_carbon_fraction_dry_mass)
  !$acc declare copyin(moss_water_initial_content)
  !$acc declare copyin(moss_water_content_min)
  !$acc declare copyin(moss_water_content_max)
  !$acc declare copyin(moss_water_drainage_threshold)
  !$acc declare copyin(moss_water_substrate_rewetting_timescale_days)
  !$acc declare copyin(peat_compaction_surface_density)
  !$acc declare copyin(peat_compaction_deep_density)
  !$acc declare copyin(peat_compaction_initial_acrotelm_depth)
  !$acc declare copyin(peat_compaction_transition_width)
  !$acc declare copyin(peat_compaction_zwt_running_mean_days)
  !$acc declare copyin(peat_compaction_zwt_smoothing_years)
  !$acc declare copyin(peat_compaction_growing_season_start_doy)
  !$acc declare copyin(peat_compaction_growing_season_end_doy)
  !$acc declare copyin(peat_compaction_timescale_years)
  !$acc declare copyin(soil_ice_impedance_exponent)
  !$acc declare copyin(peat_hydraulic_b_max)
  !$acc declare copyin(peat_hydraulic_hksat_min)
  !$acc declare copyin(peat_hydraulic_carbon_fraction)
  !$acc declare copyin(peat_hydraulic_particle_density)
  !$acc declare copyin(peat_hydraulic_bulk_density_min)
  !$acc declare copyin(peat_hydraulic_bulk_density_max)

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
   subroutine ParamsReadShared(ncid)
     !
     use ncdio_pio   , only : file_desc_t,ncd_io
     use abortutils  , only : endrun
     use shr_log_mod , only : errMsg => shr_log_errMsg
     use elm_varctl  , only : use_humhol, iulog
     use spmdMod     , only : masterproc
     !
     implicit none
     type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
     !
     character(len=32)  :: subname = 'ParamsReadShared'
     character(len=100) :: errCode = '-Error reading in CN and BGC shared params file. Var:'
     logical            :: readv ! has variable been read in or not
     real(r8)           :: tempr ! temporary to read in parameter
     character(len=100) :: tString ! temp. var for reading
     !-----------------------------------------------------------------------
     !
     ! netcdf read here
     !
     allocate(ParamsShareInst%Q10_mr               )
     allocate(ParamsShareInst%Q10_hr               )
     allocate(ParamsShareInst%minpsi               )
     allocate(ParamsShareInst%cwd_fcel             )
     allocate(ParamsShareInst%cwd_flig             )
     allocate(ParamsShareInst%froz_q10             )
     allocate(ParamsShareInst%decomp_depth_efolding)
     allocate(ParamsShareInst%mino2lim             )
     allocate(ParamsShareInst%organic_max          )
     tString='q10_mr'
     call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%Q10_mr=tempr

     tString='q10_hr'
     call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%Q10_hr=tempr


     tString='minpsi_hr'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%minpsi=tempr

     tString='cwd_fcel'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%cwd_fcel=tempr

     tString='cwd_flig'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%cwd_flig=tempr

     tString='froz_q10'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%froz_q10=tempr

     tString='decomp_depth_efolding'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%decomp_depth_efolding=tempr

     tString='mino2lim'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%mino2lim=tempr
     !ParamsShareInst%mino2lim=0.2_r8

     tString='organic_max'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%organic_max=tempr

     ! These parameters only affect runtime-gated peatland processes.  Read
     ! them from the standard ELM parameter file for HUMHOL cases, while
     ! retaining the established non-peatland values for BFB compatibility.
     if (use_humhol) then
        ! Preserve the previous generated HUMHOL default for old parameter
        ! files that predate this migration.
        soil_ice_impedance_exponent = 8._r8

        tString='moss_capillary_max_demand_fraction'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_capillary_max_demand_fraction=tempr

        tString='moss_capillary_connectivity_timescale_days'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_capillary_connectivity_timescale_days=tempr

        tString='moss_water_layer_thickness'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_water_layer_thickness=tempr

        tString='moss_water_saturated_suction'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_water_saturated_suction=tempr

        tString='moss_water_clapp_hornberger_b'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_water_clapp_hornberger_b=tempr

        tString='moss_water_hydraulic_conductivity_sat'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_water_hydraulic_conductivity_sat=tempr

        tString='moss_water_conductivity_exponent'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_water_conductivity_exponent=tempr

        tString='moss_initial_carbon'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_initial_carbon=tempr

        tString='moss_carbon_fraction_dry_mass'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_carbon_fraction_dry_mass=tempr

        tString='moss_water_initial_content'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_water_initial_content=tempr

        tString='moss_water_content_min'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_water_content_min=tempr

        tString='moss_water_content_max'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_water_content_max=tempr

        tString='moss_water_drainage_threshold'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_water_drainage_threshold=tempr

        tString='moss_water_substrate_rewetting_timescale_days'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) moss_water_substrate_rewetting_timescale_days=tempr

        tString='peat_compaction_surface_density'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_compaction_surface_density=tempr

        tString='peat_compaction_deep_density'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_compaction_deep_density=tempr

        tString='peat_compaction_initial_acrotelm_depth'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) then
           peat_compaction_initial_acrotelm_depth=tempr
        else
           ! Backward-compatible interpretation of the former fixed-profile
           ! e-folding depth while parameter files transition to the
           ! hydrologically diagnosed acrotelm boundary.
           tString='peat_compaction_efolding_depth'
           call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
           if (readv) peat_compaction_initial_acrotelm_depth=tempr
        end if

        tString='peat_compaction_transition_width'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_compaction_transition_width=tempr

        tString='peat_compaction_zwt_running_mean_days'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_compaction_zwt_running_mean_days=tempr

        tString='peat_compaction_zwt_smoothing_years'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_compaction_zwt_smoothing_years=tempr

        tString='peat_compaction_growing_season_start_doy'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_compaction_growing_season_start_doy=tempr

        tString='peat_compaction_growing_season_end_doy'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_compaction_growing_season_end_doy=tempr

        tString='peat_compaction_timescale_years'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_compaction_timescale_years=tempr

        tString='soil_ice_impedance_exponent'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) soil_ice_impedance_exponent=tempr

        tString='peat_hydraulic_b_max'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_hydraulic_b_max=tempr

        tString='peat_hydraulic_hksat_min'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_hydraulic_hksat_min=tempr

        tString='peat_hydraulic_carbon_fraction'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_hydraulic_carbon_fraction=tempr

        tString='peat_hydraulic_particle_density'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_hydraulic_particle_density=tempr

        tString='peat_hydraulic_bulk_density_min'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_hydraulic_bulk_density_min=tempr

        tString='peat_hydraulic_bulk_density_max'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_hydraulic_bulk_density_max=tempr

        if (moss_capillary_max_demand_fraction < 0._r8 .or. &
             moss_capillary_max_demand_fraction > 1._r8) then
           call endrun(msg=' ERROR: moss_capillary_max_demand_fraction must be in [0,1]'//&
                errMsg(__FILE__, __LINE__))
        end if
        if (moss_capillary_connectivity_timescale_days <= 0._r8) then
           call endrun(msg=' ERROR: moss_capillary_connectivity_timescale_days must be positive'//&
                errMsg(__FILE__, __LINE__))
        end if
        if (moss_water_layer_thickness <= 0._r8 .or. &
             moss_water_saturated_suction <= 0._r8 .or. &
             moss_water_clapp_hornberger_b <= 0._r8 .or. &
             moss_water_hydraulic_conductivity_sat <= 0._r8 .or. &
             moss_water_conductivity_exponent <= 0._r8 .or. &
             moss_initial_carbon <= 0._r8 .or. &
             moss_carbon_fraction_dry_mass <= 0._r8 .or. &
             moss_carbon_fraction_dry_mass > 1._r8 .or. &
             moss_water_initial_content < 0._r8 .or. &
             moss_water_content_min < 0._r8 .or. &
             moss_water_content_max <= moss_water_content_min .or. &
             moss_water_drainage_threshold < moss_water_content_max .or. &
             moss_water_substrate_rewetting_timescale_days <= 0._r8 .or. &
             moss_water_initial_content > moss_water_content_max) then
           call endrun(msg=' ERROR: invalid prognostic moss-water parameters'//&
                errMsg(__FILE__, __LINE__))
        end if
        if (soil_ice_impedance_exponent < 0._r8) then
           call endrun(msg=' ERROR: soil_ice_impedance_exponent must be nonnegative'//&
                errMsg(__FILE__, __LINE__))
        end if
        if (peat_hydraulic_b_max < 2.7_r8) then
           call endrun(msg=' ERROR: peat_hydraulic_b_max must be at least 2.7'//&
                errMsg(__FILE__, __LINE__))
        end if
        if (peat_hydraulic_hksat_min <= 0._r8 .or. &
             peat_hydraulic_hksat_min > 0.28_r8) then
           call endrun(msg=' ERROR: peat_hydraulic_hksat_min must be in (0,0.28] mm/s'//&
                errMsg(__FILE__, __LINE__))
        end if
        if (peat_hydraulic_carbon_fraction <= 0._r8 .or. &
             peat_hydraulic_carbon_fraction > 1._r8 .or. &
             peat_hydraulic_particle_density <= 0._r8 .or. &
             peat_hydraulic_bulk_density_min <= 0._r8 .or. &
             peat_hydraulic_bulk_density_max <= peat_hydraulic_bulk_density_min .or. &
             peat_hydraulic_bulk_density_max >= peat_hydraulic_particle_density) then
           call endrun(msg=' ERROR: invalid JULES-Peat density conversion parameters'//&
                errMsg(__FILE__, __LINE__))
        end if
        if (peat_compaction_surface_density <= 0._r8 .or. &
             peat_compaction_deep_density < peat_compaction_surface_density .or. &
             peat_compaction_initial_acrotelm_depth <= 0._r8 .or. &
             peat_compaction_transition_width <= 0._r8 .or. &
             peat_compaction_zwt_running_mean_days <= 0._r8 .or. &
             peat_compaction_zwt_smoothing_years <= 0._r8 .or. &
             peat_compaction_growing_season_start_doy < 1._r8 .or. &
             peat_compaction_growing_season_end_doy > 366._r8 .or. &
             peat_compaction_growing_season_end_doy < &
                  peat_compaction_growing_season_start_doy .or. &
             peat_compaction_timescale_years <= 0._r8) then
           call endrun(msg=' ERROR: peat compaction densities must be positive, deep density '//&
                'must be at least surface density, and hydrologic depth/time settings must be valid'//&
                errMsg(__FILE__, __LINE__))
        end if

        !$acc update device(moss_capillary_max_demand_fraction)
        !$acc update device(moss_capillary_connectivity_timescale_days)
        !$acc update device(moss_water_layer_thickness)
        !$acc update device(moss_water_saturated_suction)
        !$acc update device(moss_water_clapp_hornberger_b)
        !$acc update device(moss_water_hydraulic_conductivity_sat)
        !$acc update device(moss_water_conductivity_exponent)
        !$acc update device(moss_initial_carbon)
        !$acc update device(moss_carbon_fraction_dry_mass)
        !$acc update device(moss_water_initial_content)
        !$acc update device(moss_water_content_min)
        !$acc update device(moss_water_content_max)
        !$acc update device(moss_water_drainage_threshold)
        !$acc update device(moss_water_substrate_rewetting_timescale_days)
        !$acc update device(peat_compaction_surface_density)
        !$acc update device(peat_compaction_deep_density)
        !$acc update device(peat_compaction_initial_acrotelm_depth)
        !$acc update device(peat_compaction_transition_width)
        !$acc update device(peat_compaction_zwt_running_mean_days)
        !$acc update device(peat_compaction_zwt_smoothing_years)
        !$acc update device(peat_compaction_growing_season_start_doy)
        !$acc update device(peat_compaction_growing_season_end_doy)
        !$acc update device(peat_compaction_timescale_years)
        !$acc update device(soil_ice_impedance_exponent)
        !$acc update device(peat_hydraulic_b_max)
        !$acc update device(peat_hydraulic_hksat_min)
        !$acc update device(peat_hydraulic_carbon_fraction)
        !$acc update device(peat_hydraulic_particle_density)
        !$acc update device(peat_hydraulic_bulk_density_min)
        !$acc update device(peat_hydraulic_bulk_density_max)

        if (masterproc) then
           write(iulog,*) 'Peatland parameters read from ELM parameter file:'
           write(iulog,*) '  moss_capillary_max_demand_fraction = ', moss_capillary_max_demand_fraction
           write(iulog,*) '  moss_capillary_connectivity_timescale_days = ', &
                moss_capillary_connectivity_timescale_days
           write(iulog,*) '  moss_water_layer_thickness (m) = ', moss_water_layer_thickness
           write(iulog,*) '  moss_water saturated suction (mm) = ', &
                moss_water_saturated_suction
           write(iulog,*) '  moss_water Clapp-Hornberger b = ', &
                moss_water_clapp_hornberger_b
           write(iulog,*) '  moss_water Ksat (mm/s) and exponent = ', &
                moss_water_hydraulic_conductivity_sat, moss_water_conductivity_exponent
           write(iulog,*) '  moss initial carbon (gC/m2) and C fraction of dry mass = ', &
                moss_initial_carbon, moss_carbon_fraction_dry_mass
           write(iulog,*) '  moss initial water content (gH2O/g dry mass) = ', &
                moss_water_initial_content
           write(iulog,*) '  moss physiological water-content range (gH2O/g dry mass) = ', &
                moss_water_content_min, moss_water_content_max
           write(iulog,*) '  moss drainage threshold (gH2O/g dry mass) = ', &
                moss_water_drainage_threshold
           write(iulog,*) '  moss substrate-rewetting timescale (days) = ', &
                moss_water_substrate_rewetting_timescale_days
           write(iulog,*) '  peat_compaction_surface_density (kg C/m3) = ', peat_compaction_surface_density
           write(iulog,*) '  peat_compaction_deep_density (kg C/m3) = ', peat_compaction_deep_density
           write(iulog,*) '  peat_compaction_initial_acrotelm_depth (m) = ', &
                peat_compaction_initial_acrotelm_depth
           write(iulog,*) '  peat_compaction_transition_width (m) = ', peat_compaction_transition_width
           write(iulog,*) '  peat_compaction_zwt_running_mean_days = ', &
                peat_compaction_zwt_running_mean_days
           write(iulog,*) '  peat_compaction_zwt_smoothing_years = ', &
                peat_compaction_zwt_smoothing_years
           write(iulog,*) '  peat_compaction growing-season DOY range = ', &
                peat_compaction_growing_season_start_doy, &
                peat_compaction_growing_season_end_doy
           write(iulog,*) '  peat_compaction_timescale_years (yr) = ', peat_compaction_timescale_years
           write(iulog,*) '  soil_ice_impedance_exponent = ', soil_ice_impedance_exponent
           write(iulog,*) '  peat_hydraulic_b_max = ', peat_hydraulic_b_max
           write(iulog,*) '  peat_hydraulic_hksat_min (mm/s) = ', peat_hydraulic_hksat_min
           write(iulog,*) '  peat_hydraulic_carbon_fraction = ', peat_hydraulic_carbon_fraction
           write(iulog,*) '  peat_hydraulic_particle_density (kg/m3) = ', &
                peat_hydraulic_particle_density
           write(iulog,*) '  peat_hydraulic_bulk_density range (kg/m3) = ', &
                peat_hydraulic_bulk_density_min, peat_hydraulic_bulk_density_max
        end if
     end if

   end subroutine ParamsReadShared

end module SharedParamsMod
