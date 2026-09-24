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
  real(r8), public :: peat_compaction_surface_density = 25._r8
  real(r8), public :: peat_compaction_deep_density = 100._r8
  real(r8), public :: peat_compaction_efolding_depth = 0.25_r8
  real(r8), public :: peat_compaction_timescale_years = 1._r8
  real(r8), public :: soil_ice_impedance_exponent = 6._r8
  !$acc declare create(anoxia_wtsat)
  !$acc declare create(nlev_soildecomp_standard)
  !$acc declare copyin(moss_capillary_max_demand_fraction)
  !$acc declare copyin(moss_capillary_connectivity_timescale_days)
  !$acc declare copyin(peat_compaction_surface_density)
  !$acc declare copyin(peat_compaction_deep_density)
  !$acc declare copyin(peat_compaction_efolding_depth)
  !$acc declare copyin(peat_compaction_timescale_years)
  !$acc declare copyin(soil_ice_impedance_exponent)

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

        tString='peat_compaction_surface_density'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_compaction_surface_density=tempr

        tString='peat_compaction_deep_density'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_compaction_deep_density=tempr

        tString='peat_compaction_efolding_depth'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_compaction_efolding_depth=tempr

        tString='peat_compaction_timescale_years'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) peat_compaction_timescale_years=tempr

        tString='soil_ice_impedance_exponent'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if (readv) soil_ice_impedance_exponent=tempr

        if (moss_capillary_max_demand_fraction < 0._r8 .or. &
             moss_capillary_max_demand_fraction > 1._r8) then
           call endrun(msg=' ERROR: moss_capillary_max_demand_fraction must be in [0,1]'//&
                errMsg(__FILE__, __LINE__))
        end if
        if (moss_capillary_connectivity_timescale_days <= 0._r8) then
           call endrun(msg=' ERROR: moss_capillary_connectivity_timescale_days must be positive'//&
                errMsg(__FILE__, __LINE__))
        end if
        if (soil_ice_impedance_exponent < 0._r8) then
           call endrun(msg=' ERROR: soil_ice_impedance_exponent must be nonnegative'//&
                errMsg(__FILE__, __LINE__))
        end if
        if (peat_compaction_surface_density <= 0._r8 .or. &
             peat_compaction_deep_density < peat_compaction_surface_density .or. &
             peat_compaction_efolding_depth <= 0._r8 .or. &
             peat_compaction_timescale_years <= 0._r8) then
           call endrun(msg=' ERROR: peat compaction densities must be positive, deep density '//&
                'must be at least surface density, and depth/time scales must be positive'//&
                errMsg(__FILE__, __LINE__))
        end if

        !$acc update device(moss_capillary_max_demand_fraction)
        !$acc update device(moss_capillary_connectivity_timescale_days)
        !$acc update device(peat_compaction_surface_density)
        !$acc update device(peat_compaction_deep_density)
        !$acc update device(peat_compaction_efolding_depth)
        !$acc update device(peat_compaction_timescale_years)
        !$acc update device(soil_ice_impedance_exponent)

        if (masterproc) then
           write(iulog,*) 'Peatland parameters read from ELM parameter file:'
           write(iulog,*) '  moss_capillary_max_demand_fraction = ', moss_capillary_max_demand_fraction
           write(iulog,*) '  moss_capillary_connectivity_timescale_days = ', &
                moss_capillary_connectivity_timescale_days
           write(iulog,*) '  peat_compaction_surface_density (kg C/m3) = ', peat_compaction_surface_density
           write(iulog,*) '  peat_compaction_deep_density (kg C/m3) = ', peat_compaction_deep_density
           write(iulog,*) '  peat_compaction_efolding_depth (m) = ', peat_compaction_efolding_depth
           write(iulog,*) '  peat_compaction_timescale_years (yr) = ', peat_compaction_timescale_years
           write(iulog,*) '  soil_ice_impedance_exponent = ', soil_ice_impedance_exponent
        end if
     end if

   end subroutine ParamsReadShared

end module SharedParamsMod
