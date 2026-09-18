module MicrobeMethaneParamsMod

  !-----------------------------------------------------------------------
  ! Named parameter-file inputs for the revised microbial methane backend.
  ! Values are read only when use_microbe_methane is enabled.  No defaults
  ! are embedded here: an enabled case must use a parameter file containing
  ! the complete, versioned parameter set.
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use abortutils, only : endrun
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

  implicit none
  private
  save

  type, public :: MicrobeMethaneParamsType
     real(r8) :: mfg_biomass_min
     real(r8) :: k_acetate
     real(r8) :: acetate_prod_max
     real(r8) :: k_acetate_prod_o2
     real(r8) :: acetogenesis_max
     real(r8) :: k_acetogenesis_h2
     real(r8) :: k_acetogenesis_co2
     real(r8) :: acetogenesis_q10
     real(r8) :: h2_methanogen_growth_rate
     real(r8) :: h2_methanogen_death_rate
     real(r8) :: h2_methanogen_yield
     real(r8) :: acetate_methanogen_growth_rate
     real(r8) :: acetate_methanogen_death_rate
     real(r8) :: acetate_methanogen_yield
     real(r8) :: aerobic_methanotroph_growth_rate
     real(r8) :: aerobic_methanotroph_death_rate
     real(r8) :: aerobic_methanotroph_yield
     real(r8) :: anaerobic_methanotroph_growth_rate
     real(r8) :: anaerobic_methanotroph_death_rate
     real(r8) :: anaerobic_methanotroph_yield
     real(r8) :: dom_to_acetate_q10
     real(r8) :: k_h2_methanogenesis_h2
     real(r8) :: k_h2_methanogenesis_co2
     real(r8) :: h2_methanogenesis_q10
     real(r8) :: k_acetoclastic_methanogenesis_acetate
     real(r8) :: acetoclastic_methanogenesis_q10
     real(r8) :: acetoclastic_methanogenesis_ch4_yield
     real(r8) :: k_aerobic_oxidation_ch4
     real(r8) :: k_aerobic_oxidation_o2
     real(r8) :: aerobic_oxidation_q10
     real(r8) :: aerobic_oxidation_o2_ch4_ratio
     real(r8) :: k_anaerobic_oxidation_ch4
     real(r8) :: anaerobic_oxidation_q10
     real(r8) :: aerobic_decomp_o2_c_ratio
     real(r8) :: ch4_transport_threshold
     real(r8) :: ph_min
     real(r8) :: ph_max
     real(r8) :: ph_opt
     real(r8) :: dom_diffusivity
     real(r8) :: aqueous_gas_diffusion_multiplier
     real(r8) :: plant_transport_coefficient
     real(r8) :: h2_plant_transport_threshold
     real(r8) :: atmospheric_ch4_mixing_ratio
     real(r8) :: atmospheric_o2_mixing_ratio
     real(r8) :: atmospheric_co2_mixing_ratio
     real(r8) :: atmospheric_h2_mixing_ratio
     real(r8) :: saturation_reaction_threshold
     real(r8) :: reaction_t_ref
     real(r8) :: aom_t_ref
     real(r8) :: acetate_ph_trigger
     real(r8) :: acidification_coefficient
     real(r8) :: acetate_feedback_half_saturation
     real(r8) :: soil_water_potential_min
     real(r8) :: h2_methanogenesis_co2_inhibition_scale
     real(r8) :: aom_o2_inhibition_scale
     real(r8) :: aerobic_acetate_oxidation_rate
     real(r8) :: ch4_h2_root_efold_depth
     real(r8) :: o2_root_efold_depth
     real(r8) :: ebullition_efold_depth
     real(r8) :: transport_thaw_threshold
     real(r8) :: plant_o2_consumption_fraction
     real(r8) :: plant_co2_flux_fraction
     real(r8) :: aqueous_diffusion_t_ref
     real(r8) :: aqueous_diffusion_temperature_exponent
  end type MicrobeMethaneParamsType

  type(MicrobeMethaneParamsType), public :: MicrobeMethaneParamsInst
  logical, public :: microbe_methane_parameters_read = .false.

  public :: readMicrobeMethaneParams
  public :: assertMicrobeMethaneParamsRead

contains

  !-----------------------------------------------------------------------
  subroutine readMicrobeMethaneParams(ncid)
    use ncdio_pio, only : file_desc_t, ncd_io

    type(file_desc_t), intent(inout) :: ncid

    call read_scalar(ncid, 'microbe_methane_mfg_biomass_min', MicrobeMethaneParamsInst%mfg_biomass_min)
    call read_scalar(ncid, 'microbe_methane_k_acetate', MicrobeMethaneParamsInst%k_acetate)
    call read_scalar(ncid, 'microbe_methane_acetate_prod_max', MicrobeMethaneParamsInst%acetate_prod_max)
    call read_scalar(ncid, 'microbe_methane_k_acetate_prod_o2', MicrobeMethaneParamsInst%k_acetate_prod_o2)
    call read_scalar(ncid, 'microbe_methane_acetogenesis_max', MicrobeMethaneParamsInst%acetogenesis_max)
    call read_scalar(ncid, 'microbe_methane_k_acetogenesis_h2', MicrobeMethaneParamsInst%k_acetogenesis_h2)
    call read_scalar(ncid, 'microbe_methane_k_acetogenesis_co2', MicrobeMethaneParamsInst%k_acetogenesis_co2)
    call read_scalar(ncid, 'microbe_methane_acetogenesis_q10', MicrobeMethaneParamsInst%acetogenesis_q10)
    call read_scalar(ncid, 'microbe_methane_h2_methanogen_growth_rate', &
         MicrobeMethaneParamsInst%h2_methanogen_growth_rate)
    call read_scalar(ncid, 'microbe_methane_h2_methanogen_death_rate', &
         MicrobeMethaneParamsInst%h2_methanogen_death_rate)
    call read_scalar(ncid, 'microbe_methane_h2_methanogen_yield', &
         MicrobeMethaneParamsInst%h2_methanogen_yield)
    call read_scalar(ncid, 'microbe_methane_acetate_methanogen_growth_rate', &
         MicrobeMethaneParamsInst%acetate_methanogen_growth_rate)
    call read_scalar(ncid, 'microbe_methane_acetate_methanogen_death_rate', &
         MicrobeMethaneParamsInst%acetate_methanogen_death_rate)
    call read_scalar(ncid, 'microbe_methane_acetate_methanogen_yield', &
         MicrobeMethaneParamsInst%acetate_methanogen_yield)
    call read_scalar(ncid, 'microbe_methane_aerobic_methanotroph_growth_rate', &
         MicrobeMethaneParamsInst%aerobic_methanotroph_growth_rate)
    call read_scalar(ncid, 'microbe_methane_aerobic_methanotroph_death_rate', &
         MicrobeMethaneParamsInst%aerobic_methanotroph_death_rate)
    call read_scalar(ncid, 'microbe_methane_aerobic_methanotroph_yield', &
         MicrobeMethaneParamsInst%aerobic_methanotroph_yield)
    call read_scalar(ncid, 'microbe_methane_anaerobic_methanotroph_growth_rate', &
         MicrobeMethaneParamsInst%anaerobic_methanotroph_growth_rate)
    call read_scalar(ncid, 'microbe_methane_anaerobic_methanotroph_death_rate', &
         MicrobeMethaneParamsInst%anaerobic_methanotroph_death_rate)
    call read_scalar(ncid, 'microbe_methane_anaerobic_methanotroph_yield', &
         MicrobeMethaneParamsInst%anaerobic_methanotroph_yield)
    call read_scalar(ncid, 'microbe_methane_dom_to_acetate_q10', MicrobeMethaneParamsInst%dom_to_acetate_q10)
    call read_scalar(ncid, 'microbe_methane_k_h2_methanogenesis_h2', &
         MicrobeMethaneParamsInst%k_h2_methanogenesis_h2)
    call read_scalar(ncid, 'microbe_methane_k_h2_methanogenesis_co2', &
         MicrobeMethaneParamsInst%k_h2_methanogenesis_co2)
    call read_scalar(ncid, 'microbe_methane_h2_methanogenesis_q10', &
         MicrobeMethaneParamsInst%h2_methanogenesis_q10)
    call read_scalar(ncid, 'microbe_methane_k_acetoclastic_methanogenesis_acetate', &
         MicrobeMethaneParamsInst%k_acetoclastic_methanogenesis_acetate)
    call read_scalar(ncid, 'microbe_methane_acetoclastic_methanogenesis_q10', &
         MicrobeMethaneParamsInst%acetoclastic_methanogenesis_q10)
    call read_scalar(ncid, 'microbe_methane_acetoclastic_methanogenesis_ch4_yield', &
         MicrobeMethaneParamsInst%acetoclastic_methanogenesis_ch4_yield)
    call read_scalar(ncid, 'microbe_methane_k_aerobic_oxidation_ch4', &
         MicrobeMethaneParamsInst%k_aerobic_oxidation_ch4)
    call read_scalar(ncid, 'microbe_methane_k_aerobic_oxidation_o2', &
         MicrobeMethaneParamsInst%k_aerobic_oxidation_o2)
    call read_scalar(ncid, 'microbe_methane_aerobic_oxidation_q10', MicrobeMethaneParamsInst%aerobic_oxidation_q10)
    call read_scalar(ncid, 'microbe_methane_aerobic_oxidation_o2_ch4_ratio', &
         MicrobeMethaneParamsInst%aerobic_oxidation_o2_ch4_ratio)
    call read_scalar(ncid, 'microbe_methane_k_anaerobic_oxidation_ch4', &
         MicrobeMethaneParamsInst%k_anaerobic_oxidation_ch4)
    call read_scalar(ncid, 'microbe_methane_anaerobic_oxidation_q10', &
         MicrobeMethaneParamsInst%anaerobic_oxidation_q10)
    call read_scalar(ncid, 'microbe_methane_aerobic_decomp_o2_c_ratio', &
         MicrobeMethaneParamsInst%aerobic_decomp_o2_c_ratio)
    call read_scalar(ncid, 'microbe_methane_ch4_transport_threshold', &
         MicrobeMethaneParamsInst%ch4_transport_threshold)
    call read_scalar(ncid, 'microbe_methane_ph_min', MicrobeMethaneParamsInst%ph_min)
    call read_scalar(ncid, 'microbe_methane_ph_max', MicrobeMethaneParamsInst%ph_max)
    call read_scalar(ncid, 'microbe_methane_ph_opt', MicrobeMethaneParamsInst%ph_opt)
    call read_scalar(ncid, 'microbe_methane_dom_diffusivity', MicrobeMethaneParamsInst%dom_diffusivity)
    call read_scalar(ncid, 'microbe_methane_aqueous_gas_diffusion_multiplier', &
         MicrobeMethaneParamsInst%aqueous_gas_diffusion_multiplier)
    call read_scalar(ncid, 'microbe_methane_plant_transport_coefficient', &
         MicrobeMethaneParamsInst%plant_transport_coefficient)
    call read_scalar(ncid, 'microbe_methane_h2_plant_transport_threshold', &
         MicrobeMethaneParamsInst%h2_plant_transport_threshold)
    call read_scalar(ncid, 'microbe_methane_atmospheric_ch4_mixing_ratio', &
         MicrobeMethaneParamsInst%atmospheric_ch4_mixing_ratio)
    call read_scalar(ncid, 'microbe_methane_atmospheric_o2_mixing_ratio', &
         MicrobeMethaneParamsInst%atmospheric_o2_mixing_ratio)
    call read_scalar(ncid, 'microbe_methane_atmospheric_co2_mixing_ratio', &
         MicrobeMethaneParamsInst%atmospheric_co2_mixing_ratio)
    call read_scalar(ncid, 'microbe_methane_atmospheric_h2_mixing_ratio', &
         MicrobeMethaneParamsInst%atmospheric_h2_mixing_ratio)
    call read_scalar(ncid, 'microbe_methane_saturation_reaction_threshold', &
         MicrobeMethaneParamsInst%saturation_reaction_threshold)
    call read_scalar(ncid, 'microbe_methane_reaction_t_ref', MicrobeMethaneParamsInst%reaction_t_ref)
    call read_scalar(ncid, 'microbe_methane_aom_t_ref', MicrobeMethaneParamsInst%aom_t_ref)
    call read_scalar(ncid, 'microbe_methane_acetate_ph_trigger', MicrobeMethaneParamsInst%acetate_ph_trigger)
    call read_scalar(ncid, 'microbe_methane_acidification_coefficient', &
         MicrobeMethaneParamsInst%acidification_coefficient)
    call read_scalar(ncid, 'microbe_methane_acetate_feedback_half_saturation', &
         MicrobeMethaneParamsInst%acetate_feedback_half_saturation)
    call read_scalar(ncid, 'microbe_methane_soil_water_potential_min', &
         MicrobeMethaneParamsInst%soil_water_potential_min)
    call read_scalar(ncid, 'microbe_methane_h2_methanogenesis_co2_inhibition_scale', &
         MicrobeMethaneParamsInst%h2_methanogenesis_co2_inhibition_scale)
    call read_scalar(ncid, 'microbe_methane_aom_o2_inhibition_scale', &
         MicrobeMethaneParamsInst%aom_o2_inhibition_scale)
    call read_scalar(ncid, 'microbe_methane_aerobic_acetate_oxidation_rate', &
         MicrobeMethaneParamsInst%aerobic_acetate_oxidation_rate)
    call read_scalar(ncid, 'microbe_methane_ch4_h2_root_efold_depth', &
         MicrobeMethaneParamsInst%ch4_h2_root_efold_depth)
    call read_scalar(ncid, 'microbe_methane_o2_root_efold_depth', MicrobeMethaneParamsInst%o2_root_efold_depth)
    call read_scalar(ncid, 'microbe_methane_ebullition_efold_depth', &
         MicrobeMethaneParamsInst%ebullition_efold_depth)
    call read_scalar(ncid, 'microbe_methane_transport_thaw_threshold', &
         MicrobeMethaneParamsInst%transport_thaw_threshold)
    call read_scalar(ncid, 'microbe_methane_plant_o2_consumption_fraction', &
         MicrobeMethaneParamsInst%plant_o2_consumption_fraction)
    call read_scalar(ncid, 'microbe_methane_plant_co2_flux_fraction', &
         MicrobeMethaneParamsInst%plant_co2_flux_fraction)
    call read_scalar(ncid, 'microbe_methane_aqueous_diffusion_t_ref', &
         MicrobeMethaneParamsInst%aqueous_diffusion_t_ref)
    call read_scalar(ncid, 'microbe_methane_aqueous_diffusion_temperature_exponent', &
         MicrobeMethaneParamsInst%aqueous_diffusion_temperature_exponent)

    call validate_parameters()
    microbe_methane_parameters_read = .true.

  contains

    subroutine read_scalar(file, name, value)
      type(file_desc_t), intent(inout) :: file
      character(len=*), intent(in) :: name
      real(r8), intent(out) :: value
      logical :: readv

      call ncd_io(trim(name), value, 'read', file, readvar=readv)
      if (.not. readv) then
         call endrun(msg=' ERROR: revised methane parameter is missing: '//trim(name)//&
              errMsg(__FILE__, __LINE__))
      end if
    end subroutine read_scalar

  end subroutine readMicrobeMethaneParams

  !-----------------------------------------------------------------------
  subroutine assertMicrobeMethaneParamsRead()
    if (.not. microbe_methane_parameters_read) then
       call endrun(msg=' ERROR: revised methane parameters were not read before state initialization'//&
            errMsg(__FILE__, __LINE__))
    end if
  end subroutine assertMicrobeMethaneParamsRead

  !-----------------------------------------------------------------------
  subroutine validate_parameters()
    call require_nonnegative('mfg_biomass_min', MicrobeMethaneParamsInst%mfg_biomass_min)
    call require_positive('k_acetate', MicrobeMethaneParamsInst%k_acetate)
    call require_nonnegative('acetate_prod_max', MicrobeMethaneParamsInst%acetate_prod_max)
    call require_positive('k_acetate_prod_o2', MicrobeMethaneParamsInst%k_acetate_prod_o2)
    call require_nonnegative('acetogenesis_max', MicrobeMethaneParamsInst%acetogenesis_max)
    call require_positive('k_acetogenesis_h2', MicrobeMethaneParamsInst%k_acetogenesis_h2)
    call require_positive('k_acetogenesis_co2', MicrobeMethaneParamsInst%k_acetogenesis_co2)
    call require_positive('acetogenesis_q10', MicrobeMethaneParamsInst%acetogenesis_q10)
    call require_nonnegative('h2_methanogen_growth_rate', MicrobeMethaneParamsInst%h2_methanogen_growth_rate)
    call require_nonnegative('h2_methanogen_death_rate', MicrobeMethaneParamsInst%h2_methanogen_death_rate)
    call require_positive_fraction('h2_methanogen_yield', MicrobeMethaneParamsInst%h2_methanogen_yield)
    call require_nonnegative('acetate_methanogen_growth_rate', MicrobeMethaneParamsInst%acetate_methanogen_growth_rate)
    call require_nonnegative('acetate_methanogen_death_rate', MicrobeMethaneParamsInst%acetate_methanogen_death_rate)
    call require_positive_fraction('acetate_methanogen_yield', MicrobeMethaneParamsInst%acetate_methanogen_yield)
    call require_nonnegative('aerobic_methanotroph_growth_rate', &
         MicrobeMethaneParamsInst%aerobic_methanotroph_growth_rate)
    call require_nonnegative('aerobic_methanotroph_death_rate', &
         MicrobeMethaneParamsInst%aerobic_methanotroph_death_rate)
    call require_positive_fraction('aerobic_methanotroph_yield', &
         MicrobeMethaneParamsInst%aerobic_methanotroph_yield)
    call require_nonnegative('anaerobic_methanotroph_growth_rate', &
         MicrobeMethaneParamsInst%anaerobic_methanotroph_growth_rate)
    call require_nonnegative('anaerobic_methanotroph_death_rate', &
         MicrobeMethaneParamsInst%anaerobic_methanotroph_death_rate)
    call require_positive_fraction('anaerobic_methanotroph_yield', &
         MicrobeMethaneParamsInst%anaerobic_methanotroph_yield)
    call require_positive('dom_to_acetate_q10', MicrobeMethaneParamsInst%dom_to_acetate_q10)
    call require_positive('k_h2_methanogenesis_h2', MicrobeMethaneParamsInst%k_h2_methanogenesis_h2)
    call require_positive('k_h2_methanogenesis_co2', MicrobeMethaneParamsInst%k_h2_methanogenesis_co2)
    call require_positive('h2_methanogenesis_q10', MicrobeMethaneParamsInst%h2_methanogenesis_q10)
    call require_positive('k_acetoclastic_methanogenesis_acetate', &
         MicrobeMethaneParamsInst%k_acetoclastic_methanogenesis_acetate)
    call require_positive('acetoclastic_methanogenesis_q10', &
         MicrobeMethaneParamsInst%acetoclastic_methanogenesis_q10)
    call require_fraction('acetoclastic_methanogenesis_ch4_yield', &
         MicrobeMethaneParamsInst%acetoclastic_methanogenesis_ch4_yield)
    call require_positive('k_aerobic_oxidation_ch4', MicrobeMethaneParamsInst%k_aerobic_oxidation_ch4)
    call require_positive('k_aerobic_oxidation_o2', MicrobeMethaneParamsInst%k_aerobic_oxidation_o2)
    call require_positive('aerobic_oxidation_q10', MicrobeMethaneParamsInst%aerobic_oxidation_q10)
    call require_positive('aerobic_oxidation_o2_ch4_ratio', &
         MicrobeMethaneParamsInst%aerobic_oxidation_o2_ch4_ratio)
    call require_positive('k_anaerobic_oxidation_ch4', MicrobeMethaneParamsInst%k_anaerobic_oxidation_ch4)
    call require_positive('anaerobic_oxidation_q10', MicrobeMethaneParamsInst%anaerobic_oxidation_q10)
    call require_positive('aerobic_decomp_o2_c_ratio', MicrobeMethaneParamsInst%aerobic_decomp_o2_c_ratio)
    call require_nonnegative('ch4_transport_threshold', MicrobeMethaneParamsInst%ch4_transport_threshold)
    call require_finite('ph_min', MicrobeMethaneParamsInst%ph_min)
    call require_finite('ph_max', MicrobeMethaneParamsInst%ph_max)
    call require_finite('ph_opt', MicrobeMethaneParamsInst%ph_opt)
    if (.not. (MicrobeMethaneParamsInst%ph_min < MicrobeMethaneParamsInst%ph_opt .and. &
         MicrobeMethaneParamsInst%ph_opt < MicrobeMethaneParamsInst%ph_max)) then
       call parameter_error('pH response must satisfy ph_min < ph_opt < ph_max')
    end if
    call require_positive('dom_diffusivity', MicrobeMethaneParamsInst%dom_diffusivity)
    call require_positive('aqueous_gas_diffusion_multiplier', &
         MicrobeMethaneParamsInst%aqueous_gas_diffusion_multiplier)
    call require_nonnegative('plant_transport_coefficient', MicrobeMethaneParamsInst%plant_transport_coefficient)
    call require_nonnegative('h2_plant_transport_threshold', &
         MicrobeMethaneParamsInst%h2_plant_transport_threshold)
    call require_fraction('atmospheric_ch4_mixing_ratio', MicrobeMethaneParamsInst%atmospheric_ch4_mixing_ratio)
    call require_fraction('atmospheric_o2_mixing_ratio', MicrobeMethaneParamsInst%atmospheric_o2_mixing_ratio)
    call require_fraction('atmospheric_co2_mixing_ratio', MicrobeMethaneParamsInst%atmospheric_co2_mixing_ratio)
    call require_fraction('atmospheric_h2_mixing_ratio', MicrobeMethaneParamsInst%atmospheric_h2_mixing_ratio)
    call require_fraction('saturation_reaction_threshold', MicrobeMethaneParamsInst%saturation_reaction_threshold)
    call require_positive('reaction_t_ref', MicrobeMethaneParamsInst%reaction_t_ref)
    call require_positive('aom_t_ref', MicrobeMethaneParamsInst%aom_t_ref)
    call require_finite('acetate_ph_trigger', MicrobeMethaneParamsInst%acetate_ph_trigger)
    call require_nonnegative('acidification_coefficient', MicrobeMethaneParamsInst%acidification_coefficient)
    call require_positive('acetate_feedback_half_saturation', &
         MicrobeMethaneParamsInst%acetate_feedback_half_saturation)
    call require_finite('soil_water_potential_min', MicrobeMethaneParamsInst%soil_water_potential_min)
    call require_positive('h2_methanogenesis_co2_inhibition_scale', &
         MicrobeMethaneParamsInst%h2_methanogenesis_co2_inhibition_scale)
    call require_positive('aom_o2_inhibition_scale', MicrobeMethaneParamsInst%aom_o2_inhibition_scale)
    call require_nonnegative('aerobic_acetate_oxidation_rate', &
         MicrobeMethaneParamsInst%aerobic_acetate_oxidation_rate)
    call require_positive('ch4_h2_root_efold_depth', MicrobeMethaneParamsInst%ch4_h2_root_efold_depth)
    call require_positive('o2_root_efold_depth', MicrobeMethaneParamsInst%o2_root_efold_depth)
    call require_positive('ebullition_efold_depth', MicrobeMethaneParamsInst%ebullition_efold_depth)
    call require_positive('transport_thaw_threshold', MicrobeMethaneParamsInst%transport_thaw_threshold)
    call require_fraction('plant_o2_consumption_fraction', MicrobeMethaneParamsInst%plant_o2_consumption_fraction)
    call require_fraction('plant_co2_flux_fraction', MicrobeMethaneParamsInst%plant_co2_flux_fraction)
    call require_positive('aqueous_diffusion_t_ref', MicrobeMethaneParamsInst%aqueous_diffusion_t_ref)
    call require_positive('aqueous_diffusion_temperature_exponent', &
         MicrobeMethaneParamsInst%aqueous_diffusion_temperature_exponent)
  end subroutine validate_parameters

  !-----------------------------------------------------------------------
  subroutine require_finite(name, value)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: value
    if (.not. ieee_is_finite(value)) call parameter_error(trim(name)//' must be finite')
  end subroutine require_finite

  subroutine require_nonnegative(name, value)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: value
    call require_finite(name, value)
    if (value < 0._r8) call parameter_error(trim(name)//' must be nonnegative')
  end subroutine require_nonnegative

  subroutine require_positive(name, value)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: value
    call require_finite(name, value)
    if (value <= 0._r8) call parameter_error(trim(name)//' must be positive')
  end subroutine require_positive

  subroutine require_fraction(name, value)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: value
    call require_finite(name, value)
    if (value < 0._r8 .or. value > 1._r8) then
       call parameter_error(trim(name)//' must be in [0,1]')
    end if
  end subroutine require_fraction

  subroutine require_positive_fraction(name, value)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: value
    call require_finite(name, value)
    if (value <= 0._r8 .or. value > 1._r8) then
       call parameter_error(trim(name)//' must be in (0,1]')
    end if
  end subroutine require_positive_fraction

  subroutine parameter_error(message)
    character(len=*), intent(in) :: message
    call endrun(msg=' ERROR: invalid revised methane parameter: '//trim(message)//&
         errMsg(__FILE__, __LINE__))
  end subroutine parameter_error

end module MicrobeMethaneParamsMod
