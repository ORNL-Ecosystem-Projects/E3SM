module MicrobeMethaneMod

  ! State and lifecycle owner for the revised microbial methane backend.
  ! The ELM driver selects this backend instead of CH4Mod when enabled.
  ! DOM, bacteria, and fungi remain authoritative decomposition pools.

  use shr_kind_mod, only : r8 => shr_kind_r8
  use decompMod, only : bounds_type
  use elm_varctl, only : use_microbe_methane, use_legacy_ch4_with_microbe, &
       use_elm_microbe_methane_transport, &
       use_microbe_nonbog_lateral_gas_transport, &
       use_clm_microbe_humhol_saturation, use_clm_microbe_dom_relaxation, &
       use_microbe_aqueous_transport, use_humhol
  use pio, only : file_desc_t

  implicit none
  private
  save

  type, public :: microbe_methane_type
     ! Carbon-bearing state is g C m-3 soil; gas state is mol gas m-3 bulk soil.
     ! Partition values are concentrations, not area-weighted bulk copies.
     real(r8), pointer :: acetate_c_unsat_col(:,:) => null()
     real(r8), pointer :: acetate_c_sat_col(:,:) => null()
     real(r8), pointer :: acetate_methanogen_c_unsat_col(:,:) => null()
     real(r8), pointer :: acetate_methanogen_c_sat_col(:,:) => null()
     real(r8), pointer :: h2_methanogen_c_unsat_col(:,:) => null()
     real(r8), pointer :: h2_methanogen_c_sat_col(:,:) => null()
     real(r8), pointer :: aerobic_methanotroph_c_unsat_col(:,:) => null()
     real(r8), pointer :: aerobic_methanotroph_c_sat_col(:,:) => null()
     real(r8), pointer :: anaerobic_methanotroph_c_unsat_col(:,:) => null()
     real(r8), pointer :: anaerobic_methanotroph_c_sat_col(:,:) => null()
     real(r8), pointer :: conc_ch4_unsat_col(:,:) => null()
     real(r8), pointer :: conc_ch4_sat_col(:,:) => null()
     real(r8), pointer :: conc_o2_unsat_col(:,:) => null()
     real(r8), pointer :: conc_o2_sat_col(:,:) => null()
     real(r8), pointer :: conc_co2_unsat_col(:,:) => null()
     real(r8), pointer :: conc_co2_sat_col(:,:) => null()
     real(r8), pointer :: conc_h2_unsat_col(:,:) => null()
     real(r8), pointer :: conc_h2_sat_col(:,:) => null()
     ! Diagnostic shared-O2 limiter before ELM's 0.2 decomposition floor.
     real(r8), pointer :: o2_stress_unsat_col(:,:) => null()
     real(r8), pointer :: o2_stress_sat_col(:,:) => null()

     ! Future conservative partition remapping uses this prior-step fraction.
     real(r8), pointer :: sat_fraction_previous_col(:) => null()

     ! Adapter-owned accounting fields. DOM is deliberately excluded from
     ! additional_carbon_col because it is already in ELM's decomp pools.
     real(r8), pointer :: additional_carbon_col(:) => null()
     real(r8), pointer :: surface_carbon_flux_col(:) => null()
     ! Net lateral CH4 plus CO2 carbon transfer, positive into the column.
     real(r8), pointer :: lateral_carbon_flux_col(:) => null()
     ! Physical aqueous-transport diagnostics. Interface arrays store the
     ! downward-positive flux through the lower face of each decomposition layer.
     real(r8), pointer :: dom_advective_flux_col(:,:) => null()
     real(r8), pointer :: dom_diffusive_flux_col(:,:) => null()
     real(r8), pointer :: aqueous_carbon_export_col(:) => null()
     real(r8), pointer :: aqueous_nitrogen_export_col(:) => null()
     real(r8), pointer :: aqueous_phosphorus_export_col(:) => null()
     real(r8), pointer :: dom_bottom_carbon_export_col(:) => null()
     real(r8), pointer :: dom_bottom_nitrogen_export_col(:) => null()
     real(r8), pointer :: dom_bottom_phosphorus_export_col(:) => null()
     real(r8), pointer :: dom_carbon_transport_residual_col(:) => null()
     real(r8), pointer :: dom_nitrogen_transport_residual_col(:) => null()
     real(r8), pointer :: dom_phosphorus_transport_residual_col(:) => null()
     real(r8), pointer :: surface_ch4_flux_col(:) => null()
     real(r8), pointer :: surface_co2_flux_col(:) => null()
     ! Area-weighted CH4 surface pathways in g C m-2 s-1, positive to air.
     real(r8), pointer :: ch4_surface_diffusion_col(:) => null()
     real(r8), pointer :: ch4_surface_diffusion_unsat_col(:) => null()
     real(r8), pointer :: ch4_surface_diffusion_sat_col(:) => null()
     real(r8), pointer :: ch4_surface_aerenchyma_col(:) => null()
     real(r8), pointer :: ch4_surface_aerenchyma_unsat_col(:) => null()
     real(r8), pointer :: ch4_surface_aerenchyma_sat_col(:) => null()
     real(r8), pointer :: ch4_surface_ebullition_col(:) => null()
     real(r8), pointer :: ch4_surface_ebullition_unsat_col(:) => null()
     real(r8), pointer :: ch4_surface_ebullition_sat_col(:) => null()
     real(r8), pointer :: ch4_production_col(:) => null()
     real(r8), pointer :: ch4_oxidation_col(:) => null()
     ! Area-weighted contributions; each pair sums to its bulk diagnostic.
     real(r8), pointer :: ch4_production_unsat_col(:) => null()
     real(r8), pointer :: ch4_production_sat_col(:) => null()
     real(r8), pointer :: ch4_oxidation_unsat_col(:) => null()
     real(r8), pointer :: ch4_oxidation_sat_col(:) => null()
     real(r8), pointer :: ch4_aerobic_oxidation_col(:) => null()
     real(r8), pointer :: ch4_aerobic_oxidation_unsat_col(:) => null()
     real(r8), pointer :: ch4_aerobic_oxidation_sat_col(:) => null()
     real(r8), pointer :: ch4_aerobic_oxidation_pre_o2_col(:) => null()
     real(r8), pointer :: ch4_aerobic_oxidation_pre_o2_unsat_col(:) => null()
     real(r8), pointer :: ch4_aerobic_oxidation_pre_o2_sat_col(:) => null()
     real(r8), pointer :: ch4_anaerobic_oxidation_col(:) => null()
     real(r8), pointer :: ch4_anaerobic_oxidation_unsat_col(:) => null()
     real(r8), pointer :: ch4_anaerobic_oxidation_sat_col(:) => null()
   contains
     procedure, public  :: Init
     procedure, private :: InitAllocate
     procedure, private :: InitCold
     procedure, private :: InitHistory
     procedure, private :: ReadParams
     procedure, public  :: Repartition
     procedure, public  :: Advance
     procedure, private :: UpdateAdditionalCarbon
     procedure, public  :: Restart
  end type microbe_methane_type

contains

  subroutine Init(this, bounds)
    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds

    if (.not. use_microbe_methane .or. use_legacy_ch4_with_microbe) return
    call this%ReadParams()
    call this%InitAllocate(bounds)
    call this%InitCold(bounds)
    call this%InitHistory(bounds)
  end subroutine Init

  subroutine InitAllocate(this, bounds)
    use elm_varpar, only : nlevdecomp_full
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)

    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds
    integer :: begc, endc

    begc = bounds%begc
    endc = bounds%endc
    allocate(this%acetate_c_unsat_col(begc:endc,1:nlevdecomp_full)); this%acetate_c_unsat_col = nan
    allocate(this%acetate_c_sat_col(begc:endc,1:nlevdecomp_full)); this%acetate_c_sat_col = nan
    allocate(this%acetate_methanogen_c_unsat_col(begc:endc,1:nlevdecomp_full)); &
         this%acetate_methanogen_c_unsat_col = nan
    allocate(this%acetate_methanogen_c_sat_col(begc:endc,1:nlevdecomp_full)); &
         this%acetate_methanogen_c_sat_col = nan
    allocate(this%h2_methanogen_c_unsat_col(begc:endc,1:nlevdecomp_full)); this%h2_methanogen_c_unsat_col = nan
    allocate(this%h2_methanogen_c_sat_col(begc:endc,1:nlevdecomp_full)); this%h2_methanogen_c_sat_col = nan
    allocate(this%aerobic_methanotroph_c_unsat_col(begc:endc,1:nlevdecomp_full)); &
         this%aerobic_methanotroph_c_unsat_col = nan
    allocate(this%aerobic_methanotroph_c_sat_col(begc:endc,1:nlevdecomp_full)); &
         this%aerobic_methanotroph_c_sat_col = nan
    allocate(this%anaerobic_methanotroph_c_unsat_col(begc:endc,1:nlevdecomp_full)); &
         this%anaerobic_methanotroph_c_unsat_col = nan
    allocate(this%anaerobic_methanotroph_c_sat_col(begc:endc,1:nlevdecomp_full)); &
         this%anaerobic_methanotroph_c_sat_col = nan
    allocate(this%conc_ch4_unsat_col(begc:endc,1:nlevdecomp_full)); this%conc_ch4_unsat_col = nan
    allocate(this%conc_ch4_sat_col(begc:endc,1:nlevdecomp_full)); this%conc_ch4_sat_col = nan
    allocate(this%conc_o2_unsat_col(begc:endc,1:nlevdecomp_full)); this%conc_o2_unsat_col = nan
    allocate(this%conc_o2_sat_col(begc:endc,1:nlevdecomp_full)); this%conc_o2_sat_col = nan
    allocate(this%conc_co2_unsat_col(begc:endc,1:nlevdecomp_full)); this%conc_co2_unsat_col = nan
    allocate(this%conc_co2_sat_col(begc:endc,1:nlevdecomp_full)); this%conc_co2_sat_col = nan
    allocate(this%conc_h2_unsat_col(begc:endc,1:nlevdecomp_full)); this%conc_h2_unsat_col = nan
    allocate(this%conc_h2_sat_col(begc:endc,1:nlevdecomp_full)); this%conc_h2_sat_col = nan
    allocate(this%o2_stress_unsat_col(begc:endc,1:nlevdecomp_full)); this%o2_stress_unsat_col = nan
    allocate(this%o2_stress_sat_col(begc:endc,1:nlevdecomp_full)); this%o2_stress_sat_col = nan
    allocate(this%sat_fraction_previous_col(begc:endc)); this%sat_fraction_previous_col = nan
    allocate(this%additional_carbon_col(begc:endc)); this%additional_carbon_col = nan
    allocate(this%surface_carbon_flux_col(begc:endc)); this%surface_carbon_flux_col = nan
    allocate(this%lateral_carbon_flux_col(begc:endc)); this%lateral_carbon_flux_col = nan
    allocate(this%dom_advective_flux_col(begc:endc,1:nlevdecomp_full)); &
         this%dom_advective_flux_col = nan
    allocate(this%dom_diffusive_flux_col(begc:endc,1:nlevdecomp_full)); &
         this%dom_diffusive_flux_col = nan
    allocate(this%aqueous_carbon_export_col(begc:endc)); this%aqueous_carbon_export_col = nan
    allocate(this%aqueous_nitrogen_export_col(begc:endc)); this%aqueous_nitrogen_export_col = nan
    allocate(this%aqueous_phosphorus_export_col(begc:endc)); &
         this%aqueous_phosphorus_export_col = nan
    allocate(this%dom_bottom_carbon_export_col(begc:endc)); &
         this%dom_bottom_carbon_export_col = nan
    allocate(this%dom_bottom_nitrogen_export_col(begc:endc)); &
         this%dom_bottom_nitrogen_export_col = nan
    allocate(this%dom_bottom_phosphorus_export_col(begc:endc)); &
         this%dom_bottom_phosphorus_export_col = nan
    allocate(this%dom_carbon_transport_residual_col(begc:endc)); &
         this%dom_carbon_transport_residual_col = nan
    allocate(this%dom_nitrogen_transport_residual_col(begc:endc)); &
         this%dom_nitrogen_transport_residual_col = nan
    allocate(this%dom_phosphorus_transport_residual_col(begc:endc)); &
         this%dom_phosphorus_transport_residual_col = nan
    allocate(this%surface_ch4_flux_col(begc:endc)); this%surface_ch4_flux_col = nan
    allocate(this%surface_co2_flux_col(begc:endc)); this%surface_co2_flux_col = nan
    allocate(this%ch4_surface_diffusion_col(begc:endc)); this%ch4_surface_diffusion_col = nan
    allocate(this%ch4_surface_diffusion_unsat_col(begc:endc)); &
         this%ch4_surface_diffusion_unsat_col = nan
    allocate(this%ch4_surface_diffusion_sat_col(begc:endc)); this%ch4_surface_diffusion_sat_col = nan
    allocate(this%ch4_surface_aerenchyma_col(begc:endc)); this%ch4_surface_aerenchyma_col = nan
    allocate(this%ch4_surface_aerenchyma_unsat_col(begc:endc)); &
         this%ch4_surface_aerenchyma_unsat_col = nan
    allocate(this%ch4_surface_aerenchyma_sat_col(begc:endc)); &
         this%ch4_surface_aerenchyma_sat_col = nan
    allocate(this%ch4_surface_ebullition_col(begc:endc)); this%ch4_surface_ebullition_col = nan
    allocate(this%ch4_surface_ebullition_unsat_col(begc:endc)); &
         this%ch4_surface_ebullition_unsat_col = nan
    allocate(this%ch4_surface_ebullition_sat_col(begc:endc)); this%ch4_surface_ebullition_sat_col = nan
    allocate(this%ch4_production_col(begc:endc)); this%ch4_production_col = nan
    allocate(this%ch4_oxidation_col(begc:endc)); this%ch4_oxidation_col = nan
    allocate(this%ch4_production_unsat_col(begc:endc)); this%ch4_production_unsat_col = nan
    allocate(this%ch4_production_sat_col(begc:endc)); this%ch4_production_sat_col = nan
    allocate(this%ch4_oxidation_unsat_col(begc:endc)); this%ch4_oxidation_unsat_col = nan
    allocate(this%ch4_oxidation_sat_col(begc:endc)); this%ch4_oxidation_sat_col = nan
    allocate(this%ch4_aerobic_oxidation_col(begc:endc)); this%ch4_aerobic_oxidation_col = nan
    allocate(this%ch4_aerobic_oxidation_unsat_col(begc:endc)); this%ch4_aerobic_oxidation_unsat_col = nan
    allocate(this%ch4_aerobic_oxidation_sat_col(begc:endc)); this%ch4_aerobic_oxidation_sat_col = nan
    allocate(this%ch4_aerobic_oxidation_pre_o2_col(begc:endc)); &
         this%ch4_aerobic_oxidation_pre_o2_col = nan
    allocate(this%ch4_aerobic_oxidation_pre_o2_unsat_col(begc:endc)); &
         this%ch4_aerobic_oxidation_pre_o2_unsat_col = nan
    allocate(this%ch4_aerobic_oxidation_pre_o2_sat_col(begc:endc)); &
         this%ch4_aerobic_oxidation_pre_o2_sat_col = nan
    allocate(this%ch4_anaerobic_oxidation_col(begc:endc)); this%ch4_anaerobic_oxidation_col = nan
    allocate(this%ch4_anaerobic_oxidation_unsat_col(begc:endc)); &
         this%ch4_anaerobic_oxidation_unsat_col = nan
    allocate(this%ch4_anaerobic_oxidation_sat_col(begc:endc)); this%ch4_anaerobic_oxidation_sat_col = nan
  end subroutine InitAllocate

  subroutine ReadParams(this)
    use MicrobeMethaneParamsMod, only : assertMicrobeMethaneParamsRead
    class(microbe_methane_type) :: this
    call assertMicrobeMethaneParamsRead()
  end subroutine ReadParams

  subroutine InitCold(this, bounds)
    use elm_varpar, only : nlevdecomp
    use ColumnType, only : col_pp
    use MicrobeMethaneParamsMod, only : MicrobeMethaneParamsInst

    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds
    integer :: c
    real(r8) :: biomass_seed

    biomass_seed = MicrobeMethaneParamsInst%mfg_biomass_min
    ! History aggregation can inspect the complete allocated decomposition
    ! dimension, including inactive columns and levels below nlevdecomp. Keep
    ! those entries finite and carbon-neutral; active soil levels are seeded
    ! below.
    this%acetate_c_unsat_col = 0._r8
    this%acetate_c_sat_col = 0._r8
    this%acetate_methanogen_c_unsat_col = 0._r8
    this%acetate_methanogen_c_sat_col = 0._r8
    this%h2_methanogen_c_unsat_col = 0._r8
    this%h2_methanogen_c_sat_col = 0._r8
    this%aerobic_methanotroph_c_unsat_col = 0._r8
    this%aerobic_methanotroph_c_sat_col = 0._r8
    this%anaerobic_methanotroph_c_unsat_col = 0._r8
    this%anaerobic_methanotroph_c_sat_col = 0._r8
    this%conc_ch4_unsat_col = 0._r8
    this%conc_ch4_sat_col = 0._r8
    this%conc_o2_unsat_col = 0._r8
    this%conc_o2_sat_col = 0._r8
    this%conc_co2_unsat_col = 0._r8
    this%conc_co2_sat_col = 0._r8
    this%conc_h2_unsat_col = 0._r8
    this%conc_h2_sat_col = 0._r8
    this%o2_stress_unsat_col = 0._r8
    this%o2_stress_sat_col = 0._r8
    this%sat_fraction_previous_col = 0._r8
    this%additional_carbon_col = 0._r8
    this%surface_carbon_flux_col = 0._r8
    this%lateral_carbon_flux_col = 0._r8
    this%dom_advective_flux_col = 0._r8
    this%dom_diffusive_flux_col = 0._r8
    this%aqueous_carbon_export_col = 0._r8
    this%aqueous_nitrogen_export_col = 0._r8
    this%aqueous_phosphorus_export_col = 0._r8
    this%dom_bottom_carbon_export_col = 0._r8
    this%dom_bottom_nitrogen_export_col = 0._r8
    this%dom_bottom_phosphorus_export_col = 0._r8
    this%dom_carbon_transport_residual_col = 0._r8
    this%dom_nitrogen_transport_residual_col = 0._r8
    this%dom_phosphorus_transport_residual_col = 0._r8
    this%surface_ch4_flux_col = 0._r8
    this%surface_co2_flux_col = 0._r8
    this%ch4_surface_diffusion_col = 0._r8
    this%ch4_surface_diffusion_unsat_col = 0._r8
    this%ch4_surface_diffusion_sat_col = 0._r8
    this%ch4_surface_aerenchyma_col = 0._r8
    this%ch4_surface_aerenchyma_unsat_col = 0._r8
    this%ch4_surface_aerenchyma_sat_col = 0._r8
    this%ch4_surface_ebullition_col = 0._r8
    this%ch4_surface_ebullition_unsat_col = 0._r8
    this%ch4_surface_ebullition_sat_col = 0._r8
    this%ch4_production_col = 0._r8
    this%ch4_oxidation_col = 0._r8
    this%ch4_production_unsat_col = 0._r8
    this%ch4_production_sat_col = 0._r8
    this%ch4_oxidation_unsat_col = 0._r8
    this%ch4_oxidation_sat_col = 0._r8
    this%ch4_aerobic_oxidation_col = 0._r8
    this%ch4_aerobic_oxidation_unsat_col = 0._r8
    this%ch4_aerobic_oxidation_sat_col = 0._r8
    this%ch4_aerobic_oxidation_pre_o2_col = 0._r8
    this%ch4_aerobic_oxidation_pre_o2_unsat_col = 0._r8
    this%ch4_aerobic_oxidation_pre_o2_sat_col = 0._r8
    this%ch4_anaerobic_oxidation_col = 0._r8
    this%ch4_anaerobic_oxidation_unsat_col = 0._r8
    this%ch4_anaerobic_oxidation_sat_col = 0._r8

    do c = bounds%begc, bounds%endc
       if (col_pp%is_soil(c) .or. col_pp%is_crop(c)) then
          this%acetate_c_unsat_col(c,1:nlevdecomp) = 0._r8
          this%acetate_c_sat_col(c,1:nlevdecomp) = 0._r8
          this%acetate_methanogen_c_unsat_col(c,1:nlevdecomp) = biomass_seed
          this%acetate_methanogen_c_sat_col(c,1:nlevdecomp) = biomass_seed
          this%h2_methanogen_c_unsat_col(c,1:nlevdecomp) = biomass_seed
          this%h2_methanogen_c_sat_col(c,1:nlevdecomp) = biomass_seed
          this%aerobic_methanotroph_c_unsat_col(c,1:nlevdecomp) = biomass_seed
          this%aerobic_methanotroph_c_sat_col(c,1:nlevdecomp) = biomass_seed
          this%anaerobic_methanotroph_c_unsat_col(c,1:nlevdecomp) = biomass_seed
          this%anaerobic_methanotroph_c_sat_col(c,1:nlevdecomp) = biomass_seed
          this%conc_ch4_unsat_col(c,1:nlevdecomp) = 0._r8
          this%conc_ch4_sat_col(c,1:nlevdecomp) = 0._r8
          this%conc_o2_unsat_col(c,1:nlevdecomp) = 0._r8
          this%conc_o2_sat_col(c,1:nlevdecomp) = 0._r8
          this%conc_co2_unsat_col(c,1:nlevdecomp) = 0._r8
          this%conc_co2_sat_col(c,1:nlevdecomp) = 0._r8
          this%conc_h2_unsat_col(c,1:nlevdecomp) = 0._r8
          this%conc_h2_sat_col(c,1:nlevdecomp) = 0._r8
          this%o2_stress_unsat_col(c,1:nlevdecomp) = 0._r8
          this%o2_stress_sat_col(c,1:nlevdecomp) = 0._r8
          this%additional_carbon_col(c) = 0._r8
          this%surface_carbon_flux_col(c) = 0._r8
          this%lateral_carbon_flux_col(c) = 0._r8
          this%surface_ch4_flux_col(c) = 0._r8
          this%surface_co2_flux_col(c) = 0._r8
       end if
    end do
    call this%UpdateAdditionalCarbon(bounds)
  end subroutine InitCold

  subroutine InitHistory(this, bounds)
    use histFileMod, only : hist_addfld_decomp, hist_addfld1d

    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds

    call add_state('MM_ACETATE_C_UNSAT', 'gC/m^3', &
         'acetate carbon concentration in unsaturated subarea', this%acetate_c_unsat_col)
    call add_state('MM_ACETATE_C_SAT', 'gC/m^3', &
         'acetate carbon concentration in saturated subarea', this%acetate_c_sat_col)
    call add_state('MM_ACET_METH_C_UNSAT', 'gC/m^3', &
         'acetoclastic methanogen biomass carbon in unsaturated subarea', &
         this%acetate_methanogen_c_unsat_col)
    call add_state('MM_ACET_METH_C_SAT', 'gC/m^3', &
         'acetoclastic methanogen biomass carbon in saturated subarea', &
         this%acetate_methanogen_c_sat_col)
    call add_state('MM_H2_METH_C_UNSAT', 'gC/m^3', &
         'hydrogenotrophic methanogen biomass carbon in unsaturated subarea', &
         this%h2_methanogen_c_unsat_col)
    call add_state('MM_H2_METH_C_SAT', 'gC/m^3', &
         'hydrogenotrophic methanogen biomass carbon in saturated subarea', &
         this%h2_methanogen_c_sat_col)
    call add_state('MM_AER_METHANOTROPH_C_UNSAT', 'gC/m^3', &
         'aerobic methanotroph biomass carbon in unsaturated subarea', &
         this%aerobic_methanotroph_c_unsat_col)
    call add_state('MM_AER_METHANOTROPH_C_SAT', 'gC/m^3', &
         'aerobic methanotroph biomass carbon in saturated subarea', &
         this%aerobic_methanotroph_c_sat_col)
    call add_state('MM_ANAER_METHANOTROPH_C_UNSAT', 'gC/m^3', &
         'anaerobic methanotroph biomass carbon in unsaturated subarea', &
         this%anaerobic_methanotroph_c_unsat_col)
    call add_state('MM_ANAER_METHANOTROPH_C_SAT', 'gC/m^3', &
         'anaerobic methanotroph biomass carbon in saturated subarea', &
         this%anaerobic_methanotroph_c_sat_col)
    call add_state('MM_CONC_CH4_UNSAT', 'mol/m^3', &
         'bulk-soil methane inventory density in unsaturated subarea', this%conc_ch4_unsat_col)
    call add_state('MM_CONC_CH4_SAT', 'mol/m^3', &
         'bulk-soil methane inventory density in saturated subarea', this%conc_ch4_sat_col)
    call add_state('MM_CONC_O2_UNSAT', 'mol/m^3', &
         'bulk-soil oxygen inventory density in unsaturated subarea', this%conc_o2_unsat_col)
    call add_state('MM_CONC_O2_SAT', 'mol/m^3', &
         'bulk-soil oxygen inventory density in saturated subarea', this%conc_o2_sat_col)
    call add_state('MM_CONC_CO2_UNSAT', 'mol/m^3', &
         'bulk-soil carbon dioxide inventory density in unsaturated subarea', this%conc_co2_unsat_col)
    call add_state('MM_CONC_CO2_SAT', 'mol/m^3', &
         'bulk-soil carbon dioxide inventory density in saturated subarea', this%conc_co2_sat_col)
    call add_state('MM_CONC_H2_UNSAT', 'mol/m^3', &
         'bulk-soil hydrogen inventory density in unsaturated subarea', this%conc_h2_unsat_col)
    call add_state('MM_CONC_H2_SAT', 'mol/m^3', &
         'bulk-soil hydrogen inventory density in saturated subarea', this%conc_h2_sat_col)
    call hist_addfld_decomp(fname='MM_O2_STRESS_UNSAT', units='1', type2d='levdcmp', &
         avgflag='A', &
         long_name='raw shared oxygen stress in unsaturated subarea before decomposition floor', &
         ptr_col=this%o2_stress_unsat_col, default='inactive')
    call hist_addfld_decomp(fname='MM_O2_STRESS_SAT', units='1', type2d='levdcmp', &
         avgflag='A', &
         long_name='raw shared oxygen stress in saturated subarea before decomposition floor', &
         ptr_col=this%o2_stress_sat_col, default='inactive')
    call hist_addfld1d(fname='MM_ADDITIONAL_C', units='gC/m^2', avgflag='A', &
         long_name='revised methane carbon storage excluding authoritative DOM', &
         ptr_col=this%additional_carbon_col, default='inactive')
    call hist_addfld1d(fname='MM_SURFACE_C_FLUX', units='gC/m^2/s', avgflag='A', &
         long_name='net revised methane CH4 plus CO2 carbon flux; positive to atmosphere', &
         ptr_col=this%surface_carbon_flux_col, default='inactive')
    call hist_addfld1d(fname='MM_LATERAL_C_FLUX', units='gC/m^2/s', avgflag='A', &
         long_name='net revised methane CH4 plus CO2 lateral carbon transfer; positive into column', &
         ptr_col=this%lateral_carbon_flux_col, default='inactive')
    call hist_addfld_decomp(fname='MM_DOM_ADV_FLUX', units='gC/m^2/s', type2d='levdcmp', &
         avgflag='A', long_name='downward DOM carbon advective flux at lower layer face', &
         ptr_col=this%dom_advective_flux_col, default='inactive')
    call hist_addfld_decomp(fname='MM_DOM_DIFF_FLUX', units='gC/m^2/s', type2d='levdcmp', &
         avgflag='A', long_name='downward DOM carbon diffusive flux at lower layer face', &
         ptr_col=this%dom_diffusive_flux_col, default='inactive')
    call add_aqueous_flux('MM_AQUEOUS_C_EXPORT', 'gC/m^2/s', &
         'DOM plus acetate aqueous carbon export', this%aqueous_carbon_export_col)
    call add_aqueous_flux('MM_AQUEOUS_N_EXPORT', 'gN/m^2/s', &
         'DOM aqueous nitrogen export', this%aqueous_nitrogen_export_col)
    call add_aqueous_flux('MM_AQUEOUS_P_EXPORT', 'gP/m^2/s', &
         'DOM aqueous phosphorus export', this%aqueous_phosphorus_export_col)
    call add_aqueous_flux('MM_DOM_BOTTOM_C_EXPORT', 'gC/m^2/s', &
         'DOM carbon export through the bottom boundary', this%dom_bottom_carbon_export_col)
    call add_aqueous_flux('MM_DOM_BOTTOM_N_EXPORT', 'gN/m^2/s', &
         'DOM nitrogen export through the bottom boundary', this%dom_bottom_nitrogen_export_col)
    call add_aqueous_flux('MM_DOM_BOTTOM_P_EXPORT', 'gP/m^2/s', &
         'DOM phosphorus export through the bottom boundary', this%dom_bottom_phosphorus_export_col)
    call add_aqueous_flux('MM_DOM_C_TRANSPORT_RESIDUAL', 'gC/m^2', &
         'DOM carbon aqueous-transport inventory residual per timestep', &
         this%dom_carbon_transport_residual_col)
    call add_aqueous_flux('MM_DOM_N_TRANSPORT_RESIDUAL', 'gN/m^2', &
         'DOM nitrogen aqueous-transport inventory residual per timestep', &
         this%dom_nitrogen_transport_residual_col)
    call add_aqueous_flux('MM_DOM_P_TRANSPORT_RESIDUAL', 'gP/m^2', &
         'DOM phosphorus aqueous-transport inventory residual per timestep', &
         this%dom_phosphorus_transport_residual_col)
    call hist_addfld1d(fname='MM_SURFACE_CH4_FLUX', units='kgC/m^2/s', avgflag='A', &
         long_name='net revised methane CH4 carbon flux; positive to atmosphere', &
         ptr_col=this%surface_ch4_flux_col, default='inactive')
    call hist_addfld1d(fname='MM_SURFACE_CO2_FLUX', units='gC/m^2/s', avgflag='A', &
         long_name='net revised methane CO2 carbon flux; positive to atmosphere', &
         ptr_col=this%surface_co2_flux_col, default='inactive')
    call add_ch4_surface_path('MM_CH4_SURF_DIFF', &
         'surface CH4 diffusion', this%ch4_surface_diffusion_col)
    call add_ch4_surface_path('MM_CH4_SURF_DIFF_UNSAT', &
         'area-weighted unsaturated-subarea surface CH4 diffusion', &
         this%ch4_surface_diffusion_unsat_col)
    call add_ch4_surface_path('MM_CH4_SURF_DIFF_SAT', &
         'area-weighted saturated-subarea surface CH4 diffusion', &
         this%ch4_surface_diffusion_sat_col)
    call add_ch4_surface_path('MM_CH4_SURF_AERE', &
         'surface CH4 aerenchyma transport', this%ch4_surface_aerenchyma_col)
    call add_ch4_surface_path('MM_CH4_SURF_AERE_UNSAT', &
         'area-weighted unsaturated-subarea surface CH4 aerenchyma transport', &
         this%ch4_surface_aerenchyma_unsat_col)
    call add_ch4_surface_path('MM_CH4_SURF_AERE_SAT', &
         'area-weighted saturated-subarea surface CH4 aerenchyma transport', &
         this%ch4_surface_aerenchyma_sat_col)
    call add_ch4_surface_path('MM_CH4_SURF_EBUL', &
         'surface CH4 ebullition', this%ch4_surface_ebullition_col)
    call add_ch4_surface_path('MM_CH4_SURF_EBUL_UNSAT', &
         'area-weighted unsaturated-subarea surface CH4 ebullition', &
         this%ch4_surface_ebullition_unsat_col)
    call add_ch4_surface_path('MM_CH4_SURF_EBUL_SAT', &
         'area-weighted saturated-subarea surface CH4 ebullition', &
         this%ch4_surface_ebullition_sat_col)
    call hist_addfld1d(fname='MM_CH4_PROD', units='gC/m^2/s', avgflag='A', &
         long_name='gross revised methane production integrated over the soil column', &
         ptr_col=this%ch4_production_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_OXID', units='gC/m^2/s', avgflag='A', &
         long_name='gross revised methane oxidation integrated over the soil column', &
         ptr_col=this%ch4_oxidation_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_PROD_UNSAT', units='gC/m^2/s', avgflag='A', &
         long_name='area-weighted gross methane production from unsaturated subarea', &
         ptr_col=this%ch4_production_unsat_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_PROD_SAT', units='gC/m^2/s', avgflag='A', &
         long_name='area-weighted gross methane production from saturated subarea', &
         ptr_col=this%ch4_production_sat_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_OXID_UNSAT', units='gC/m^2/s', avgflag='A', &
         long_name='area-weighted gross methane oxidation from unsaturated subarea', &
         ptr_col=this%ch4_oxidation_unsat_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_OXID_SAT', units='gC/m^2/s', avgflag='A', &
         long_name='area-weighted gross methane oxidation from saturated subarea', &
         ptr_col=this%ch4_oxidation_sat_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_OXID_AER', units='gC/m^2/s', avgflag='A', &
         long_name='gross aerobic methane oxidation integrated over the soil column', &
         ptr_col=this%ch4_aerobic_oxidation_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_OXID_AER_UNSAT', units='gC/m^2/s', avgflag='A', &
         long_name='area-weighted aerobic methane oxidation from unsaturated subarea', &
         ptr_col=this%ch4_aerobic_oxidation_unsat_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_OXID_AER_SAT', units='gC/m^2/s', avgflag='A', &
         long_name='area-weighted aerobic methane oxidation from saturated subarea', &
         ptr_col=this%ch4_aerobic_oxidation_sat_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_OXID_AER_PREO2', units='gC/m^2/s', avgflag='A', &
         long_name='aerobic methane oxidation after CH4 limitation and before shared O2 limitation', &
         ptr_col=this%ch4_aerobic_oxidation_pre_o2_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_OXID_AER_PREO2_UNSAT', units='gC/m^2/s', avgflag='A', &
         long_name='area-weighted pre-O2-limit aerobic methane oxidation from unsaturated subarea', &
         ptr_col=this%ch4_aerobic_oxidation_pre_o2_unsat_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_OXID_AER_PREO2_SAT', units='gC/m^2/s', avgflag='A', &
         long_name='area-weighted pre-O2-limit aerobic methane oxidation from saturated subarea', &
         ptr_col=this%ch4_aerobic_oxidation_pre_o2_sat_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_OXID_AOM', units='gC/m^2/s', avgflag='A', &
         long_name='gross anaerobic methane oxidation integrated over the soil column', &
         ptr_col=this%ch4_anaerobic_oxidation_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_OXID_AOM_UNSAT', units='gC/m^2/s', avgflag='A', &
         long_name='area-weighted anaerobic methane oxidation from unsaturated subarea', &
         ptr_col=this%ch4_anaerobic_oxidation_unsat_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_OXID_AOM_SAT', units='gC/m^2/s', avgflag='A', &
         long_name='area-weighted anaerobic methane oxidation from saturated subarea', &
         ptr_col=this%ch4_anaerobic_oxidation_sat_col, default='inactive')
    call hist_addfld1d(fname='MM_SAT_FRACTION', units='1', avgflag='A', &
         long_name='revised methane saturated-area fraction', &
         ptr_col=this%sat_fraction_previous_col, default='inactive')

  contains
    subroutine add_state(name, units, long_name, field)
      character(len=*), intent(in) :: name, units, long_name
      real(r8), pointer, intent(inout) :: field(:,:)
      call hist_addfld_decomp(fname=name, units=units, type2d='levdcmp', &
           avgflag='A', long_name=long_name, ptr_col=field, default='inactive')
    end subroutine add_state

    subroutine add_ch4_surface_path(name, long_name, field)
      character(len=*), intent(in) :: name, long_name
      real(r8), pointer, intent(in) :: field(:)
      call hist_addfld1d(fname=name, units='gC/m^2/s', avgflag='A', &
           long_name=trim(long_name)//'; positive to atmosphere', &
           ptr_col=field, default='inactive')
    end subroutine add_ch4_surface_path

    subroutine add_aqueous_flux(name, units, long_name, field)
      character(len=*), intent(in) :: name, units, long_name
      real(r8), pointer, intent(inout) :: field(:)
      call hist_addfld1d(fname=name, units=units, avgflag='A', &
           long_name=long_name, ptr_col=field, default='inactive')
    end subroutine add_aqueous_flux
  end subroutine InitHistory

  subroutine Repartition(this, bounds, num_soilc, filter_soilc, saturated_fraction)
    use elm_varpar, only : nlevdecomp
    use MicrobeGasTransportMod, only : repartitionMicrobeMethaneScalar

    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc
    integer, intent(in) :: filter_soilc(:)
    real(r8), intent(in) :: saturated_fraction(bounds%begc:bounds%endc)
    integer :: c, fc, j
    real(r8) :: old_fraction, new_fraction

    if (.not. use_microbe_methane .or. use_legacy_ch4_with_microbe) return

    do fc = 1, num_soilc
       c = filter_soilc(fc)
       old_fraction = this%sat_fraction_previous_col(c)
       new_fraction = min(1._r8, max(0._r8, saturated_fraction(c)))
       do j = 1, nlevdecomp
          call repartition_pair(this%acetate_c_unsat_col(c,j), this%acetate_c_sat_col(c,j))
          call repartition_pair(this%acetate_methanogen_c_unsat_col(c,j), &
               this%acetate_methanogen_c_sat_col(c,j))
          call repartition_pair(this%h2_methanogen_c_unsat_col(c,j), &
               this%h2_methanogen_c_sat_col(c,j))
          call repartition_pair(this%aerobic_methanotroph_c_unsat_col(c,j), &
               this%aerobic_methanotroph_c_sat_col(c,j))
          call repartition_pair(this%anaerobic_methanotroph_c_unsat_col(c,j), &
               this%anaerobic_methanotroph_c_sat_col(c,j))
          call repartition_pair(this%conc_ch4_unsat_col(c,j), this%conc_ch4_sat_col(c,j))
          call repartition_pair(this%conc_o2_unsat_col(c,j), this%conc_o2_sat_col(c,j))
          call repartition_pair(this%conc_co2_unsat_col(c,j), this%conc_co2_sat_col(c,j))
          call repartition_pair(this%conc_h2_unsat_col(c,j), this%conc_h2_sat_col(c,j))
       end do
       this%sat_fraction_previous_col(c) = new_fraction
    end do

  contains
    subroutine repartition_pair(unsaturated_concentration, saturated_concentration)
      real(r8), intent(inout) :: unsaturated_concentration, saturated_concentration
      real(r8) :: repartitioned_unsaturated, repartitioned_saturated

      call repartitionMicrobeMethaneScalar(old_fraction, new_fraction, &
           unsaturated_concentration, saturated_concentration, &
           repartitioned_unsaturated, repartitioned_saturated)
      unsaturated_concentration = repartitioned_unsaturated
      saturated_concentration = repartitioned_saturated
    end subroutine repartition_pair
  end subroutine Repartition

  subroutine Advance(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       dt, atm2lnd_vars, &
       col_es, col_ws, col_wf, chemstate_vars, soilstate_vars, soilhydrology_vars, &
       ground_conductance_patch, ch4_vars, col_cs, col_cf, veg_cf, col_ns, col_nf, col_ps)
    ! Mutable ELM boundary for the pure revised-methane kernels. Every update
    ! for one column is staged locally, checked, and then committed once.
    use abortutils, only : endrun
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use elm_varctl, only : iulog
    use elm_varpar, only : nlevdecomp, i_dom
    use elm_varcon, only : denh2o, denice, tfrz, d_con_w, d_con_g, catomw, rgas
    use elm_varcon, only : c_h_inv, kh_theta, kh_tbase
    use ColumnType, only : col_pp
    use TopounitType, only : top_pp
    use ColumnDataType, only : column_energy_state, column_water_state, column_water_flux
    use ColumnDataType, only : column_carbon_state, column_carbon_flux
    use ColumnDataType, only : column_nitrogen_state, column_nitrogen_flux
    use ColumnDataType, only : column_phosphorus_state
    use VegetationType, only : veg_pp
    use VegetationDataType, only : vegetation_carbon_flux
    use pftvarcon, only : noveg
    use atm2lndType, only : atm2lnd_type
    use ChemStateType, only : chemstate_type
    use SoilStateType, only : soilstate_type
    use SoilHydrologyType, only : soilhydrology_type
    use CH4Mod, only : ch4_type, CH4ParamsInst
    use SharedParamsMod, only : ParamsShareInst
    use subgridAveMod, only : p2c
    use MicrobeDecompMod, only : MicrobeDecompParamsInst
    use MicrobeMethaneParamsMod, only : MicrobeMethaneParamsInst
    use MicrobeMethaneReactionMod, only : microbe_methane_reaction_state_type
    use MicrobeMethaneReactionMod, only : microbe_methane_reaction_environment_type
    use MicrobeMethaneStateUpdateMod, only : microbe_methane_reaction_transaction_type
    use MicrobeMethaneStateUpdateMod, only : advanceMicrobeMethaneReactionLayer
    use MicrobeMethaneStateUpdateMod, only : advanceMicrobeMethaneAcetateTransport
    use MicrobeMethaneStateUpdateMod, only : advanceMicrobeMethaneDOMRelaxation
    use MicrobeMethaneStateUpdateMod, only : advanceMicrobeAqueousTracerTransport
    use MicrobeMethaneStateUpdateMod, only : advanceMicrobeMethaneGasTransport
    use MicrobeMethaneStateUpdateMod, only : microbeMethaneColumnAdditionalCarbon
    use MicrobeMethaneStateUpdateMod, only : microbeMethaneSurfaceCarbonFlux
    use MicrobeMethaneStateUpdateMod, only : microbeMethaneCH4SurfaceFluxKgC
    use MicrobeMethaneStateUpdateMod, only : microbeMethaneCO2Correction
    use MicrobeGasTransportMod, only : repartitionMicrobeMethaneScalar
    use MicrobeGasTransportMod, only : microbeMethaneEffectiveAqueousDiffusivity
    use MicrobeGasTransportMod, only : computeMicrobeTopounitLateralDiffusion
    use MicrobeGasTransportMod, only : microbe_gas_ch4, microbe_gas_o2
    use MicrobeGasTransportMod, only : microbe_gas_co2, microbe_gas_h2
    use MicrobeGasTransportMod, only : microbe_gas_count
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc
    integer, intent(in) :: filter_soilc(:)
    integer, intent(in) :: num_soilp
    integer, intent(in) :: filter_soilp(:)
    real(r8), intent(in) :: dt
    type(atm2lnd_type), intent(in) :: atm2lnd_vars
    type(column_energy_state), intent(in) :: col_es
    type(column_water_state), intent(in) :: col_ws
    type(column_water_flux), intent(in) :: col_wf
    type(chemstate_type), intent(in) :: chemstate_vars
    type(soilstate_type), intent(in) :: soilstate_vars
    type(soilhydrology_type), intent(in) :: soilhydrology_vars
    real(r8), intent(in) :: ground_conductance_patch(bounds%begp:)
    type(ch4_type), intent(inout) :: ch4_vars
    type(column_carbon_state), intent(inout) :: col_cs
    type(column_carbon_flux), intent(in) :: col_cf
    type(vegetation_carbon_flux), intent(in) :: veg_cf
    type(column_nitrogen_state), intent(inout) :: col_ns
    type(column_nitrogen_flux), intent(in) :: col_nf
    type(column_phosphorus_state), intent(inout) :: col_ps

    type(microbe_methane_reaction_state_type) :: unsaturated_state(nlevdecomp)
    type(microbe_methane_reaction_state_type) :: saturated_state(nlevdecomp)
    type(microbe_methane_reaction_state_type) :: unsaturated_work(nlevdecomp)
    type(microbe_methane_reaction_state_type) :: saturated_work(nlevdecomp)
    type(microbe_methane_reaction_state_type) :: unsaturated_candidate(nlevdecomp)
    type(microbe_methane_reaction_state_type) :: saturated_candidate(nlevdecomp)
    type(microbe_methane_reaction_environment_type) :: unsaturated_environment
    type(microbe_methane_reaction_environment_type) :: saturated_environment
    type(microbe_methane_reaction_transaction_type) :: reaction(nlevdecomp)
    real(r8) :: dom_c(nlevdecomp), dom_n(nlevdecomp), dom_p(nlevdecomp)
    real(r8) :: updated_dom_c(nlevdecomp), updated_dom_n(nlevdecomp)
    real(r8) :: updated_dom_p(nlevdecomp), dom_relaxation_rate(nlevdecomp)
    real(r8) :: liquid_fraction(nlevdecomp), porosity(nlevdecomp)
    real(r8) :: water_flux(0:nlevdecomp)
    real(r8) :: dom_diffusion_conductivity(nlevdecomp)
    real(r8) :: acetate_diffusion_conductivity(nlevdecomp)
    real(r8) :: solute_advective_flux(0:nlevdecomp)
    real(r8) :: solute_diffusive_flux(0:nlevdecomp)
    real(r8) :: solute_tendency(nlevdecomp)
    real(r8) :: acetate_concentration(nlevdecomp), updated_acetate(nlevdecomp)
    real(r8) :: mineral_n(nlevdecomp), mineral_p(nlevdecomp)
    real(r8) :: layer_thickness(nlevdecomp), layer_depth(nlevdecomp)
    real(r8) :: unsaturated_diffusivity(nlevdecomp,microbe_gas_count)
    real(r8) :: saturated_diffusivity(nlevdecomp,microbe_gas_count)
    real(r8) :: unsaturated_transport_capacity(nlevdecomp,microbe_gas_count)
    real(r8) :: saturated_transport_capacity(nlevdecomp,microbe_gas_count)
    real(r8) :: unsaturated_acetate_diffusivity(nlevdecomp)
    real(r8) :: saturated_acetate_diffusivity(nlevdecomp)
    real(r8) :: surface_equilibrium(microbe_gas_count)
    real(r8) :: aerenchyma_equilibrium(nlevdecomp,microbe_gas_count)
    real(r8) :: aerenchyma_exchange_rate(nlevdecomp,microbe_gas_count)
    real(r8) :: aerenchyma_minimum_emission(nlevdecomp,microbe_gas_count)
    real(r8) :: unsaturated_surface_conductance(microbe_gas_count)
    real(r8) :: saturated_surface_conductance(microbe_gas_count)
    real(r8) :: ch4_ebullition_threshold(nlevdecomp)
    real(r8) :: unsaturated_ebullition_activation(nlevdecomp)
    real(r8) :: saturated_ebullition_activation(nlevdecomp)
    real(r8) :: interface_flux(0:nlevdecomp,microbe_gas_count)
    real(r8) :: aerenchyma_flux(nlevdecomp,microbe_gas_count)
    real(r8) :: ch4_ebullition_loss(nlevdecomp)
    real(r8) :: unsaturated_surface_diffusive_flux(microbe_gas_count)
    real(r8) :: unsaturated_surface_aerenchyma_flux(microbe_gas_count)
    real(r8) :: unsaturated_surface_flux(microbe_gas_count)
    real(r8) :: saturated_surface_flux(microbe_gas_count)
    real(r8) :: unsaturated_surface_ebullition_flux
    real(r8) :: acetate_interface_flux(0:nlevdecomp)
    real(r8) :: acetate_tendency(nlevdecomp)
    real(r8) :: ground_conductance_col(bounds%begc:bounds%endc)
    real(r8) :: root_fraction_col(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: root_respiration_col_vr(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: partial_pressure(microbe_gas_count)
    real(r8) :: liquid_saturation, thawed_fraction, moisture_scalar
    real(r8) :: saturation_scalar, fraction, old_fraction, depth_scale
    real(r8) :: bulk_surface_flux(microbe_gas_count)
    real(r8) :: saturated_surface_diffusive_flux(microbe_gas_count)
    real(r8) :: saturated_surface_aerenchyma_flux(microbe_gas_count)
    real(r8) :: saturated_surface_ebullition_flux
    real(r8) :: unsaturated_gas_residual, saturated_gas_residual
    real(r8) :: unsaturated_acetate_residual, saturated_acetate_residual
    real(r8) :: dom_carbon_residual, dom_nitrogen_residual, dom_phosphorus_residual
    real(r8) :: dom_carbon_export, dom_nitrogen_export, dom_phosphorus_export
    real(r8) :: unsaturated_acetate_export, saturated_acetate_export
    real(r8) :: aqueous_temperature_scalar, relative_liquid_saturation
    logical :: column_valid
    logical :: unsaturated_acetate_valid, saturated_acetate_valid
    logical :: unsaturated_gas_valid, saturated_gas_valid
    logical :: dom_relaxation_valid
    logical :: dom_carbon_transport_valid, dom_nitrogen_transport_valid
    logical :: dom_phosphorus_transport_valid
    logical :: aerenchyma_allows_influx(microbe_gas_count)
    integer :: c, fc, fp, g, gas, j, p
    character(len=512) :: message

    if (.not. use_microbe_methane .or. use_legacy_ch4_with_microbe) return
    if (dt <= 0._r8) then
       call endrun(msg=' ERROR: revised methane adapter requires a positive timestep'//&
            errMsg(__FILE__, __LINE__))
    end if

    this%surface_carbon_flux_col(bounds%begc:bounds%endc) = 0._r8
    this%lateral_carbon_flux_col(bounds%begc:bounds%endc) = 0._r8
    this%dom_advective_flux_col(bounds%begc:bounds%endc,:) = 0._r8
    this%dom_diffusive_flux_col(bounds%begc:bounds%endc,:) = 0._r8
    this%aqueous_carbon_export_col(bounds%begc:bounds%endc) = 0._r8
    this%aqueous_nitrogen_export_col(bounds%begc:bounds%endc) = 0._r8
    this%aqueous_phosphorus_export_col(bounds%begc:bounds%endc) = 0._r8
    this%dom_bottom_carbon_export_col(bounds%begc:bounds%endc) = 0._r8
    this%dom_bottom_nitrogen_export_col(bounds%begc:bounds%endc) = 0._r8
    this%dom_bottom_phosphorus_export_col(bounds%begc:bounds%endc) = 0._r8
    this%dom_carbon_transport_residual_col(bounds%begc:bounds%endc) = 0._r8
    this%dom_nitrogen_transport_residual_col(bounds%begc:bounds%endc) = 0._r8
    this%dom_phosphorus_transport_residual_col(bounds%begc:bounds%endc) = 0._r8
    this%surface_ch4_flux_col(bounds%begc:bounds%endc) = 0._r8
    this%surface_co2_flux_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_surface_diffusion_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_surface_diffusion_unsat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_surface_diffusion_sat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_surface_aerenchyma_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_surface_aerenchyma_unsat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_surface_aerenchyma_sat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_surface_ebullition_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_surface_ebullition_unsat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_surface_ebullition_sat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_production_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_oxidation_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_production_unsat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_production_sat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_oxidation_unsat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_oxidation_sat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_aerobic_oxidation_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_aerobic_oxidation_unsat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_aerobic_oxidation_sat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_aerobic_oxidation_pre_o2_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_aerobic_oxidation_pre_o2_unsat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_aerobic_oxidation_pre_o2_sat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_anaerobic_oxidation_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_anaerobic_oxidation_unsat_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_anaerobic_oxidation_sat_col(bounds%begc:bounds%endc) = 0._r8

    call p2c(bounds, num_soilc, filter_soilc, &
         ground_conductance_patch(bounds%begp:bounds%endp), &
         ground_conductance_col(bounds%begc:bounds%endc))
    ! Legacy CH4 computes this column average inside its solver. Revised mode
    ! bypasses that solver, so it must establish its own plant-transport
    ! input directly from the authoritative patch root profile.
    call p2c(bounds, nlevdecomp, &
         soilstate_vars%rootfr_patch(bounds%begp:bounds%endp,1:nlevdecomp), &
         root_fraction_col(bounds%begc:bounds%endc,1:nlevdecomp), 0)

    ! CH4Mod normally constructs this profile inside ch4_prod. Revised mode
    ! bypasses that routine, so col_cf%rr_vr remains at its initialization
    ! fill value. Reproduce the legacy patch-to-column weighting here from
    ! the authoritative root respiration flux and root profile.
    root_respiration_col_vr = 0._r8
    do fp = 1, num_soilp
       p = filter_soilp(fp)
       c = veg_pp%column(p)
       if (.not. col_pp%is_fates(c)) then
          if (veg_pp%wtcol(p) > 0._r8 .and. veg_pp%itype(p) /= noveg) then
             do j = 1, nlevdecomp
                root_respiration_col_vr(c,j) = root_respiration_col_vr(c,j) + &
                     veg_cf%rr(p) * soilstate_vars%rootfr_patch(p,j) * veg_pp%wtcol(p)
             end do
          end if
       end if
    end do

    do fc = 1, num_soilc
       c = filter_soilc(fc)
       g = col_pp%gridcell(c)
       fraction = clampUnitInterval(max(soilhydrology_vars%fsat_col(c), &
            col_ws%frac_h2osfc(c)))
       if (use_clm_microbe_humhol_saturation) fraction = 0.99_r8
       old_fraction = clampUnitInterval(this%sat_fraction_previous_col(c))
       ! ELM hydrology reports interface water flux in mm s-1, positive
       ! downward. Convert once to the m s-1 required by solute transport.
       ! The top value is the net soil-water boundary and can be negative
       ! when ground evaporation exceeds infiltration. Dissolved organic
       ! matter must not leave with water vapor, so only the downward clean-
       ! water infiltration component is admitted here. A future surface-
       ! water/runoff solute pathway must use an explicit aqueous surface pool.
       water_flux = 1.e-3_r8 * col_wf%qflx_adv(c,0:nlevdecomp)
       water_flux(0) = max(0._r8, water_flux(0))

       partial_pressure(microbe_gas_ch4) = atmosphericPartialPressure( &
            atm2lnd_vars%forc_pch4_grc(g), atm2lnd_vars%forc_pbot_downscaled_col(c), &
            MicrobeMethaneParamsInst%atmospheric_ch4_mixing_ratio)
       partial_pressure(microbe_gas_o2) = atmosphericPartialPressure( &
            atm2lnd_vars%forc_po2_grc(g), atm2lnd_vars%forc_pbot_downscaled_col(c), &
            MicrobeMethaneParamsInst%atmospheric_o2_mixing_ratio)
       partial_pressure(microbe_gas_co2) = atmosphericPartialPressure( &
            atm2lnd_vars%forc_pco2_grc(g), atm2lnd_vars%forc_pbot_downscaled_col(c), &
            MicrobeMethaneParamsInst%atmospheric_co2_mixing_ratio)
       partial_pressure(microbe_gas_h2) = atmosphericPartialPressure(0._r8, &
            atm2lnd_vars%forc_pbot_downscaled_col(c), &
            MicrobeMethaneParamsInst%atmospheric_h2_mixing_ratio)
       if (any(.not. ieee_is_finite(partial_pressure)) .or. &
            any(partial_pressure < 0._r8)) then
          write(message,'(a,i0)') ' ERROR: invalid revised methane atmospheric boundary for column ', c
          call endrun(msg=trim(message)//errMsg(__FILE__, __LINE__))
       end if

       ! CH4 and H2 plant transport are one-way emission pathways. Their
       ! thresholds must not replace the atmospheric Henry-law equilibrium or
       ! create an artificial atmospheric source when soil concentrations are low.
       aerenchyma_minimum_emission = 0._r8
       aerenchyma_minimum_emission(:,microbe_gas_ch4) = &
            1.e-3_r8 * MicrobeMethaneParamsInst%ch4_transport_threshold
       aerenchyma_minimum_emission(:,microbe_gas_h2) = &
            1.e-3_r8 * MicrobeMethaneParamsInst%h2_plant_transport_threshold
       aerenchyma_allows_influx = .true.
       aerenchyma_allows_influx(microbe_gas_ch4) = .false.
       aerenchyma_allows_influx(microbe_gas_h2) = .false.

       call gatherColumnState(c, unsaturated_state, saturated_state)
       call repartitionColumnState(old_fraction, fraction, unsaturated_state, saturated_state)
       unsaturated_work = unsaturated_state
       saturated_work = saturated_state
       column_valid = .true.

       do j = 1, nlevdecomp
          layer_thickness(j) = col_pp%dz(c,j)
          layer_depth(j) = col_pp%z(c,j)
          dom_c(j) = col_cs%decomp_cpools_vr(c,j,i_dom)
          dom_n(j) = col_ns%decomp_npools_vr(c,j,i_dom)
          dom_p(j) = col_ps%decomp_ppools_vr(c,j,i_dom)
          mineral_n(j) = col_ns%smin_nh4_vr(c,j)
          mineral_p(j) = col_ps%solutionp_vr(c,j)

          liquid_saturation = layerLiquidSaturation(c, j)
          thawed_fraction = layerThawedFraction(c, j)
          liquid_fraction(j) = max(0._r8, col_ws%h2osoi_liq(c,j)) / &
               (denh2o * max(layer_thickness(j), tiny(1._r8)))
          porosity(j) = max(soilstate_vars%watsat_col(c,j), tiny(1._r8))
          relative_liquid_saturation = clampUnitInterval(liquid_fraction(j) / porosity(j))
          aqueous_temperature_scalar = &
               (max(col_es%t_soisno(c,j), tiny(1._r8)) / &
               MicrobeMethaneParamsInst%aqueous_diffusion_t_ref) ** &
               MicrobeMethaneParamsInst%aqueous_diffusion_temperature_exponent
          ! These are theta*D_eff [m2 s-1], the conductivity multiplying
          ! a porewater concentration gradient. Mechanical dispersion uses
          ! alpha*abs(q); both terms shut down as liquid water freezes.
          dom_diffusion_conductivity(j) = thawed_fraction * &
               (liquid_fraction(j) * &
               MicrobeMethaneParamsInst%aqueous_dom_molecular_diffusivity * &
               relative_liquid_saturation ** &
               MicrobeMethaneParamsInst%aqueous_solute_tortuosity_exponent * &
               aqueous_temperature_scalar + &
               MicrobeMethaneParamsInst%aqueous_solute_dispersivity * &
               0.5_r8 * (abs(water_flux(j-1)) + abs(water_flux(j))))
          acetate_diffusion_conductivity(j) = thawed_fraction * &
               (liquid_fraction(j) * &
               MicrobeMethaneParamsInst%aqueous_acetate_molecular_diffusivity * &
               relative_liquid_saturation ** &
               MicrobeMethaneParamsInst%aqueous_solute_tortuosity_exponent * &
               aqueous_temperature_scalar + &
               MicrobeMethaneParamsInst%aqueous_solute_dispersivity * &
               0.5_r8 * (abs(water_flux(j-1)) + abs(water_flux(j))))
          moisture_scalar = soilMoistureResponse(soilstate_vars%soilpsi_col(c,j), &
               soilstate_vars%sucsat_col(c,j), &
               MicrobeMethaneParamsInst%soil_water_potential_min)
          saturation_scalar = clampUnitInterval(liquid_saturation / &
               max(MicrobeMethaneParamsInst%saturation_reaction_threshold, tiny(1._r8)))

          unsaturated_environment%soil_temperature = col_es%t_soisno(c,j)
          ! Native ELM allocates chemstate soil pH but does not currently
          ! populate it. External chemistry backends that do populate it are
          ! rejected by the revised-methane configuration gate. Use the named
          ! named parameter-file fallback until Phase 4 adds a spatial
          ! soil-pH input; acetate feedback can still lower the effective pH
          ! inside the reaction kernel.
          unsaturated_environment%soil_ph = MicrobeMethaneParamsInst%soil_ph_fallback
          unsaturated_environment%dom_fermentation_scalar = &
               moisture_scalar * saturation_scalar
          unsaturated_environment%aerobic_acetate_oxidation_scalar = &
               moisture_scalar * (1._r8 - saturation_scalar)
          unsaturated_environment%elm_aerobic_o2_demand = &
               elmAerobicOxygenDemand(c, j, root_respiration_col_vr(c,j))
          saturated_environment%soil_temperature = col_es%t_soisno(c,j)
          saturated_environment%soil_ph = MicrobeMethaneParamsInst%soil_ph_fallback
          saturated_environment%dom_fermentation_scalar = moisture_scalar
          saturated_environment%aerobic_acetate_oxidation_scalar = 0._r8
          saturated_environment%elm_aerobic_o2_demand = &
               unsaturated_environment%elm_aerobic_o2_demand

          call advanceMicrobeMethaneReactionLayer(dom_c(j), dom_n(j), dom_p(j), &
               mineral_n(j), mineral_p(j), fraction, unsaturated_work(j), &
               saturated_work(j), unsaturated_environment, saturated_environment, &
               MicrobeMethaneParamsInst, MicrobeDecompParamsInst%cn_dom, &
               MicrobeDecompParamsInst%cp_dom, dt, reaction(j))
          if (.not. reaction(j)%valid) column_valid = .false.
          this%o2_stress_unsat_col(c,j) = reaction(j)%unsaturated_rates%oxygen_stress
          this%o2_stress_sat_col(c,j) = reaction(j)%saturated_rates%oxygen_stress
          dom_c(j) = reaction(j)%dom_c
          dom_n(j) = reaction(j)%dom_n
          dom_p(j) = reaction(j)%dom_p
          mineral_n(j) = reaction(j)%mineral_n
          mineral_p(j) = reaction(j)%mineral_p
          unsaturated_work(j) = reaction(j)%unsaturated_state
          saturated_work(j) = reaction(j)%saturated_state

          do gas = 1, microbe_gas_count
             call gasTransportProperties(c, j, gas, .false., &
                  unsaturated_transport_capacity(j,gas), &
                  unsaturated_diffusivity(j,gas))
             call gasTransportProperties(c, j, gas, .true., &
                  saturated_transport_capacity(j,gas), &
                  saturated_diffusivity(j,gas))
             aerenchyma_equilibrium(j,gas) = henryEquilibriumConcentration( &
                  partial_pressure(gas), col_es%t_soisno(c,j), gas)
          end do
          if (use_elm_microbe_methane_transport) then
             unsaturated_acetate_diffusivity(j) = &
                  MicrobeMethaneParamsInst%dom_diffusivity * liquid_saturation * thawed_fraction
             saturated_acetate_diffusivity(j) = &
                  MicrobeMethaneParamsInst%dom_diffusivity * thawed_fraction
          else
             ! CLM-Microbe applies the same aqueous Fickian coefficient to
             ! both area partitions and scales it only by temperature.
             unsaturated_acetate_diffusivity(j) = legacyFickianAcetateDiffusivity( &
                  col_es%t_soisno(c,j))
             saturated_acetate_diffusivity(j) = unsaturated_acetate_diffusivity(j)
          end if

          ! The legacy coefficient is dimensionally m s-1: division by a
          ! finite layer depth converts it to the s-1 rate required by the
          ! plant exchange kernel.
          depth_scale = max(layer_depth(j), 0.5_r8 * layer_thickness(j), tiny(1._r8))
          aerenchyma_exchange_rate(j,:) = MicrobeMethaneParamsInst%plant_transport_coefficient * &
               validRootFraction(root_fraction_col(c,j)) / depth_scale
          aerenchyma_exchange_rate(j,microbe_gas_ch4) = &
               aerenchyma_exchange_rate(j,microbe_gas_ch4) * &
               exp(-layer_depth(j) / MicrobeMethaneParamsInst%ch4_h2_root_efold_depth)
          aerenchyma_exchange_rate(j,microbe_gas_h2) = &
               aerenchyma_exchange_rate(j,microbe_gas_h2) * &
               exp(-layer_depth(j) / MicrobeMethaneParamsInst%ch4_h2_root_efold_depth)
          aerenchyma_exchange_rate(j,microbe_gas_o2) = &
               aerenchyma_exchange_rate(j,microbe_gas_o2) * &
               exp(-layer_depth(j) / MicrobeMethaneParamsInst%o2_root_efold_depth) * &
               (1._r8 - MicrobeMethaneParamsInst%plant_o2_consumption_fraction)
          aerenchyma_exchange_rate(j,microbe_gas_co2) = &
               aerenchyma_exchange_rate(j,microbe_gas_co2) * &
               exp(-layer_depth(j) / MicrobeMethaneParamsInst%ch4_h2_root_efold_depth) * &
               MicrobeMethaneParamsInst%plant_co2_flux_fraction
          if (col_es%t_soisno(c,j) < MicrobeMethaneParamsInst%transport_thaw_threshold) then
             aerenchyma_exchange_rate(j,:) = 0._r8
          end if
          ch4_ebullition_threshold(j) = &
               1.e-3_r8 * MicrobeMethaneParamsInst%ch4_transport_threshold
          unsaturated_ebullition_activation(j) = thawed_fraction * liquid_saturation * &
               exp(-layer_depth(j) / MicrobeMethaneParamsInst%ebullition_efold_depth)
          saturated_ebullition_activation(j) = thawed_fraction * &
               exp(-layer_depth(j) / MicrobeMethaneParamsInst%ebullition_efold_depth)
          dom_relaxation_rate(j) = MicrobeMethaneParamsInst%dom_relaxation_rate * &
               (max(col_es%t_soisno(c,j), tiny(1._r8)) / &
               MicrobeMethaneParamsInst%aqueous_diffusion_t_ref) ** &
               MicrobeMethaneParamsInst%aqueous_diffusion_temperature_exponent
       end do

       dom_relaxation_valid = .true.
       dom_carbon_transport_valid = .true.
       dom_nitrogen_transport_valid = .true.
       dom_phosphorus_transport_valid = .true.
       dom_carbon_residual = 0._r8
       dom_nitrogen_residual = 0._r8
       dom_phosphorus_residual = 0._r8
       dom_carbon_export = 0._r8
       dom_nitrogen_export = 0._r8
       dom_phosphorus_export = 0._r8
       if (use_clm_microbe_dom_relaxation) then
          call advanceMicrobeMethaneDOMRelaxation(dom_c, dom_n, dom_p, &
               layer_thickness, dom_relaxation_rate, dt, updated_dom_c, &
               updated_dom_n, updated_dom_p, dom_carbon_residual, &
               dom_nitrogen_residual, dom_phosphorus_residual, dom_relaxation_valid)
          dom_c = updated_dom_c
          dom_n = updated_dom_n
          dom_p = updated_dom_p
       else if (use_microbe_aqueous_transport) then
          call advanceMicrobeAqueousTracerTransport(dom_c, layer_thickness, &
               liquid_fraction, dom_diffusion_conductivity, water_flux, &
               MicrobeMethaneParamsInst%aqueous_dom_mobile_fraction, &
               MicrobeMethaneParamsInst%aqueous_solute_min_liquid_fraction, dt, &
               updated_dom_c, solute_advective_flux, solute_diffusive_flux, &
               solute_tendency, dom_carbon_export, dom_carbon_residual, &
               dom_carbon_transport_valid)
          this%dom_advective_flux_col(c,1:nlevdecomp) = &
               solute_advective_flux(1:nlevdecomp)
          this%dom_diffusive_flux_col(c,1:nlevdecomp) = &
               solute_diffusive_flux(1:nlevdecomp)
          this%dom_bottom_carbon_export_col(c) = &
               solute_advective_flux(nlevdecomp) + solute_diffusive_flux(nlevdecomp)
          call advanceMicrobeAqueousTracerTransport(dom_n, layer_thickness, &
               liquid_fraction, dom_diffusion_conductivity, water_flux, &
               MicrobeMethaneParamsInst%aqueous_dom_mobile_fraction, &
               MicrobeMethaneParamsInst%aqueous_solute_min_liquid_fraction, dt, &
               updated_dom_n, solute_advective_flux, solute_diffusive_flux, &
               solute_tendency, dom_nitrogen_export, dom_nitrogen_residual, &
               dom_nitrogen_transport_valid)
          this%dom_bottom_nitrogen_export_col(c) = &
               solute_advective_flux(nlevdecomp) + solute_diffusive_flux(nlevdecomp)
          call advanceMicrobeAqueousTracerTransport(dom_p, layer_thickness, &
               liquid_fraction, dom_diffusion_conductivity, water_flux, &
               MicrobeMethaneParamsInst%aqueous_dom_mobile_fraction, &
               MicrobeMethaneParamsInst%aqueous_solute_min_liquid_fraction, dt, &
               updated_dom_p, solute_advective_flux, solute_diffusive_flux, &
               solute_tendency, dom_phosphorus_export, dom_phosphorus_residual, &
               dom_phosphorus_transport_valid)
          this%dom_bottom_phosphorus_export_col(c) = &
               solute_advective_flux(nlevdecomp) + solute_diffusive_flux(nlevdecomp)
          dom_c = updated_dom_c
          dom_n = updated_dom_n
          dom_p = updated_dom_p
          this%aqueous_carbon_export_col(c) = dom_carbon_export
          this%aqueous_nitrogen_export_col(c) = dom_nitrogen_export
          this%aqueous_phosphorus_export_col(c) = dom_phosphorus_export
          this%dom_carbon_transport_residual_col(c) = dom_carbon_residual
          this%dom_nitrogen_transport_residual_col(c) = dom_nitrogen_residual
          this%dom_phosphorus_transport_residual_col(c) = dom_phosphorus_residual
       end if
       if (use_clm_microbe_dom_relaxation .or. use_microbe_aqueous_transport) then
          do j = 1, nlevdecomp
             ! DOM is one authoritative ELM pool; reaction partitions receive
             ! the same post-transport value but do not own duplicate storage.
             unsaturated_work(j)%dom_c = dom_c(j)
             saturated_work(j)%dom_c = dom_c(j)
          end do
       end if
       column_valid = column_valid .and. dom_relaxation_valid .and. &
            dom_carbon_transport_valid .and. dom_nitrogen_transport_valid .and. &
            dom_phosphorus_transport_valid

       do gas = 1, microbe_gas_count
          if (use_elm_microbe_methane_transport) then
             ! The optional ELM mapping uses a gas-equivalent mobile
             ! concentration. Henry-law capacities map that concentration
             ! back to the conserved bulk-soil inventories.
             surface_equilibrium(gas) = atmosphericGasConcentration( &
                  partial_pressure(gas), atm2lnd_vars%forc_t_downscaled_col(c))
             unsaturated_surface_conductance(gas) = surfaceConductance(c, gas, .false., &
                  ground_conductance_col(c), unsaturated_diffusivity(1,gas))
             saturated_surface_conductance(gas) = surfaceConductance(c, gas, .true., &
                  ground_conductance_col(c), saturated_diffusivity(1,gas))
          else
             ! The default CLM-Microbe mapping transports the stored molar
             ! concentration directly and couples its top layer to the
             ! atmospheric Henry-law aqueous equilibrium with a Fickian
             ! half-layer conductance. This retains conservative implicit
             ! numerics instead of the source's timestep reset.
             surface_equilibrium(gas) = henryEquilibriumConcentration( &
                  partial_pressure(gas), col_es%t_soisno(c,1), gas)
             unsaturated_surface_conductance(gas) = legacyFickianSurfaceConductance( &
                  unsaturated_diffusivity(1,gas), layer_thickness(1), &
                  col_es%t_soisno(c,1))
             saturated_surface_conductance(gas) = legacyFickianSurfaceConductance( &
                  saturated_diffusivity(1,gas), layer_thickness(1), &
                  col_es%t_soisno(c,1))
          end if
       end do

       if (use_microbe_aqueous_transport) then
          do j = 1, nlevdecomp
             acetate_concentration(j) = unsaturated_work(j)%acetate_c
          end do
          call advanceMicrobeAqueousTracerTransport(acetate_concentration, &
               layer_thickness, liquid_fraction, acetate_diffusion_conductivity, &
               water_flux, MicrobeMethaneParamsInst%aqueous_acetate_mobile_fraction, &
               MicrobeMethaneParamsInst%aqueous_solute_min_liquid_fraction, dt, &
               updated_acetate, solute_advective_flux, solute_diffusive_flux, &
               solute_tendency, unsaturated_acetate_export, &
               unsaturated_acetate_residual, unsaturated_acetate_valid)
          unsaturated_candidate = unsaturated_work
          do j = 1, nlevdecomp
             unsaturated_candidate(j)%acetate_c = updated_acetate(j)
          end do

          do j = 1, nlevdecomp
             acetate_concentration(j) = saturated_work(j)%acetate_c
          end do
          call advanceMicrobeAqueousTracerTransport(acetate_concentration, &
               layer_thickness, liquid_fraction, acetate_diffusion_conductivity, &
               water_flux, MicrobeMethaneParamsInst%aqueous_acetate_mobile_fraction, &
               MicrobeMethaneParamsInst%aqueous_solute_min_liquid_fraction, dt, &
               updated_acetate, solute_advective_flux, solute_diffusive_flux, &
               solute_tendency, saturated_acetate_export, &
               saturated_acetate_residual, saturated_acetate_valid)
          saturated_candidate = saturated_work
          do j = 1, nlevdecomp
             saturated_candidate(j)%acetate_c = updated_acetate(j)
          end do
          this%aqueous_carbon_export_col(c) = this%aqueous_carbon_export_col(c) + &
               (1._r8 - fraction) * unsaturated_acetate_export + &
               fraction * saturated_acetate_export
       else
          unsaturated_acetate_export = 0._r8
          saturated_acetate_export = 0._r8
          call advanceMicrobeMethaneAcetateTransport(unsaturated_work, layer_thickness, &
               unsaturated_acetate_diffusivity, dt, unsaturated_candidate, &
               acetate_interface_flux, acetate_tendency, unsaturated_acetate_residual, &
               unsaturated_acetate_valid)
          call advanceMicrobeMethaneAcetateTransport(saturated_work, layer_thickness, &
               saturated_acetate_diffusivity, dt, saturated_candidate, &
               acetate_interface_flux, acetate_tendency, saturated_acetate_residual, &
               saturated_acetate_valid)
       end if
       column_valid = column_valid .and. unsaturated_acetate_valid .and. &
            saturated_acetate_valid

       call advanceMicrobeMethaneGasTransport(unsaturated_candidate, layer_thickness, &
            unsaturated_diffusivity, unsaturated_transport_capacity, surface_equilibrium, &
            unsaturated_surface_conductance, aerenchyma_equilibrium, &
            aerenchyma_exchange_rate, aerenchyma_minimum_emission, &
            aerenchyma_allows_influx, ch4_ebullition_threshold, &
            unsaturated_ebullition_activation, dt, unsaturated_work, interface_flux, &
            aerenchyma_flux, ch4_ebullition_loss, unsaturated_surface_diffusive_flux, &
            unsaturated_surface_aerenchyma_flux, unsaturated_surface_ebullition_flux, &
            unsaturated_surface_flux, unsaturated_gas_residual, unsaturated_gas_valid)
       call advanceMicrobeMethaneGasTransport(saturated_candidate, layer_thickness, &
            saturated_diffusivity, saturated_transport_capacity, surface_equilibrium, &
            saturated_surface_conductance, &
            aerenchyma_equilibrium, aerenchyma_exchange_rate, &
            aerenchyma_minimum_emission, aerenchyma_allows_influx, &
            ch4_ebullition_threshold, &
            saturated_ebullition_activation, dt, saturated_work, interface_flux, &
            aerenchyma_flux, ch4_ebullition_loss, saturated_surface_diffusive_flux, &
            saturated_surface_aerenchyma_flux, saturated_surface_ebullition_flux, &
            saturated_surface_flux, saturated_gas_residual, saturated_gas_valid)
       column_valid = column_valid .and. unsaturated_gas_valid .and. saturated_gas_valid
       bulk_surface_flux = (1._r8 - fraction) * unsaturated_surface_flux + &
            fraction * saturated_surface_flux

       if (.not. column_valid) then
          do j = 1, nlevdecomp
             if (.not. reaction(j)%valid) then
                write(iulog,*) 'invalid revised methane reaction layer: ', j
                write(iulog,*) 'reaction C/N/P residuals: ', reaction(j)%carbon_residual, &
                     reaction(j)%nitrogen_residual, reaction(j)%phosphorus_residual
                write(iulog,*) 'reaction DOM C/N/P and mineral N/P: ', reaction(j)%dom_c, &
                     reaction(j)%dom_n, reaction(j)%dom_p, reaction(j)%mineral_n, &
                     reaction(j)%mineral_p
                call reportReactionState('unsaturated reaction state', &
                     reaction(j)%unsaturated_state)
                call reportReactionState('saturated reaction state', &
                     reaction(j)%saturated_state)
             end if
          end do
          write(iulog,*) 'revised methane atmospheric partial pressures: ', partial_pressure
          write(iulog,*) 'revised methane surface equilibrium concentrations: ', surface_equilibrium
          write(iulog,*) 'revised methane unsaturated surface conductance: ', &
               unsaturated_surface_conductance
          write(iulog,*) 'revised methane saturated surface conductance: ', &
               saturated_surface_conductance
          write(iulog,*) 'revised methane top aerenchyma equilibrium: ', &
               aerenchyma_equilibrium(1,:)
          write(iulog,*) 'revised methane top aerenchyma exchange rate: ', &
               aerenchyma_exchange_rate(1,:)
          write(iulog,*) 'revised methane bulk surface flux: ', bulk_surface_flux
          write(iulog,*) 'revised methane DOM transport validity C/N/P: ', &
               dom_carbon_transport_valid, dom_nitrogen_transport_valid, &
               dom_phosphorus_transport_valid
          write(iulog,*) 'revised methane unsaturated top state before gas transport: ', &
               unsaturated_candidate(1)%conc_ch4, unsaturated_candidate(1)%conc_o2, &
               unsaturated_candidate(1)%conc_co2, unsaturated_candidate(1)%conc_h2
          write(iulog,*) 'revised methane unsaturated top state after gas transport: ', &
               unsaturated_work(1)%conc_ch4, unsaturated_work(1)%conc_o2, &
               unsaturated_work(1)%conc_co2, unsaturated_work(1)%conc_h2
          write(message,'(a,i0,6(a,l1),7(a,es12.4))') &
               ' ERROR: revised methane column transaction failed validation for column ', c, &
               '; reaction=', all(reaction(:)%valid), &
               '; dom_relaxation=', dom_relaxation_valid, &
               '; acetate_unsat=', unsaturated_acetate_valid, &
               '; acetate_sat=', saturated_acetate_valid, &
               '; gas_unsat=', unsaturated_gas_valid, &
               '; gas_sat=', saturated_gas_valid, &
               '; dom_C_residual=', dom_carbon_residual, &
               '; dom_N_residual=', dom_nitrogen_residual, &
               '; dom_P_residual=', dom_phosphorus_residual, &
               '; acetate_residual_unsat=', unsaturated_acetate_residual, &
               '; acetate_residual_sat=', saturated_acetate_residual, &
               '; gas_residual_unsat=', unsaturated_gas_residual, &
               '; gas_residual_sat=', saturated_gas_residual
          call endrun(msg=trim(message)//errMsg(__FILE__, __LINE__))
       end if

       call commitColumnState(c, unsaturated_work, saturated_work)
       do j = 1, nlevdecomp
          col_cs%decomp_cpools_vr(c,j,i_dom) = dom_c(j)
          col_ns%decomp_npools_vr(c,j,i_dom) = dom_n(j)
          col_ps%decomp_ppools_vr(c,j,i_dom) = dom_p(j)
          col_ns%smin_nh4_vr(c,j) = mineral_n(j)
          col_ns%sminn_vr(c,j) = mineral_n(j) + col_ns%smin_no3_vr(c,j)
          col_ps%solutionp_vr(c,j) = mineral_p(j)
       end do
       this%sat_fraction_previous_col(c) = fraction
       this%additional_carbon_col(c) = microbeMethaneColumnAdditionalCarbon( &
            fraction, unsaturated_work, saturated_work, layer_thickness)
       this%surface_carbon_flux_col(c) = microbeMethaneSurfaceCarbonFlux(bulk_surface_flux)
       this%surface_ch4_flux_col(c) = microbeMethaneCH4SurfaceFluxKgC(bulk_surface_flux)
       this%surface_co2_flux_col(c) = microbeMethaneCO2Correction(bulk_surface_flux)
       this%ch4_surface_diffusion_unsat_col(c) = catomw * (1._r8 - fraction) * &
            unsaturated_surface_diffusive_flux(microbe_gas_ch4)
       this%ch4_surface_diffusion_sat_col(c) = catomw * fraction * &
            saturated_surface_diffusive_flux(microbe_gas_ch4)
       this%ch4_surface_diffusion_col(c) = this%ch4_surface_diffusion_unsat_col(c) + &
            this%ch4_surface_diffusion_sat_col(c)
       this%ch4_surface_aerenchyma_unsat_col(c) = catomw * (1._r8 - fraction) * &
            unsaturated_surface_aerenchyma_flux(microbe_gas_ch4)
       this%ch4_surface_aerenchyma_sat_col(c) = catomw * fraction * &
            saturated_surface_aerenchyma_flux(microbe_gas_ch4)
       this%ch4_surface_aerenchyma_col(c) = this%ch4_surface_aerenchyma_unsat_col(c) + &
            this%ch4_surface_aerenchyma_sat_col(c)
       this%ch4_surface_ebullition_unsat_col(c) = catomw * (1._r8 - fraction) * &
            unsaturated_surface_ebullition_flux
       this%ch4_surface_ebullition_sat_col(c) = catomw * fraction * &
            saturated_surface_ebullition_flux
       this%ch4_surface_ebullition_col(c) = this%ch4_surface_ebullition_unsat_col(c) + &
            this%ch4_surface_ebullition_sat_col(c)
       do j = 1, nlevdecomp
          this%ch4_production_unsat_col(c) = this%ch4_production_unsat_col(c) + &
               catomw * layer_thickness(j) * (1._r8 - fraction) * &
               (MicrobeMethaneParamsInst%acetoclastic_methanogenesis_ch4_yield * &
               (1._r8 - MicrobeMethaneParamsInst%acetate_methanogen_yield) * &
               reaction(j)%unsaturated_rates%acetoclastic_methanogenesis_c + &
               reaction(j)%unsaturated_rates%hydrogenotrophic_methanogenesis_c)
          this%ch4_production_sat_col(c) = this%ch4_production_sat_col(c) + &
               catomw * layer_thickness(j) * fraction * &
               (MicrobeMethaneParamsInst%acetoclastic_methanogenesis_ch4_yield * &
               (1._r8 - MicrobeMethaneParamsInst%acetate_methanogen_yield) * &
               reaction(j)%saturated_rates%acetoclastic_methanogenesis_c + &
               reaction(j)%saturated_rates%hydrogenotrophic_methanogenesis_c)
          this%ch4_aerobic_oxidation_unsat_col(c) = &
               this%ch4_aerobic_oxidation_unsat_col(c) + &
               catomw * layer_thickness(j) * (1._r8 - fraction) * &
               reaction(j)%unsaturated_rates%aerobic_methane_oxidation_c
          this%ch4_aerobic_oxidation_sat_col(c) = &
               this%ch4_aerobic_oxidation_sat_col(c) + &
               catomw * layer_thickness(j) * fraction * &
               reaction(j)%saturated_rates%aerobic_methane_oxidation_c
          this%ch4_aerobic_oxidation_pre_o2_unsat_col(c) = &
               this%ch4_aerobic_oxidation_pre_o2_unsat_col(c) + &
               catomw * layer_thickness(j) * (1._r8 - fraction) * &
               reaction(j)%unsaturated_rates%aerobic_methane_oxidation_pre_o2_c
          this%ch4_aerobic_oxidation_pre_o2_sat_col(c) = &
               this%ch4_aerobic_oxidation_pre_o2_sat_col(c) + &
               catomw * layer_thickness(j) * fraction * &
               reaction(j)%saturated_rates%aerobic_methane_oxidation_pre_o2_c
          this%ch4_anaerobic_oxidation_unsat_col(c) = &
               this%ch4_anaerobic_oxidation_unsat_col(c) + &
               catomw * layer_thickness(j) * (1._r8 - fraction) * &
               reaction(j)%unsaturated_rates%anaerobic_methane_oxidation_c
          this%ch4_anaerobic_oxidation_sat_col(c) = &
               this%ch4_anaerobic_oxidation_sat_col(c) + &
               catomw * layer_thickness(j) * fraction * &
               reaction(j)%saturated_rates%anaerobic_methane_oxidation_c
       end do
       this%ch4_production_col(c) = this%ch4_production_unsat_col(c) + &
            this%ch4_production_sat_col(c)
       this%ch4_aerobic_oxidation_col(c) = this%ch4_aerobic_oxidation_unsat_col(c) + &
            this%ch4_aerobic_oxidation_sat_col(c)
       this%ch4_aerobic_oxidation_pre_o2_col(c) = &
            this%ch4_aerobic_oxidation_pre_o2_unsat_col(c) + &
            this%ch4_aerobic_oxidation_pre_o2_sat_col(c)
       this%ch4_anaerobic_oxidation_col(c) = this%ch4_anaerobic_oxidation_unsat_col(c) + &
            this%ch4_anaerobic_oxidation_sat_col(c)
       this%ch4_oxidation_unsat_col(c) = this%ch4_aerobic_oxidation_unsat_col(c) + &
            this%ch4_anaerobic_oxidation_unsat_col(c)
       this%ch4_oxidation_sat_col(c) = this%ch4_aerobic_oxidation_sat_col(c) + &
            this%ch4_anaerobic_oxidation_sat_col(c)
       this%ch4_oxidation_col(c) = this%ch4_oxidation_unsat_col(c) + &
            this%ch4_oxidation_sat_col(c)
       call syncLegacyOxygenBridge(c, fraction, unsaturated_work, saturated_work, reaction)
    end do

    ! CLM-SPRUCE's Fickian backend also exchanged gas laterally between its
    ! two hard-coded hummock/hollow columns. Apply the generalized topounit
    ! graph transaction only in the parity transport mode; the optional ELM
    ! multiphase mapping intentionally remains the pre-existing alternative.
    if (use_humhol .and. .not. use_elm_microbe_methane_transport) then
       call applyTopounitLateralGasDiffusion()
    end if

  contains

    subroutine applyTopounitLateralGasDiffusion()
      real(r8) :: lateral_concentration(num_soilc,nlevdecomp)
      real(r8) :: lateral_updated(num_soilc,nlevdecomp)
      real(r8) :: lateral_storage_weight(num_soilc,nlevdecomp)
      real(r8) :: lateral_layer_top(num_soilc,nlevdecomp)
      real(r8) :: lateral_layer_bottom(num_soilc,nlevdecomp)
      real(r8) :: lateral_diffusivity(num_soilc,nlevdecomp)
      real(r8) :: lateral_edge_distance(num_soilc)
      real(r8) :: additional_carbon_before(num_soilc)
      real(r8) :: lateral_residual, partition_fraction, updated_additional_carbon
      integer :: lateral_edge_source(num_soilc), lateral_edge_target(num_soilc)
      integer :: soil_fc_by_topounit(bounds%begt:bounds%endt)
      integer :: edge_count, edge, source_fc, target_fc, target_t, t, partition
      logical :: lateral_valid, duplicate_edge

      if (num_soilc <= 0) return
      soil_fc_by_topounit = 0
      do source_fc = 1, num_soilc
         c = filter_soilc(source_fc)
         t = col_pp%topounit(c)
         if (t < bounds%begt .or. t > bounds%endt) cycle
         if (soil_fc_by_topounit(t) /= 0 .and. &
              soil_fc_by_topounit(t) /= source_fc) then
            call endrun(msg=' ERROR: revised methane lateral transport requires one '//&
                 'soil column per connected topounit'//errMsg(__FILE__, __LINE__))
         end if
         soil_fc_by_topounit(t) = source_fc
      end do

      edge_count = 0
      do t = bounds%begt, bounds%endt
         if (.not. top_pp%active(t) .or. soil_fc_by_topounit(t) == 0) cycle
         target_t = top_pp%regional_target_ti(t)
         if (target_t < bounds%begt .or. target_t > bounds%endt .or. &
              target_t == t) cycle
         if (.not. top_pp%active(target_t) .or. &
              top_pp%gridcell(target_t) /= top_pp%gridcell(t)) cycle
         ! Standalone SPRUCE includes a large non-bog boardwalk/fen topounit.
         ! Keep its legacy-compatible gas edge by default, but allow a
         ! diagnostic to retain bog hollow--hummock exchange while excluding
         ! gas transfer into or out of the non-bog unit.
         if (.not. use_microbe_nonbog_lateral_gas_transport) then
            if (.not. top_pp%is_bog(t) .or. .not. top_pp%is_bog(target_t)) cycle
         end if
         if (soil_fc_by_topounit(target_t) == 0 .or. &
              top_pp%lateral_dist(t) <= 0._r8) cycle
         source_fc = soil_fc_by_topounit(t)
         target_fc = soil_fc_by_topounit(target_t)
         duplicate_edge = .false.
         do edge = 1, edge_count
            if ((lateral_edge_source(edge) == source_fc .and. &
                 lateral_edge_target(edge) == target_fc) .or. &
                 (lateral_edge_source(edge) == target_fc .and. &
                 lateral_edge_target(edge) == source_fc)) then
               duplicate_edge = .true.
            end if
         end do
         if (duplicate_edge) cycle
         edge_count = edge_count + 1
         lateral_edge_source(edge_count) = source_fc
         lateral_edge_target(edge_count) = target_fc
         lateral_edge_distance(edge_count) = top_pp%lateral_dist(t)
      end do
      if (edge_count == 0) return

      do source_fc = 1, num_soilc
         c = filter_soilc(source_fc)
         t = col_pp%topounit(c)
         additional_carbon_before(source_fc) = this%additional_carbon_col(c)
         do j = 1, nlevdecomp
            lateral_layer_top(source_fc,j) = top_pp%elevation(t) - col_pp%zi(c,j-1)
            lateral_layer_bottom(source_fc,j) = top_pp%elevation(t) - col_pp%zi(c,j)
         end do
      end do

      do partition = 1, 2
         do gas = 1, microbe_gas_count
            do source_fc = 1, num_soilc
               c = filter_soilc(source_fc)
               fraction = clampUnitInterval(this%sat_fraction_previous_col(c))
               if (partition == 1) then
                  partition_fraction = 1._r8 - fraction
               else
                  partition_fraction = fraction
               end if
               do j = 1, nlevdecomp
                  lateral_storage_weight(source_fc,j) = &
                       max(0._r8, top_pp%wtgcell(col_pp%topounit(c))) * &
                       max(0._r8, col_pp%wttopounit(c)) * partition_fraction * &
                       col_pp%dz(c,j)
                  lateral_diffusivity(source_fc,j) = &
                       legacyFickianLateralGasDiffusivity(gas, col_es%t_soisno(c,j))
                  select case (gas)
                  case (microbe_gas_ch4)
                     if (partition == 1) then
                        lateral_concentration(source_fc,j) = this%conc_ch4_unsat_col(c,j)
                     else
                        lateral_concentration(source_fc,j) = this%conc_ch4_sat_col(c,j)
                     end if
                  case (microbe_gas_o2)
                     if (partition == 1) then
                        lateral_concentration(source_fc,j) = this%conc_o2_unsat_col(c,j)
                     else
                        lateral_concentration(source_fc,j) = this%conc_o2_sat_col(c,j)
                     end if
                  case (microbe_gas_co2)
                     if (partition == 1) then
                        lateral_concentration(source_fc,j) = this%conc_co2_unsat_col(c,j)
                     else
                        lateral_concentration(source_fc,j) = this%conc_co2_sat_col(c,j)
                     end if
                  case (microbe_gas_h2)
                     if (partition == 1) then
                        lateral_concentration(source_fc,j) = this%conc_h2_unsat_col(c,j)
                     else
                        lateral_concentration(source_fc,j) = this%conc_h2_sat_col(c,j)
                     end if
                  end select
               end do
            end do

            call computeMicrobeTopounitLateralDiffusion(lateral_concentration, &
                 lateral_storage_weight, lateral_layer_top, lateral_layer_bottom, &
                 lateral_diffusivity, lateral_edge_source(1:edge_count), &
                 lateral_edge_target(1:edge_count), &
                 lateral_edge_distance(1:edge_count), dt, lateral_updated, &
                 lateral_residual, lateral_valid)
            if (.not. lateral_valid) then
               write(message,'(a,i0,a,i0,a,es24.16)') &
                    ' ERROR: invalid revised methane lateral transport partition ', &
                    partition, ' gas ', gas, ' residual ', lateral_residual
               call endrun(msg=trim(message)//errMsg(__FILE__, __LINE__))
            end if

            do source_fc = 1, num_soilc
               c = filter_soilc(source_fc)
               do j = 1, nlevdecomp
                  select case (gas)
                  case (microbe_gas_ch4)
                     if (partition == 1) then
                        this%conc_ch4_unsat_col(c,j) = lateral_updated(source_fc,j)
                     else
                        this%conc_ch4_sat_col(c,j) = lateral_updated(source_fc,j)
                     end if
                  case (microbe_gas_o2)
                     if (partition == 1) then
                        this%conc_o2_unsat_col(c,j) = lateral_updated(source_fc,j)
                     else
                        this%conc_o2_sat_col(c,j) = lateral_updated(source_fc,j)
                     end if
                  case (microbe_gas_co2)
                     if (partition == 1) then
                        this%conc_co2_unsat_col(c,j) = lateral_updated(source_fc,j)
                     else
                        this%conc_co2_sat_col(c,j) = lateral_updated(source_fc,j)
                     end if
                  case (microbe_gas_h2)
                     if (partition == 1) then
                        this%conc_h2_unsat_col(c,j) = lateral_updated(source_fc,j)
                     else
                        this%conc_h2_sat_col(c,j) = lateral_updated(source_fc,j)
                     end if
                  end select
               end do
            end do
         end do
      end do

      ! Lateral exchange is internal to the gridcell. Record its signed
      ! column-level carbon transfer separately so ELM's per-column balance
      ! check can distinguish it from an atmospheric surface flux.
      do source_fc = 1, num_soilc
         c = filter_soilc(source_fc)
         fraction = clampUnitInterval(this%sat_fraction_previous_col(c))
         call gatherColumnState(c, unsaturated_state, saturated_state)
         do j = 1, nlevdecomp
            layer_thickness(j) = col_pp%dz(c,j)
            ch4_vars%conc_o2_unsat_col(c,j) = max(0._r8, unsaturated_state(j)%conc_o2)
            ch4_vars%conc_o2_sat_col(c,j) = max(0._r8, saturated_state(j)%conc_o2)
         end do
         updated_additional_carbon = microbeMethaneColumnAdditionalCarbon( &
              fraction, unsaturated_state, saturated_state, layer_thickness)
         this%lateral_carbon_flux_col(c) = &
              (updated_additional_carbon - additional_carbon_before(source_fc)) / dt
         this%additional_carbon_col(c) = updated_additional_carbon
      end do
    end subroutine applyTopounitLateralGasDiffusion

    subroutine reportReactionState(label, state)
      character(len=*), intent(in) :: label
      type(microbe_methane_reaction_state_type), intent(in) :: state

      write(iulog,*) trim(label)//' DOM/acetate/guild C: ', state%dom_c, &
           state%acetate_c, state%acetate_methanogen_c, state%h2_methanogen_c, &
           state%aerobic_methanotroph_c, state%anaerobic_methanotroph_c
      write(iulog,*) trim(label)//' CH4/O2/CO2/H2: ', state%conc_ch4, &
           state%conc_o2, state%conc_co2, state%conc_h2
    end subroutine reportReactionState

    real(r8) function elmAerobicOxygenDemand(column, layer, root_respiration) result(potential_demand)
      ! CH4Mod-compatible potential O2 demand [mol O2 m-3 s-1]. The CTC
      ! potential HR still contains the preceding o_scalar, so divide it off
      ! before forming the new shared stress.
      integer, intent(in) :: column, layer
      real(r8), intent(in) :: root_respiration

      potential_demand = max(0._r8, col_cf%phr_vr(column,layer)) / catomw
      if (col_cf%o_scalar(column,layer) > tiny(1._r8)) then
         potential_demand = potential_demand / col_cf%o_scalar(column,layer)
      end if
      potential_demand = potential_demand + &
           max(0._r8, root_respiration) / &
           (catomw * max(col_pp%dz(column,layer), tiny(1._r8)))
      potential_demand = potential_demand + &
           max(0._r8, col_nf%pot_f_nit_vr(column,layer)) * 2._r8 / 14._r8
    end function elmAerobicOxygenDemand

    subroutine syncLegacyOxygenBridge(column, saturated_fraction, unsaturated, saturated, reactions)
      ! Decomposition and nitrification still consume the CH4Mod oxygen
      ! interface.  In revised mode ch4_vars is only a compatibility carrier;
      ! the legacy CH4 solver is not executed.  Publish the revised O2 state
      ! and a lagged estimate of aerobic demand for the next ELM timestep.
      integer, intent(in) :: column
      real(r8), intent(in) :: saturated_fraction
      type(microbe_methane_reaction_state_type), intent(in) :: unsaturated(nlevdecomp)
      type(microbe_methane_reaction_state_type), intent(in) :: saturated(nlevdecomp)
      type(microbe_methane_reaction_transaction_type), intent(in) :: reactions(nlevdecomp)
      real(r8) :: potential_demand, stress_unsaturated, stress_saturated
      integer :: layer

      ch4_vars%finundated_col(column) = saturated_fraction
      do layer = 1, nlevdecomp
         potential_demand = elmAerobicOxygenDemand(column, layer, &
              root_respiration_col_vr(column,layer))
         stress_unsaturated = reactions(layer)%unsaturated_rates%oxygen_stress
         stress_saturated = reactions(layer)%saturated_rates%oxygen_stress

         ! The conservative reaction/transport transactions permit roundoff
         ! down to -state_tolerance. NitrifDenitrifMod subsequently raises
         ! this compatibility concentration to a fractional power, for which
         ! even a tiny negative value is undefined. Publish the physical
         ! nonnegative concentration without altering the transaction state
         ! or its mass-balance accounting.
         ch4_vars%conc_o2_unsat_col(column,layer) = &
              max(0._r8, unsaturated(layer)%conc_o2)
         ch4_vars%conc_o2_sat_col(column,layer) = &
              max(0._r8, saturated(layer)%conc_o2)
         ch4_vars%o2stress_unsat_col(column,layer) = stress_unsaturated
         ch4_vars%o2stress_sat_col(column,layer) = stress_saturated
         ch4_vars%o2_decomp_depth_unsat_col(column,layer) = &
              potential_demand * stress_unsaturated
         ch4_vars%o2_decomp_depth_sat_col(column,layer) = &
              potential_demand * stress_saturated
      end do
    end subroutine syncLegacyOxygenBridge

    subroutine gatherColumnState(column, unsaturated, saturated)
      integer, intent(in) :: column
      type(microbe_methane_reaction_state_type), intent(out) :: unsaturated(nlevdecomp)
      type(microbe_methane_reaction_state_type), intent(out) :: saturated(nlevdecomp)
      integer :: layer
      do layer = 1, nlevdecomp
         unsaturated(layer)%acetate_c = this%acetate_c_unsat_col(column,layer)
         saturated(layer)%acetate_c = this%acetate_c_sat_col(column,layer)
         unsaturated(layer)%acetate_methanogen_c = &
              this%acetate_methanogen_c_unsat_col(column,layer)
         saturated(layer)%acetate_methanogen_c = &
              this%acetate_methanogen_c_sat_col(column,layer)
         unsaturated(layer)%h2_methanogen_c = this%h2_methanogen_c_unsat_col(column,layer)
         saturated(layer)%h2_methanogen_c = this%h2_methanogen_c_sat_col(column,layer)
         unsaturated(layer)%aerobic_methanotroph_c = &
              this%aerobic_methanotroph_c_unsat_col(column,layer)
         saturated(layer)%aerobic_methanotroph_c = &
              this%aerobic_methanotroph_c_sat_col(column,layer)
         unsaturated(layer)%anaerobic_methanotroph_c = &
              this%anaerobic_methanotroph_c_unsat_col(column,layer)
         saturated(layer)%anaerobic_methanotroph_c = &
              this%anaerobic_methanotroph_c_sat_col(column,layer)
         unsaturated(layer)%conc_ch4 = this%conc_ch4_unsat_col(column,layer)
         saturated(layer)%conc_ch4 = this%conc_ch4_sat_col(column,layer)
         unsaturated(layer)%conc_o2 = this%conc_o2_unsat_col(column,layer)
         saturated(layer)%conc_o2 = this%conc_o2_sat_col(column,layer)
         unsaturated(layer)%conc_co2 = this%conc_co2_unsat_col(column,layer)
         saturated(layer)%conc_co2 = this%conc_co2_sat_col(column,layer)
         unsaturated(layer)%conc_h2 = this%conc_h2_unsat_col(column,layer)
         saturated(layer)%conc_h2 = this%conc_h2_sat_col(column,layer)
      end do
    end subroutine gatherColumnState

    subroutine repartitionColumnState(old_sat, new_sat, unsaturated, saturated)
      real(r8), intent(in) :: old_sat, new_sat
      type(microbe_methane_reaction_state_type), intent(inout) :: unsaturated(nlevdecomp)
      type(microbe_methane_reaction_state_type), intent(inout) :: saturated(nlevdecomp)
      integer :: layer
      do layer = 1, nlevdecomp
         call repartitionLocalPair(old_sat, new_sat, unsaturated(layer)%acetate_c, &
              saturated(layer)%acetate_c)
         call repartitionLocalPair(old_sat, new_sat, &
              unsaturated(layer)%acetate_methanogen_c, &
              saturated(layer)%acetate_methanogen_c)
         call repartitionLocalPair(old_sat, new_sat, unsaturated(layer)%h2_methanogen_c, &
              saturated(layer)%h2_methanogen_c)
         call repartitionLocalPair(old_sat, new_sat, &
              unsaturated(layer)%aerobic_methanotroph_c, &
              saturated(layer)%aerobic_methanotroph_c)
         call repartitionLocalPair(old_sat, new_sat, &
              unsaturated(layer)%anaerobic_methanotroph_c, &
              saturated(layer)%anaerobic_methanotroph_c)
         call repartitionLocalPair(old_sat, new_sat, unsaturated(layer)%conc_ch4, &
              saturated(layer)%conc_ch4)
         call repartitionLocalPair(old_sat, new_sat, unsaturated(layer)%conc_o2, &
              saturated(layer)%conc_o2)
         call repartitionLocalPair(old_sat, new_sat, unsaturated(layer)%conc_co2, &
              saturated(layer)%conc_co2)
         call repartitionLocalPair(old_sat, new_sat, unsaturated(layer)%conc_h2, &
              saturated(layer)%conc_h2)
      end do
    end subroutine repartitionColumnState

    subroutine repartitionLocalPair(old_sat, new_sat, unsaturated_value, saturated_value)
      real(r8), intent(in) :: old_sat, new_sat
      real(r8), intent(inout) :: unsaturated_value, saturated_value
      real(r8) :: new_unsaturated, new_saturated
      call repartitionMicrobeMethaneScalar(old_sat, new_sat, unsaturated_value, &
           saturated_value, new_unsaturated, new_saturated)
      unsaturated_value = new_unsaturated
      saturated_value = new_saturated
    end subroutine repartitionLocalPair

    subroutine commitColumnState(column, unsaturated, saturated)
      integer, intent(in) :: column
      type(microbe_methane_reaction_state_type), intent(in) :: unsaturated(nlevdecomp)
      type(microbe_methane_reaction_state_type), intent(in) :: saturated(nlevdecomp)
      integer :: layer
      do layer = 1, nlevdecomp
         this%acetate_c_unsat_col(column,layer) = unsaturated(layer)%acetate_c
         this%acetate_c_sat_col(column,layer) = saturated(layer)%acetate_c
         this%acetate_methanogen_c_unsat_col(column,layer) = &
              unsaturated(layer)%acetate_methanogen_c
         this%acetate_methanogen_c_sat_col(column,layer) = &
              saturated(layer)%acetate_methanogen_c
         this%h2_methanogen_c_unsat_col(column,layer) = unsaturated(layer)%h2_methanogen_c
         this%h2_methanogen_c_sat_col(column,layer) = saturated(layer)%h2_methanogen_c
         this%aerobic_methanotroph_c_unsat_col(column,layer) = &
              unsaturated(layer)%aerobic_methanotroph_c
         this%aerobic_methanotroph_c_sat_col(column,layer) = &
              saturated(layer)%aerobic_methanotroph_c
         this%anaerobic_methanotroph_c_unsat_col(column,layer) = &
              unsaturated(layer)%anaerobic_methanotroph_c
         this%anaerobic_methanotroph_c_sat_col(column,layer) = &
              saturated(layer)%anaerobic_methanotroph_c
         this%conc_ch4_unsat_col(column,layer) = unsaturated(layer)%conc_ch4
         this%conc_ch4_sat_col(column,layer) = saturated(layer)%conc_ch4
         this%conc_o2_unsat_col(column,layer) = unsaturated(layer)%conc_o2
         this%conc_o2_sat_col(column,layer) = saturated(layer)%conc_o2
         this%conc_co2_unsat_col(column,layer) = unsaturated(layer)%conc_co2
         this%conc_co2_sat_col(column,layer) = saturated(layer)%conc_co2
         this%conc_h2_unsat_col(column,layer) = unsaturated(layer)%conc_h2
         this%conc_h2_sat_col(column,layer) = saturated(layer)%conc_h2
      end do
    end subroutine commitColumnState

    real(r8) function layerLiquidSaturation(column, layer) result(value)
      integer, intent(in) :: column, layer
      real(r8) :: porosity, liquid_volume
      porosity = max(soilstate_vars%watsat_col(column,layer), tiny(1._r8))
      liquid_volume = max(0._r8, col_ws%h2osoi_liq(column,layer)) / &
           (denh2o * max(col_pp%dz(column,layer), tiny(1._r8)))
      value = clampUnitInterval(liquid_volume / porosity)
    end function layerLiquidSaturation

    real(r8) function layerThawedFraction(column, layer) result(value)
      integer, intent(in) :: column, layer
      real(r8) :: liquid_volume, ice_volume
      liquid_volume = max(0._r8, col_ws%h2osoi_liq(column,layer)) / denh2o
      ice_volume = max(0._r8, col_ws%h2osoi_ice(column,layer)) / denice
      if (liquid_volume + ice_volume > tiny(1._r8)) then
         value = clampUnitInterval(liquid_volume / (liquid_volume + ice_volume))
      else
         value = 0._r8
      end if
    end function layerThawedFraction

    pure real(r8) function soilMoistureResponse(soil_water_potential, saturated_suction, &
         minimum_potential) result(value)
      real(r8), intent(in) :: soil_water_potential, saturated_suction, minimum_potential
      real(r8) :: wet_potential, bounded_potential, denominator
      wet_potential = -abs(saturated_suction) * 9.8e-6_r8
      if (minimum_potential >= 0._r8 .or. wet_potential >= 0._r8 .or. &
           minimum_potential >= wet_potential) then
         value = 0._r8
         return
      end if
      bounded_potential = min(wet_potential, max(minimum_potential, soil_water_potential))
      denominator = log(minimum_potential / wet_potential)
      value = clampUnitInterval(log(minimum_potential / bounded_potential) / denominator)
    end function soilMoistureResponse

    pure real(r8) function atmosphericPartialPressure(forcing, pressure, fallback_mixing_ratio) &
         result(value)
      real(r8), intent(in) :: forcing, pressure, fallback_mixing_ratio
      if (ieee_is_finite(forcing) .and. forcing > 0._r8) then
         value = forcing
      else if (ieee_is_finite(pressure) .and. pressure > 0._r8) then
         value = pressure * fallback_mixing_ratio
      else
         value = -1._r8
      end if
    end function atmosphericPartialPressure

    pure real(r8) function validRootFraction(value_in) result(value)
      real(r8), intent(in) :: value_in
      if (ieee_is_finite(value_in) .and. value_in >= 0._r8 .and. value_in <= 1._r8) then
         value = value_in
      else
         ! Missing patch roots denote a genuinely unvegetated column for this
         ! pathway, not an arbitrarily large aerenchyma conductance.
         value = 0._r8
      end if
    end function validRootFraction

    subroutine gasTransportProperties(column, layer, gas_index, saturated_partition, &
         capacity, transport_diffusivity)
      ! Recast ELM's CH4Mod gas/aqueous diffusion formulation onto the
      ! revised backend's conserved bulk-soil gas inventories. capacity maps
      ! gas-equivalent mobile concentration to bulk inventory. Aqueous
      ! diffusivity is multiplied by Henry solubility so every layer's flux
      ! is driven by the same gas-equivalent concentration.
      integer, intent(in) :: column, layer, gas_index
      logical, intent(in) :: saturated_partition
      real(r8), intent(out) :: capacity, transport_diffusivity
      real(r8) :: porosity, water_filled_fraction, liquid_fraction
      real(r8) :: air_fraction, relative_air_fraction, thawed_fraction
      real(r8) :: henry_solubility, organic_fraction, soil_structure_factor
      real(r8) :: aqueous_diffusivity
      logical :: aqueous_layer

      if (.not. use_elm_microbe_methane_transport) then
         ! Reproduce the active CLM-Microbe Fickian coefficient and its
         ! direct concentration-gradient basis in both area partitions.
         ! The source labels Fick_D_w as cm2 s-1 but multiplies it by 1e-3
         ! in the vertical equation; preserve that executed arithmetic here.
         capacity = 1._r8
         transport_diffusivity = legacyFickianGasDiffusivity(gas_index, &
              col_es%t_soisno(column,layer))
         return
      end if

      porosity = max(soilstate_vars%watsat_col(column,layer), tiny(1._r8))
      water_filled_fraction = min(porosity, max(0._r8, &
           col_ws%h2osoi_vol(column,layer)))
      thawed_fraction = layerThawedFraction(column, layer)
      liquid_fraction = water_filled_fraction * thawed_fraction
      air_fraction = max(0._r8, porosity - water_filled_fraction)
      henry_solubility = henryDimensionlessSolubility( &
           col_es%t_soisno(column,layer), gas_index)
      aqueous_layer = saturated_partition .or. water_filled_fraction >= &
           CH4ParamsInst%f_sat * porosity

      if (.not. aqueous_layer .and. air_fraction > 0._r8) then
         relative_air_fraction = clampUnitInterval(air_fraction / porosity)
         if (ParamsShareInst%organic_max > 0._r8) then
            organic_fraction = clampUnitInterval( &
                 soilstate_vars%cellorg_col(column,layer) / ParamsShareInst%organic_max)
         else
            organic_fraction = 1._r8
         end if
         soil_structure_factor = organic_fraction * &
              relative_air_fraction**(10._r8/3._r8) / porosity**2 + &
              (1._r8 - organic_fraction) * air_fraction**2 * &
              relative_air_fraction**(3._r8 / &
              max(soilstate_vars%bsw_col(column,layer), tiny(1._r8)))
         capacity = air_fraction + liquid_fraction * henry_solubility
         transport_diffusivity = referenceGasDiffusivity(gas_index, &
              col_es%t_soisno(column,layer)) * soil_structure_factor * &
              CH4ParamsInst%scale_factor_gasdiff
      else
         ! The saturated subarea and layers below the unsaturated-subarea
         ! water table carry gas in liquid-filled pore space.
         liquid_fraction = porosity * thawed_fraction
         capacity = liquid_fraction * henry_solubility
         aqueous_diffusivity = microbeMethaneEffectiveAqueousDiffusivity( &
              referenceAqueousDiffusivity(gas_index), col_es%t_soisno(column,layer), &
              porosity**CH4ParamsInst%satpow * thawed_fraction, &
              MicrobeMethaneParamsInst) * CH4ParamsInst%scale_factor_liqdiff
         transport_diffusivity = aqueous_diffusivity * henry_solubility
      end if

      capacity = max(capacity, tiny(1._r8))
      transport_diffusivity = max(0._r8, transport_diffusivity)
    end subroutine gasTransportProperties

    pure real(r8) function legacyFickianGasDiffusivity(gas_index, temperature) result(value)
      integer, intent(in) :: gas_index
      real(r8), intent(in) :: temperature
      real(r8) :: source_fick_d_w

      select case (gas_index)
      case (microbe_gas_ch4)
         source_fick_d_w = 1.49e-5_r8
      case (microbe_gas_o2)
         source_fick_d_w = 2.10e-5_r8
      case (microbe_gas_co2)
         source_fick_d_w = 1.92e-5_r8
      case (microbe_gas_h2)
         source_fick_d_w = 4.50e-5_r8
      case default
         source_fick_d_w = 0._r8
      end select
      value = source_fick_d_w * 1.e-3_r8 * &
           MicrobeMethaneParamsInst%aqueous_gas_diffusion_multiplier * &
           (max(0._r8, temperature) / &
           MicrobeMethaneParamsInst%aqueous_diffusion_t_ref) ** &
           MicrobeMethaneParamsInst%aqueous_diffusion_temperature_exponent
    end function legacyFickianGasDiffusivity

    pure real(r8) function legacyFickianLateralGasDiffusivity(gas_index, temperature) &
         result(value)
      integer, intent(in) :: gas_index
      real(r8), intent(in) :: temperature
      ! The source's lateral HUM_HOL equation uses 1e-4 where its vertical
      ! equation uses 1e-3. Preserve that factor-of-ten distinction.
      value = 0.1_r8 * legacyFickianGasDiffusivity(gas_index, temperature)
    end function legacyFickianLateralGasDiffusivity

    pure real(r8) function legacyFickianAcetateDiffusivity(temperature) result(value)
      real(r8), intent(in) :: temperature
      value = MicrobeMethaneParamsInst%dom_diffusivity * &
           (max(0._r8, temperature) / &
           MicrobeMethaneParamsInst%aqueous_diffusion_t_ref) ** &
           MicrobeMethaneParamsInst%aqueous_diffusion_temperature_exponent
    end function legacyFickianAcetateDiffusivity

    pure real(r8) function legacyFickianSurfaceConductance(top_diffusivity, &
         top_layer_thickness, top_temperature) result(value)
      real(r8), intent(in) :: top_diffusivity, top_layer_thickness, top_temperature
      value = 0._r8
      if (top_temperature < MicrobeMethaneParamsInst%transport_thaw_threshold) return
      if (top_diffusivity <= 0._r8 .or. top_layer_thickness <= 0._r8) return
      value = top_diffusivity / max(0.5_r8 * top_layer_thickness, tiny(1._r8))
    end function legacyFickianSurfaceConductance

    pure real(r8) function referenceAqueousDiffusivity(gas_index) result(value)
      integer, intent(in) :: gas_index
      real(r8) :: temperature_c
      temperature_c = MicrobeMethaneParamsInst%aqueous_diffusion_t_ref - 273.15_r8
      if (gas_index <= microbe_gas_co2) then
         value = max(0._r8, d_con_w(gas_index,1) + &
              d_con_w(gas_index,2) * temperature_c + &
              d_con_w(gas_index,3) * temperature_c**2) * 1.e-9_r8
      else
         ! ELM's shared table currently stops at CO2. This non-tunable H2
         ! physical constant is the CLM-SPRUCE Fick_D_w(4) value in SI.
         value = 4.5e-9_r8
      end if
    end function referenceAqueousDiffusivity

    pure real(r8) function henryEquilibriumConcentration(pressure_pa, temperature, &
         gas_index) result(value)
      real(r8), intent(in) :: pressure_pa, temperature
      integer, intent(in) :: gas_index
      real(r8) :: coefficient, reference_henry
      if (temperature <= 0._r8) then
         value = 0._r8
         return
      end if
      if (gas_index <= microbe_gas_co2) then
         coefficient = c_h_inv(gas_index)
         reference_henry = kh_theta(gas_index)
      else
         ! H2 counterparts of ELM's shared CH4/O2/CO2 constants.
         coefficient = 500._r8
         reference_henry = 1282.1_r8
      end if
      value = max(0._r8, pressure_pa) / 101325._r8 * 1000._r8 / &
           (reference_henry * exp(-coefficient * &
           (1._r8 / temperature - 1._r8 / kh_tbase)))
    end function henryEquilibriumConcentration

    pure real(r8) function henryDimensionlessSolubility(temperature, gas_index) result(value)
      real(r8), intent(in) :: temperature
      integer, intent(in) :: gas_index
      real(r8), parameter :: rgas_liter_atm = 0.0821_r8
      real(r8) :: coefficient, inverse_henry
      if (temperature <= 0._r8) then
         value = 0._r8
         return
      end if
      if (gas_index <= microbe_gas_co2) then
         coefficient = c_h_inv(gas_index)
         inverse_henry = kh_theta(gas_index) * exp(-coefficient * &
              (1._r8 / temperature - 1._r8 / kh_tbase))
      else
         coefficient = 500._r8
         inverse_henry = 1282.1_r8 * exp(-coefficient * &
              (1._r8 / temperature - 1._r8 / kh_tbase))
      end if
      value = max(0._r8, temperature * rgas_liter_atm / inverse_henry)
    end function henryDimensionlessSolubility

    pure real(r8) function atmosphericGasConcentration(pressure_pa, temperature) result(value)
      real(r8), intent(in) :: pressure_pa, temperature
      if (pressure_pa <= 0._r8 .or. temperature <= 0._r8) then
         value = 0._r8
      else
         ! rgas is J K-1 kmol-1; convert the result from kmol to mol.
         value = 1000._r8 * pressure_pa / (rgas * temperature)
      end if
    end function atmosphericGasConcentration

    real(r8) function surfaceConductance(column, gas_index, saturated_partition, &
         atmospheric_conductance, top_diffusivity) result(value)
      integer, intent(in) :: column, gas_index
      logical, intent(in) :: saturated_partition
      real(r8), intent(in) :: atmospheric_conductance, top_diffusivity
      real(r8) :: resistance, snow_diffusivity, pond_diffusivity
      real(r8) :: air_fraction, water_fraction, ice_fraction, fluid_fraction
      real(r8) :: pond_depth, henry_solubility
      integer :: snow_layer
      value = 0._r8
      if (atmospheric_conductance <= 0._r8 .or. top_diffusivity <= 0._r8) return
      resistance = 1._r8 / atmospheric_conductance + &
           0.5_r8 * col_pp%dz(column,1) / top_diffusivity

      do snow_layer = col_pp%snl(column) + 1, 0
         if (col_pp%dz(column,snow_layer) <= 0._r8) cycle
         ice_fraction = max(0._r8, col_ws%h2osoi_ice(column,snow_layer)) / &
              (denice * col_pp%dz(column,snow_layer))
         water_fraction = max(0._r8, col_ws%h2osoi_liq(column,snow_layer)) / &
              (denh2o * col_pp%dz(column,snow_layer))
         air_fraction = max(0._r8, 1._r8 - ice_fraction - water_fraction)
         fluid_fraction = air_fraction + water_fraction
         if (air_fraction > 0.05_r8 .and. fluid_fraction > tiny(1._r8)) then
            snow_diffusivity = referenceGasDiffusivity(gas_index, &
                 col_es%t_soisno(column,snow_layer)) * &
                 (air_fraction / fluid_fraction)**(10._r8/3._r8) / fluid_fraction**2 * &
                 CH4ParamsInst%scale_factor_gasdiff
         else
            henry_solubility = henryDimensionlessSolubility( &
                 col_es%t_soisno(column,snow_layer), gas_index)
            snow_diffusivity = microbeMethaneEffectiveAqueousDiffusivity( &
                 referenceAqueousDiffusivity(gas_index), &
                 col_es%t_soisno(column,snow_layer), &
                 water_fraction**CH4ParamsInst%satpow, MicrobeMethaneParamsInst) * &
                 CH4ParamsInst%scale_factor_liqdiff * henry_solubility
         end if
         if (snow_diffusivity <= 0._r8) return
         resistance = resistance + col_pp%dz(column,snow_layer) / snow_diffusivity
      end do

      if (saturated_partition .and. col_ws%frac_h2osfc(column) > 0._r8 .and. &
           col_ws%h2osfc(column) > 0._r8) then
         if (col_es%t_h2osfc(column) < MicrobeMethaneParamsInst%transport_thaw_threshold) return
         pond_depth = col_ws%h2osfc(column) / denh2o / col_ws%frac_h2osfc(column)
         henry_solubility = henryDimensionlessSolubility( &
              col_es%t_h2osfc(column), gas_index)
         pond_diffusivity = microbeMethaneEffectiveAqueousDiffusivity( &
              referenceAqueousDiffusivity(gas_index), col_es%t_h2osfc(column), &
              1._r8, MicrobeMethaneParamsInst) * &
              CH4ParamsInst%scale_factor_liqdiff * henry_solubility
         if (pond_diffusivity <= 0._r8) return
         resistance = resistance + pond_depth / pond_diffusivity
      end if
      if (resistance > 0._r8) value = 1._r8 / resistance
    end function surfaceConductance

    pure real(r8) function referenceGasDiffusivity(gas_index, temperature) result(value)
      integer, intent(in) :: gas_index
      real(r8), intent(in) :: temperature
      real(r8) :: temperature_c
      temperature_c = temperature - tfrz
      if (gas_index <= microbe_gas_co2) then
         value = max(0._r8, d_con_g(gas_index,1) + &
              d_con_g(gas_index,2) * temperature_c) * 1.e-4_r8
      else
         ! H2 molecular diffusivity in air near standard temperature (m2 s-1).
         value = 6.11e-5_r8 * (max(temperature, 1._r8) / 298._r8)**1.75_r8
      end if
    end function referenceGasDiffusivity

    pure real(r8) function clampUnitInterval(value) result(clamped)
      real(r8), intent(in) :: value
      clamped = min(1._r8, max(0._r8, value))
    end function clampUnitInterval

  end subroutine Advance

  subroutine UpdateAdditionalCarbon(this, bounds)
    ! Keep the derived storage term valid before the first balance check and
    ! immediately after restart. DOM is not included because ELM already
    ! carries it in the authoritative decomposition pools.
    use elm_varpar, only : nlevdecomp
    use ColumnType, only : col_pp
    use MicrobeMethaneReactionMod, only : microbe_methane_reaction_state_type
    use MicrobeMethaneStateUpdateMod, only : microbeMethaneColumnAdditionalCarbon

    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds
    type(microbe_methane_reaction_state_type) :: unsaturated_state(nlevdecomp)
    type(microbe_methane_reaction_state_type) :: saturated_state(nlevdecomp)
    real(r8) :: layer_thickness(nlevdecomp)
    integer :: c, j

    this%additional_carbon_col(bounds%begc:bounds%endc) = 0._r8
    do c = bounds%begc, bounds%endc
       if (.not. (col_pp%is_soil(c) .or. col_pp%is_crop(c))) cycle
       do j = 1, nlevdecomp
          unsaturated_state(j)%acetate_c = this%acetate_c_unsat_col(c,j)
          saturated_state(j)%acetate_c = this%acetate_c_sat_col(c,j)
          unsaturated_state(j)%acetate_methanogen_c = &
               this%acetate_methanogen_c_unsat_col(c,j)
          saturated_state(j)%acetate_methanogen_c = &
               this%acetate_methanogen_c_sat_col(c,j)
          unsaturated_state(j)%h2_methanogen_c = this%h2_methanogen_c_unsat_col(c,j)
          saturated_state(j)%h2_methanogen_c = this%h2_methanogen_c_sat_col(c,j)
          unsaturated_state(j)%aerobic_methanotroph_c = &
               this%aerobic_methanotroph_c_unsat_col(c,j)
          saturated_state(j)%aerobic_methanotroph_c = &
               this%aerobic_methanotroph_c_sat_col(c,j)
          unsaturated_state(j)%anaerobic_methanotroph_c = &
               this%anaerobic_methanotroph_c_unsat_col(c,j)
          saturated_state(j)%anaerobic_methanotroph_c = &
               this%anaerobic_methanotroph_c_sat_col(c,j)
          unsaturated_state(j)%conc_ch4 = this%conc_ch4_unsat_col(c,j)
          saturated_state(j)%conc_ch4 = this%conc_ch4_sat_col(c,j)
          unsaturated_state(j)%conc_o2 = this%conc_o2_unsat_col(c,j)
          saturated_state(j)%conc_o2 = this%conc_o2_sat_col(c,j)
          unsaturated_state(j)%conc_co2 = this%conc_co2_unsat_col(c,j)
          saturated_state(j)%conc_co2 = this%conc_co2_sat_col(c,j)
          unsaturated_state(j)%conc_h2 = this%conc_h2_unsat_col(c,j)
          saturated_state(j)%conc_h2 = this%conc_h2_sat_col(c,j)
          layer_thickness(j) = col_pp%dz(c,j)
       end do
       this%additional_carbon_col(c) = microbeMethaneColumnAdditionalCarbon( &
            this%sat_fraction_previous_col(c), unsaturated_state, saturated_state, &
            layer_thickness)
    end do
  end subroutine UpdateAdditionalCarbon

  subroutine Restart(this, bounds, ncid, flag)
    use ncdio_pio, only : ncd_double
    use restUtilMod, only : restartvar
    use elm_varctl, only : nsrest, nsrStartup
    use abortutils, only : endrun
    use shr_log_mod, only : errMsg => shr_log_errMsg

    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds
    type(file_desc_t), intent(inout) :: ncid
    character(len=*), intent(in) :: flag
    logical :: readvar

    if (.not. use_microbe_methane .or. use_legacy_ch4_with_microbe) return
    call restart_state('MM_ACETATE_C_UNSAT', 'gC/m^3', this%acetate_c_unsat_col)
    call restart_state('MM_ACETATE_C_SAT', 'gC/m^3', this%acetate_c_sat_col)
    call restart_state('MM_ACET_METH_C_UNSAT', 'gC/m^3', this%acetate_methanogen_c_unsat_col)
    call restart_state('MM_ACET_METH_C_SAT', 'gC/m^3', this%acetate_methanogen_c_sat_col)
    call restart_state('MM_H2_METH_C_UNSAT', 'gC/m^3', this%h2_methanogen_c_unsat_col)
    call restart_state('MM_H2_METH_C_SAT', 'gC/m^3', this%h2_methanogen_c_sat_col)
    call restart_state('MM_AER_METHANOTROPH_C_UNSAT', 'gC/m^3', this%aerobic_methanotroph_c_unsat_col)
    call restart_state('MM_AER_METHANOTROPH_C_SAT', 'gC/m^3', this%aerobic_methanotroph_c_sat_col)
    call restart_state('MM_ANAER_METHANOTROPH_C_UNSAT', 'gC/m^3', this%anaerobic_methanotroph_c_unsat_col)
    call restart_state('MM_ANAER_METHANOTROPH_C_SAT', 'gC/m^3', this%anaerobic_methanotroph_c_sat_col)
    call restart_state('MM_CONC_CH4_UNSAT', 'mol/m^3', this%conc_ch4_unsat_col)
    call restart_state('MM_CONC_CH4_SAT', 'mol/m^3', this%conc_ch4_sat_col)
    call restart_state('MM_CONC_O2_UNSAT', 'mol/m^3', this%conc_o2_unsat_col)
    call restart_state('MM_CONC_O2_SAT', 'mol/m^3', this%conc_o2_sat_col)
    call restart_state('MM_CONC_CO2_UNSAT', 'mol/m^3', this%conc_co2_unsat_col)
    call restart_state('MM_CONC_CO2_SAT', 'mol/m^3', this%conc_co2_sat_col)
    call restart_state('MM_CONC_H2_UNSAT', 'mol/m^3', this%conc_h2_unsat_col)
    call restart_state('MM_CONC_H2_SAT', 'mol/m^3', this%conc_h2_sat_col)
    call restartvar(ncid=ncid, flag=flag, varname='MM_SAT_FRACTION_PREVIOUS', xtype=ncd_double, &
         dim1name='column', long_name='previous revised-methane saturated-area fraction', units='1', &
         readvar=readvar, interpinic_flag='interp', data=this%sat_fraction_previous_col)
    call require_state_on_restart('MM_SAT_FRACTION_PREVIOUS', readvar)
    if (flag == 'read') call this%UpdateAdditionalCarbon(bounds)

  contains
    subroutine restart_state(name, units, field)
      character(len=*), intent(in) :: name, units
      real(r8), pointer, intent(inout) :: field(:,:)
      call restartvar(ncid=ncid, flag=flag, varname=name, xtype=ncd_double, &
           dim1name='column', dim2name='levgrnd', switchdim=.true., &
           long_name='revised microbial methane prognostic state', units=units, &
           readvar=readvar, interpinic_flag='interp', data=field)
      call require_state_on_restart(name, readvar)
    end subroutine restart_state

    subroutine require_state_on_restart(name, was_read)
      character(len=*), intent(in) :: name
      logical, intent(in) :: was_read

      if (flag == 'read' .and. nsrest /= nsrStartup .and. .not. was_read) then
         call endrun(msg=' ERROR: revised methane restart is missing required state: '//trim(name)//&
              '. A legacy restart cannot be continued or branched as revised methane.'//&
              errMsg(__FILE__, __LINE__))
      end if
    end subroutine require_state_on_restart
  end subroutine Restart

end module MicrobeMethaneMod
