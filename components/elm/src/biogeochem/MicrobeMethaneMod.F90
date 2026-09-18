module MicrobeMethaneMod

  ! State and lifecycle owner for the revised microbial methane backend.
  ! The ELM driver selects this backend instead of CH4Mod when enabled.
  ! DOM, bacteria, and fungi remain authoritative decomposition pools.

  use shr_kind_mod, only : r8 => shr_kind_r8
  use decompMod, only : bounds_type
  use elm_varctl, only : use_microbe_methane
  use pio, only : file_desc_t

  implicit none
  private
  save

  type, public :: microbe_methane_type
     ! Carbon-bearing state is g C m-3 soil; gas state is mol gas m-3.
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

     ! Future conservative partition remapping uses this prior-step fraction.
     real(r8), pointer :: sat_fraction_previous_col(:) => null()

     ! Adapter-owned accounting fields. DOM is deliberately excluded from
     ! additional_carbon_col because it is already in ELM's decomp pools.
     real(r8), pointer :: additional_carbon_col(:) => null()
     real(r8), pointer :: surface_carbon_flux_col(:) => null()
     real(r8), pointer :: surface_ch4_flux_col(:) => null()
     real(r8), pointer :: surface_co2_flux_col(:) => null()
     real(r8), pointer :: ch4_production_col(:) => null()
     real(r8), pointer :: ch4_oxidation_col(:) => null()
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

    if (.not. use_microbe_methane) return
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
    allocate(this%sat_fraction_previous_col(begc:endc)); this%sat_fraction_previous_col = nan
    allocate(this%additional_carbon_col(begc:endc)); this%additional_carbon_col = nan
    allocate(this%surface_carbon_flux_col(begc:endc)); this%surface_carbon_flux_col = nan
    allocate(this%surface_ch4_flux_col(begc:endc)); this%surface_ch4_flux_col = nan
    allocate(this%surface_co2_flux_col(begc:endc)); this%surface_co2_flux_col = nan
    allocate(this%ch4_production_col(begc:endc)); this%ch4_production_col = nan
    allocate(this%ch4_oxidation_col(begc:endc)); this%ch4_oxidation_col = nan
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
    this%sat_fraction_previous_col = 0._r8
    this%additional_carbon_col = 0._r8
    this%surface_carbon_flux_col = 0._r8
    this%surface_ch4_flux_col = 0._r8
    this%surface_co2_flux_col = 0._r8
    this%ch4_production_col = 0._r8
    this%ch4_oxidation_col = 0._r8

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
          this%sat_fraction_previous_col(c) = 0._r8
          this%additional_carbon_col(c) = 0._r8
          this%surface_carbon_flux_col(c) = 0._r8
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
         'dissolved methane concentration in unsaturated subarea', this%conc_ch4_unsat_col)
    call add_state('MM_CONC_CH4_SAT', 'mol/m^3', &
         'dissolved methane concentration in saturated subarea', this%conc_ch4_sat_col)
    call add_state('MM_CONC_O2_UNSAT', 'mol/m^3', &
         'dissolved oxygen concentration in unsaturated subarea', this%conc_o2_unsat_col)
    call add_state('MM_CONC_O2_SAT', 'mol/m^3', &
         'dissolved oxygen concentration in saturated subarea', this%conc_o2_sat_col)
    call add_state('MM_CONC_CO2_UNSAT', 'mol/m^3', &
         'dissolved carbon dioxide concentration in unsaturated subarea', this%conc_co2_unsat_col)
    call add_state('MM_CONC_CO2_SAT', 'mol/m^3', &
         'dissolved carbon dioxide concentration in saturated subarea', this%conc_co2_sat_col)
    call add_state('MM_CONC_H2_UNSAT', 'mol/m^3', &
         'dissolved hydrogen concentration in unsaturated subarea', this%conc_h2_unsat_col)
    call add_state('MM_CONC_H2_SAT', 'mol/m^3', &
         'dissolved hydrogen concentration in saturated subarea', this%conc_h2_sat_col)
    call hist_addfld1d(fname='MM_ADDITIONAL_C', units='gC/m^2', avgflag='A', &
         long_name='revised methane carbon storage excluding authoritative DOM', &
         ptr_col=this%additional_carbon_col, default='inactive')
    call hist_addfld1d(fname='MM_SURFACE_C_FLUX', units='gC/m^2/s', avgflag='A', &
         long_name='net revised methane CH4 plus CO2 carbon flux; positive to atmosphere', &
         ptr_col=this%surface_carbon_flux_col, default='inactive')
    call hist_addfld1d(fname='MM_SURFACE_CH4_FLUX', units='kgC/m^2/s', avgflag='A', &
         long_name='net revised methane CH4 carbon flux; positive to atmosphere', &
         ptr_col=this%surface_ch4_flux_col, default='inactive')
    call hist_addfld1d(fname='MM_SURFACE_CO2_FLUX', units='gC/m^2/s', avgflag='A', &
         long_name='net revised methane CO2 carbon flux; positive to atmosphere', &
         ptr_col=this%surface_co2_flux_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_PROD', units='gC/m^2/s', avgflag='A', &
         long_name='gross revised methane production integrated over the soil column', &
         ptr_col=this%ch4_production_col, default='inactive')
    call hist_addfld1d(fname='MM_CH4_OXID', units='gC/m^2/s', avgflag='A', &
         long_name='gross revised methane oxidation integrated over the soil column', &
         ptr_col=this%ch4_oxidation_col, default='inactive')

  contains
    subroutine add_state(name, units, long_name, field)
      character(len=*), intent(in) :: name, units, long_name
      real(r8), pointer, intent(inout) :: field(:,:)
      call hist_addfld_decomp(fname=name, units=units, type2d='levdcmp', &
           avgflag='A', long_name=long_name, ptr_col=field, default='inactive')
    end subroutine add_state
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

    if (.not. use_microbe_methane) return

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

  subroutine Advance(this, bounds, num_soilc, filter_soilc, dt, atm2lnd_vars, &
       col_es, col_ws, chemstate_vars, soilstate_vars, soilhydrology_vars, &
       ground_conductance_patch, ch4_vars, col_cs, col_cf, col_ns, col_nf, col_ps)
    ! Mutable ELM boundary for the pure revised-methane kernels. Every update
    ! for one column is staged locally, checked, and then committed once.
    use abortutils, only : endrun
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use elm_varctl, only : iulog
    use elm_varpar, only : nlevdecomp, i_dom
    use elm_varcon, only : denh2o, denice, tfrz, d_con_w, d_con_g, catomw
    use elm_varcon, only : c_h_inv, kh_theta, kh_tbase
    use ColumnType, only : col_pp
    use ColumnDataType, only : column_energy_state, column_water_state
    use ColumnDataType, only : column_carbon_state, column_carbon_flux
    use ColumnDataType, only : column_nitrogen_state, column_nitrogen_flux
    use ColumnDataType, only : column_phosphorus_state
    use atm2lndType, only : atm2lnd_type
    use ChemStateType, only : chemstate_type
    use SoilStateType, only : soilstate_type
    use SoilHydrologyType, only : soilhydrology_type
    use CH4Mod, only : ch4_type
    use subgridAveMod, only : p2c
    use MicrobeDecompMod, only : MicrobeDecompParamsInst
    use MicrobeMethaneParamsMod, only : MicrobeMethaneParamsInst
    use MicrobeMethaneReactionMod, only : microbe_methane_reaction_state_type
    use MicrobeMethaneReactionMod, only : microbe_methane_reaction_environment_type
    use MicrobeMethaneStateUpdateMod, only : microbe_methane_reaction_transaction_type
    use MicrobeMethaneStateUpdateMod, only : advanceMicrobeMethaneReactionLayer
    use MicrobeMethaneStateUpdateMod, only : advanceMicrobeMethaneAcetateTransport
    use MicrobeMethaneStateUpdateMod, only : advanceMicrobeMethaneGasTransport
    use MicrobeMethaneStateUpdateMod, only : microbeMethaneColumnAdditionalCarbon
    use MicrobeMethaneStateUpdateMod, only : microbeMethaneSurfaceCarbonFlux
    use MicrobeMethaneStateUpdateMod, only : microbeMethaneCH4SurfaceFluxKgC
    use MicrobeMethaneStateUpdateMod, only : microbeMethaneCO2Correction
    use MicrobeGasTransportMod, only : repartitionMicrobeMethaneScalar
    use MicrobeGasTransportMod, only : microbeMethaneEffectiveAqueousDiffusivity
    use MicrobeGasTransportMod, only : microbe_gas_ch4, microbe_gas_o2
    use MicrobeGasTransportMod, only : microbe_gas_co2, microbe_gas_h2
    use MicrobeGasTransportMod, only : microbe_gas_count
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc
    integer, intent(in) :: filter_soilc(:)
    real(r8), intent(in) :: dt
    type(atm2lnd_type), intent(in) :: atm2lnd_vars
    type(column_energy_state), intent(in) :: col_es
    type(column_water_state), intent(in) :: col_ws
    type(chemstate_type), intent(in) :: chemstate_vars
    type(soilstate_type), intent(in) :: soilstate_vars
    type(soilhydrology_type), intent(in) :: soilhydrology_vars
    real(r8), intent(in) :: ground_conductance_patch(bounds%begp:)
    type(ch4_type), intent(inout) :: ch4_vars
    type(column_carbon_state), intent(inout) :: col_cs
    type(column_carbon_flux), intent(in) :: col_cf
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
    real(r8) :: mineral_n(nlevdecomp), mineral_p(nlevdecomp)
    real(r8) :: layer_thickness(nlevdecomp), layer_depth(nlevdecomp)
    real(r8) :: unsaturated_diffusivity(nlevdecomp,microbe_gas_count)
    real(r8) :: saturated_diffusivity(nlevdecomp,microbe_gas_count)
    real(r8) :: unsaturated_acetate_diffusivity(nlevdecomp)
    real(r8) :: saturated_acetate_diffusivity(nlevdecomp)
    real(r8) :: surface_equilibrium(microbe_gas_count)
    real(r8) :: aerenchyma_equilibrium(nlevdecomp,microbe_gas_count)
    real(r8) :: aerenchyma_exchange_rate(nlevdecomp,microbe_gas_count)
    real(r8) :: unsaturated_surface_conductance(microbe_gas_count)
    real(r8) :: saturated_surface_conductance(microbe_gas_count)
    real(r8) :: ch4_ebullition_threshold(nlevdecomp)
    real(r8) :: unsaturated_ebullition_activation(nlevdecomp)
    real(r8) :: saturated_ebullition_activation(nlevdecomp)
    real(r8) :: interface_flux(0:nlevdecomp,microbe_gas_count)
    real(r8) :: aerenchyma_flux(nlevdecomp,microbe_gas_count)
    real(r8) :: ch4_ebullition_loss(nlevdecomp)
    real(r8) :: surface_diffusive_flux(microbe_gas_count)
    real(r8) :: surface_aerenchyma_flux(microbe_gas_count)
    real(r8) :: unsaturated_surface_flux(microbe_gas_count)
    real(r8) :: saturated_surface_flux(microbe_gas_count)
    real(r8) :: surface_ebullition_flux
    real(r8) :: acetate_interface_flux(0:nlevdecomp)
    real(r8) :: acetate_tendency(nlevdecomp)
    real(r8) :: ground_conductance_col(bounds%begc:bounds%endc)
    real(r8) :: root_fraction_col(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: partial_pressure(microbe_gas_count)
    real(r8) :: liquid_saturation, thawed_fraction, moisture_scalar
    real(r8) :: saturation_scalar, fraction, old_fraction, depth_scale
    real(r8) :: bulk_surface_flux(microbe_gas_count)
    real(r8) :: unused_surface_diffusive_flux(microbe_gas_count)
    real(r8) :: unused_surface_aerenchyma_flux(microbe_gas_count)
    real(r8) :: unused_surface_ebullition_flux
    real(r8) :: unsaturated_gas_residual, saturated_gas_residual
    real(r8) :: unsaturated_acetate_residual, saturated_acetate_residual
    logical :: column_valid
    logical :: unsaturated_acetate_valid, saturated_acetate_valid
    logical :: unsaturated_gas_valid, saturated_gas_valid
    integer :: c, fc, g, gas, j
    character(len=512) :: message

    if (.not. use_microbe_methane) return
    if (dt <= 0._r8) then
       call endrun(msg=' ERROR: revised methane adapter requires a positive timestep'//&
            errMsg(__FILE__, __LINE__))
    end if

    this%surface_carbon_flux_col(bounds%begc:bounds%endc) = 0._r8
    this%surface_ch4_flux_col(bounds%begc:bounds%endc) = 0._r8
    this%surface_co2_flux_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_production_col(bounds%begc:bounds%endc) = 0._r8
    this%ch4_oxidation_col(bounds%begc:bounds%endc) = 0._r8

    call p2c(bounds, num_soilc, filter_soilc, &
         ground_conductance_patch(bounds%begp:bounds%endp), &
         ground_conductance_col(bounds%begc:bounds%endc))
    ! Legacy CH4 computes this column average inside its solver. Revised mode
    ! bypasses that solver, so it must establish its own plant-transport
    ! input directly from the authoritative patch root profile.
    call p2c(bounds, nlevdecomp, &
         soilstate_vars%rootfr_patch(bounds%begp:bounds%endp,1:nlevdecomp), &
         root_fraction_col(bounds%begc:bounds%endc,1:nlevdecomp), 0)

    do fc = 1, num_soilc
       c = filter_soilc(fc)
       g = col_pp%gridcell(c)
       fraction = clampUnitInterval(max(soilhydrology_vars%fsat_col(c), &
            col_ws%frac_h2osfc(c)))
       old_fraction = clampUnitInterval(this%sat_fraction_previous_col(c))

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

       call gatherColumnState(c, unsaturated_state, saturated_state)
       call repartitionColumnState(old_fraction, fraction, unsaturated_state, saturated_state)
       call consumeELMAerobicOxygen(c, fraction, unsaturated_state, saturated_state)
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
          moisture_scalar = soilMoistureResponse(soilstate_vars%soilpsi_col(c,j), &
               soilstate_vars%sucsat_col(c,j), &
               MicrobeMethaneParamsInst%soil_water_potential_min)
          saturation_scalar = clampUnitInterval(liquid_saturation / &
               max(MicrobeMethaneParamsInst%saturation_reaction_threshold, tiny(1._r8)))

          unsaturated_environment%soil_temperature = col_es%t_soisno(c,j)
          ! Native ELM allocates chemstate soil pH but does not currently
          ! populate it. External chemistry backends that do populate it are
          ! rejected by the revised-methane configuration gate. Use the named
          ! optimum as the explicit native-ELM fallback until Phase 4 adds a
          ! spatial soil-pH input; acetate feedback can still lower the
          ! effective pH inside the reaction kernel.
          unsaturated_environment%soil_ph = MicrobeMethaneParamsInst%ph_opt
          unsaturated_environment%dom_fermentation_scalar = &
               moisture_scalar * saturation_scalar
          unsaturated_environment%aerobic_acetate_oxidation_scalar = &
               moisture_scalar * (1._r8 - saturation_scalar)
          saturated_environment%soil_temperature = col_es%t_soisno(c,j)
          saturated_environment%soil_ph = MicrobeMethaneParamsInst%ph_opt
          saturated_environment%dom_fermentation_scalar = moisture_scalar
          saturated_environment%aerobic_acetate_oxidation_scalar = 0._r8

          call advanceMicrobeMethaneReactionLayer(dom_c(j), dom_n(j), dom_p(j), &
               mineral_n(j), mineral_p(j), fraction, unsaturated_work(j), &
               saturated_work(j), unsaturated_environment, saturated_environment, &
               MicrobeMethaneParamsInst, MicrobeDecompParamsInst%cn_dom, &
               MicrobeDecompParamsInst%cp_dom, dt, reaction(j))
          if (.not. reaction(j)%valid) column_valid = .false.
          dom_c(j) = reaction(j)%dom_c
          dom_n(j) = reaction(j)%dom_n
          dom_p(j) = reaction(j)%dom_p
          mineral_n(j) = reaction(j)%mineral_n
          mineral_p(j) = reaction(j)%mineral_p
          unsaturated_work(j) = reaction(j)%unsaturated_state
          saturated_work(j) = reaction(j)%saturated_state

          do gas = 1, microbe_gas_count
             unsaturated_diffusivity(j,gas) = &
                  microbeMethaneEffectiveAqueousDiffusivity( &
                  referenceAqueousDiffusivity(gas), col_es%t_soisno(c,j), &
                  liquid_saturation * thawed_fraction, MicrobeMethaneParamsInst)
             saturated_diffusivity(j,gas) = &
                  microbeMethaneEffectiveAqueousDiffusivity( &
                  referenceAqueousDiffusivity(gas), col_es%t_soisno(c,j), &
                  thawed_fraction, MicrobeMethaneParamsInst)
             aerenchyma_equilibrium(j,gas) = henryEquilibriumConcentration( &
                  partial_pressure(gas), col_es%t_soisno(c,j), gas)
          end do
          unsaturated_acetate_diffusivity(j) = &
               MicrobeMethaneParamsInst%dom_diffusivity * liquid_saturation * thawed_fraction
          saturated_acetate_diffusivity(j) = &
               MicrobeMethaneParamsInst%dom_diffusivity * thawed_fraction

          ! The legacy coefficient is dimensionally m s-1: division by a
          ! finite layer depth converts it to the s-1 rate required by the
          ! signed exchange kernel.
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
          aerenchyma_equilibrium(j,microbe_gas_ch4) = max( &
               aerenchyma_equilibrium(j,microbe_gas_ch4), &
               1.e-3_r8 * MicrobeMethaneParamsInst%ch4_transport_threshold)
          aerenchyma_equilibrium(j,microbe_gas_h2) = max( &
               aerenchyma_equilibrium(j,microbe_gas_h2), &
               1.e-3_r8 * MicrobeMethaneParamsInst%h2_plant_transport_threshold)
          ch4_ebullition_threshold(j) = &
               1.e-3_r8 * MicrobeMethaneParamsInst%ch4_transport_threshold
          unsaturated_ebullition_activation(j) = thawed_fraction * liquid_saturation * &
               exp(-layer_depth(j) / MicrobeMethaneParamsInst%ebullition_efold_depth)
          saturated_ebullition_activation(j) = thawed_fraction * &
               exp(-layer_depth(j) / MicrobeMethaneParamsInst%ebullition_efold_depth)
       end do

       do gas = 1, microbe_gas_count
          surface_equilibrium(gas) = henryEquilibriumConcentration( &
               partial_pressure(gas), col_es%t_soisno(c,1), gas)
          unsaturated_surface_conductance(gas) = surfaceConductance(c, gas, .false., &
               ground_conductance_col(c), unsaturated_diffusivity(1,gas))
          saturated_surface_conductance(gas) = surfaceConductance(c, gas, .true., &
               ground_conductance_col(c), saturated_diffusivity(1,gas))
       end do

       call advanceMicrobeMethaneAcetateTransport(unsaturated_work, layer_thickness, &
            unsaturated_acetate_diffusivity, dt, unsaturated_candidate, &
            acetate_interface_flux, acetate_tendency, unsaturated_acetate_residual, &
            unsaturated_acetate_valid)
       call advanceMicrobeMethaneAcetateTransport(saturated_work, layer_thickness, &
            saturated_acetate_diffusivity, dt, saturated_candidate, &
            acetate_interface_flux, acetate_tendency, saturated_acetate_residual, &
            saturated_acetate_valid)
       column_valid = column_valid .and. unsaturated_acetate_valid .and. &
            saturated_acetate_valid

       call advanceMicrobeMethaneGasTransport(unsaturated_candidate, layer_thickness, &
            unsaturated_diffusivity, surface_equilibrium, &
            unsaturated_surface_conductance, aerenchyma_equilibrium, &
            aerenchyma_exchange_rate, ch4_ebullition_threshold, &
            unsaturated_ebullition_activation, dt, unsaturated_work, interface_flux, &
            aerenchyma_flux, ch4_ebullition_loss, surface_diffusive_flux, &
            surface_aerenchyma_flux, surface_ebullition_flux, &
            unsaturated_surface_flux, unsaturated_gas_residual, unsaturated_gas_valid)
       call advanceMicrobeMethaneGasTransport(saturated_candidate, layer_thickness, &
            saturated_diffusivity, surface_equilibrium, saturated_surface_conductance, &
            aerenchyma_equilibrium, aerenchyma_exchange_rate, ch4_ebullition_threshold, &
            saturated_ebullition_activation, dt, saturated_work, interface_flux, &
            aerenchyma_flux, ch4_ebullition_loss, unused_surface_diffusive_flux, &
            unused_surface_aerenchyma_flux, unused_surface_ebullition_flux, &
            saturated_surface_flux, saturated_gas_residual, saturated_gas_valid)
       column_valid = column_valid .and. unsaturated_gas_valid .and. saturated_gas_valid

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
          write(iulog,*) 'revised methane unsaturated top state before gas transport: ', &
               unsaturated_candidate(1)%conc_ch4, unsaturated_candidate(1)%conc_o2, &
               unsaturated_candidate(1)%conc_co2, unsaturated_candidate(1)%conc_h2
          write(iulog,*) 'revised methane unsaturated top state after gas transport: ', &
               unsaturated_work(1)%conc_ch4, unsaturated_work(1)%conc_o2, &
               unsaturated_work(1)%conc_co2, unsaturated_work(1)%conc_h2
          write(message,'(a,i0,5(a,l1),4(a,es12.4))') &
               ' ERROR: revised methane column transaction failed validation for column ', c, &
               '; reaction=', all(reaction(:)%valid), &
               '; acetate_unsat=', unsaturated_acetate_valid, &
               '; acetate_sat=', saturated_acetate_valid, &
               '; gas_unsat=', unsaturated_gas_valid, &
               '; gas_sat=', saturated_gas_valid, &
               '; acetate_residual_unsat=', unsaturated_acetate_residual, &
               '; acetate_residual_sat=', saturated_acetate_residual, &
               '; gas_residual_unsat=', unsaturated_gas_residual, &
               '; gas_residual_sat=', saturated_gas_residual
          call endrun(msg=trim(message)//errMsg(__FILE__, __LINE__))
       end if

       bulk_surface_flux = (1._r8 - fraction) * unsaturated_surface_flux + &
            fraction * saturated_surface_flux
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
       do j = 1, nlevdecomp
          this%ch4_production_col(c) = this%ch4_production_col(c) + catomw * &
               layer_thickness(j) * ((1._r8 - fraction) * &
               (reaction(j)%unsaturated_rates%acetoclastic_methanogenesis_c + &
               reaction(j)%unsaturated_rates%hydrogenotrophic_methanogenesis_c) + &
               fraction * (reaction(j)%saturated_rates%acetoclastic_methanogenesis_c + &
               reaction(j)%saturated_rates%hydrogenotrophic_methanogenesis_c))
          this%ch4_oxidation_col(c) = this%ch4_oxidation_col(c) + catomw * &
               layer_thickness(j) * ((1._r8 - fraction) * &
               (reaction(j)%unsaturated_rates%aerobic_methane_oxidation_c + &
               reaction(j)%unsaturated_rates%anaerobic_methane_oxidation_c) + &
               fraction * (reaction(j)%saturated_rates%aerobic_methane_oxidation_c + &
               reaction(j)%saturated_rates%anaerobic_methane_oxidation_c))
       end do
       call syncLegacyOxygenBridge(c, fraction, unsaturated_work, saturated_work)
    end do

  contains

    subroutine reportReactionState(label, state)
      character(len=*), intent(in) :: label
      type(microbe_methane_reaction_state_type), intent(in) :: state

      write(iulog,*) trim(label)//' DOM/acetate/guild C: ', state%dom_c, &
           state%acetate_c, state%acetate_methanogen_c, state%h2_methanogen_c, &
           state%aerobic_methanotroph_c, state%anaerobic_methanotroph_c
      write(iulog,*) trim(label)//' CH4/O2/CO2/H2: ', state%conc_ch4, &
           state%conc_o2, state%conc_co2, state%conc_h2
    end subroutine reportReactionState

    subroutine consumeELMAerobicOxygen(column, saturated_fraction, unsaturated, saturated)
      ! Standard decomposition, root respiration, and allocation have all
      ! resolved before this adapter runs. Remove their actual aerobic demand,
      ! including nitrification, before revised methane reactions compete for
      ! the remaining prognostic O2 inventory.
      integer, intent(in) :: column
      real(r8), intent(in) :: saturated_fraction
      type(microbe_methane_reaction_state_type), intent(inout) :: unsaturated(nlevdecomp)
      type(microbe_methane_reaction_state_type), intent(inout) :: saturated(nlevdecomp)
      real(r8) :: bulk_inventory, oxygen_demand_rate, oxygen_consumed
      real(r8) :: remaining_fraction
      integer :: layer

      do layer = 1, nlevdecomp
         bulk_inventory = (1._r8 - saturated_fraction) * &
              max(0._r8, unsaturated(layer)%conc_o2) + saturated_fraction * &
              max(0._r8, saturated(layer)%conc_o2)
         oxygen_demand_rate = max(0._r8, col_cf%hr_vr(column,layer)) / catomw
         oxygen_demand_rate = oxygen_demand_rate + &
              max(0._r8, col_cf%rr_vr(column,layer)) / &
              (catomw * max(col_pp%dz(column,layer), tiny(1._r8)))
         oxygen_demand_rate = oxygen_demand_rate + &
              max(0._r8, col_nf%f_nit_vr(column,layer)) * (2._r8 / 14._r8)
         oxygen_consumed = oxygen_demand_rate * dt
         if (bulk_inventory > tiny(1._r8)) then
            remaining_fraction = max(0._r8, 1._r8 - oxygen_consumed / bulk_inventory)
            unsaturated(layer)%conc_o2 = unsaturated(layer)%conc_o2 * remaining_fraction
            saturated(layer)%conc_o2 = saturated(layer)%conc_o2 * remaining_fraction
         else
            unsaturated(layer)%conc_o2 = 0._r8
            saturated(layer)%conc_o2 = 0._r8
         end if
      end do
    end subroutine consumeELMAerobicOxygen

    subroutine syncLegacyOxygenBridge(column, saturated_fraction, unsaturated, saturated)
      ! Decomposition and nitrification still consume the CH4Mod oxygen
      ! interface.  In revised mode ch4_vars is only a compatibility carrier;
      ! the legacy CH4 solver is not executed.  Publish the revised O2 state
      ! and a lagged estimate of aerobic demand for the next ELM timestep.
      integer, intent(in) :: column
      real(r8), intent(in) :: saturated_fraction
      type(microbe_methane_reaction_state_type), intent(in) :: unsaturated(nlevdecomp)
      type(microbe_methane_reaction_state_type), intent(in) :: saturated(nlevdecomp)
      real(r8) :: potential_demand, stress_unsaturated, stress_saturated
      integer :: layer

      ch4_vars%finundated_col(column) = saturated_fraction
      do layer = 1, nlevdecomp
         potential_demand = max(0._r8, col_cf%phr_vr(column,layer)) / catomw
         if (col_cf%o_scalar(column,layer) > tiny(1._r8)) then
            potential_demand = potential_demand / col_cf%o_scalar(column,layer)
         end if
         potential_demand = potential_demand + &
              max(0._r8, col_cf%rr_vr(column,layer)) / &
              (catomw * max(layer_thickness(layer), tiny(1._r8)))
         potential_demand = potential_demand + &
              max(0._r8, col_nf%f_nit_vr(column,layer)) * 2._r8 / 14._r8

         if (potential_demand > tiny(1._r8)) then
            stress_unsaturated = min(unsaturated(layer)%conc_o2 / dt / potential_demand, 1._r8)
            stress_saturated = min(saturated(layer)%conc_o2 / dt / potential_demand, 1._r8)
         else
            stress_unsaturated = 1._r8
            stress_saturated = 1._r8
         end if

         ch4_vars%conc_o2_unsat_col(column,layer) = unsaturated(layer)%conc_o2
         ch4_vars%conc_o2_sat_col(column,layer) = saturated(layer)%conc_o2
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

    real(r8) function surfaceConductance(column, gas_index, saturated_partition, &
         atmospheric_conductance, top_diffusivity) result(value)
      integer, intent(in) :: column, gas_index
      logical, intent(in) :: saturated_partition
      real(r8), intent(in) :: atmospheric_conductance, top_diffusivity
      real(r8) :: resistance, snow_diffusivity, pond_diffusivity
      real(r8) :: air_fraction, water_fraction, ice_fraction, fluid_fraction
      real(r8) :: temperature_c, pond_depth
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
                 (air_fraction / fluid_fraction)**(10._r8/3._r8) / fluid_fraction**2
         else
            snow_diffusivity = microbeMethaneEffectiveAqueousDiffusivity( &
                 referenceAqueousDiffusivity(gas_index), &
                 col_es%t_soisno(column,snow_layer), water_fraction, &
                 MicrobeMethaneParamsInst)
         end if
         if (snow_diffusivity <= 0._r8) return
         resistance = resistance + col_pp%dz(column,snow_layer) / snow_diffusivity
      end do

      if (saturated_partition .and. col_ws%frac_h2osfc(column) > 0._r8 .and. &
           col_ws%h2osfc(column) > 0._r8) then
         if (col_es%t_h2osfc(column) < MicrobeMethaneParamsInst%transport_thaw_threshold) return
         pond_depth = col_ws%h2osfc(column) / denh2o / col_ws%frac_h2osfc(column)
         pond_diffusivity = microbeMethaneEffectiveAqueousDiffusivity( &
              referenceAqueousDiffusivity(gas_index), col_es%t_h2osfc(column), &
              1._r8, MicrobeMethaneParamsInst)
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

    if (.not. use_microbe_methane) return
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
