module MicrobeMethaneStateUpdateMod

  ! Conservative state-update boundary for the revised microbial methane
  ! backend. The routines in this module are pure: they stage a complete
  ! update in caller-owned temporaries and never mutate ELM state directly.
  ! The ELM adapter validates the returned transaction before committing it.

  use shr_kind_mod, only : r8 => shr_kind_r8
  use elm_varcon, only : catomw
  use MicrobeMethaneParamsMod, only : MicrobeMethaneParamsType
  use MicrobeMethaneReactionMod, only : microbe_methane_reaction_state_type
  use MicrobeMethaneReactionMod, only : microbe_methane_reaction_environment_type
  use MicrobeMethaneReactionMod, only : microbe_methane_reaction_rates_type
  use MicrobeMethaneReactionMod, only : microbe_methane_reaction_tendencies_type
  use MicrobeMethaneReactionMod, only : computeMicrobeMethaneReactionTendencies
  use MicrobeMethaneReactionMod, only : assembleMicrobeMethaneReactionTendencies
  use MicrobeGasTransportMod, only : microbe_gas_ch4, microbe_gas_co2
  use MicrobeGasTransportMod, only : microbe_gas_count
  use MicrobeGasTransportMod, only : computeMicrobeGasTransport
  use MicrobeGasTransportMod, only : computeMicrobeVerticalDiffusion

  implicit none
  private
  save

  real(r8), parameter :: state_tolerance = 1.e-12_r8

  type, public :: microbe_methane_nutrient_tendencies_type
     ! Organic tendencies apply to the authoritative DOM pool. Equal and
     ! opposite mineral tendencies keep represented N and P closed.
     real(r8) :: dom_n = 0._r8
     real(r8) :: dom_p = 0._r8
     real(r8) :: mineral_n = 0._r8
     real(r8) :: mineral_p = 0._r8
  end type microbe_methane_nutrient_tendencies_type

  type, public :: microbe_methane_reaction_transaction_type
     type(microbe_methane_reaction_state_type) :: unsaturated_state
     type(microbe_methane_reaction_state_type) :: saturated_state
     type(microbe_methane_reaction_rates_type) :: unsaturated_rates
     type(microbe_methane_reaction_rates_type) :: saturated_rates
     type(microbe_methane_reaction_tendencies_type) :: unsaturated_tendencies
     type(microbe_methane_reaction_tendencies_type) :: saturated_tendencies
     type(microbe_methane_nutrient_tendencies_type) :: nutrient_tendencies
     real(r8) :: dom_c = 0._r8
     real(r8) :: dom_n = 0._r8
     real(r8) :: dom_p = 0._r8
     real(r8) :: mineral_n = 0._r8
     real(r8) :: mineral_p = 0._r8
     ! Per-step residuals in g C|N|P m-3. Reactions are a closed transaction.
     real(r8) :: carbon_residual = 0._r8
     real(r8) :: nitrogen_residual = 0._r8
     real(r8) :: phosphorus_residual = 0._r8
     logical :: valid = .false.
  end type microbe_methane_reaction_transaction_type

  public :: advanceMicrobeMethaneReactionLayer
  public :: advanceMicrobeMethaneGasTransport
  public :: advanceMicrobeMethaneAcetateTransport
  public :: microbeMethaneAdditionalCarbonDensity
  public :: microbeMethaneColumnAdditionalCarbon
  public :: microbeMethaneSurfaceCarbonFlux
  public :: microbeMethaneCH4SurfaceFluxKgC
  public :: microbeMethaneCO2Correction

contains

  pure subroutine advanceMicrobeMethaneReactionLayer(dom_c, dom_n, dom_p, mineral_n, &
       mineral_p, saturated_fraction, unsaturated_state, saturated_state, &
       unsaturated_environment, saturated_environment, parameters, cn_dom, cp_dom, &
       dt, transaction)
    real(r8), intent(in) :: dom_c, dom_n, dom_p
    real(r8), intent(in) :: mineral_n, mineral_p
    real(r8), intent(in) :: saturated_fraction
    type(microbe_methane_reaction_state_type), intent(in) :: unsaturated_state
    type(microbe_methane_reaction_state_type), intent(in) :: saturated_state
    type(microbe_methane_reaction_environment_type), intent(in) :: unsaturated_environment
    type(microbe_methane_reaction_environment_type), intent(in) :: saturated_environment
    type(MicrobeMethaneParamsType), intent(in) :: parameters
    real(r8), intent(in) :: cn_dom, cp_dom, dt
    type(microbe_methane_reaction_transaction_type), intent(out) :: transaction
    type(microbe_methane_reaction_state_type) :: unsaturated_work, saturated_work
    real(r8) :: fraction, available_dom_c, dom_c_tendency
    real(r8) :: initial_carbon, final_carbon
    real(r8) :: initial_nitrogen, final_nitrogen
    real(r8) :: initial_phosphorus, final_phosphorus

    transaction = microbe_methane_reaction_transaction_type()
    if (dt <= 0._r8 .or. cn_dom <= 0._r8 .or. cp_dom <= 0._r8) return
    fraction = clampUnitInterval(saturated_fraction)

    ! DOM is one bulk ELM pool. Restrict the carbon presented to both
    ! partition kernels to the amount carrying the selected DOM N and P.
    ! Because the two accepted rates are area weighted below, each may safely
    ! limit against this same concentration without double consuming it.
    if (cn_dom > 0._r8 .and. cp_dom > 0._r8) then
       available_dom_c = min(max(0._r8, dom_c), &
            max(0._r8, dom_n) * cn_dom, max(0._r8, dom_p) * cp_dom)
    else
       available_dom_c = 0._r8
    end if
    unsaturated_work = unsaturated_state
    saturated_work = saturated_state
    unsaturated_work%dom_c = available_dom_c
    saturated_work%dom_c = available_dom_c

    call computeMicrobeMethaneReactionTendencies(unsaturated_work, &
         unsaturated_environment, parameters, dt, transaction%unsaturated_rates, &
         transaction%unsaturated_tendencies)
    call computeMicrobeMethaneReactionTendencies(saturated_work, &
         saturated_environment, parameters, dt, transaction%saturated_rates, &
         transaction%saturated_tendencies)

    ! Guild mortality returns carbon-only biomass to DOM. If that produces a
    ! net DOM gain, immobilize N and P and proportionally limit mortality to
    ! the available mineral inventories before assembling final tendencies.
    call limitMortalityForDOMNutrients(fraction, mineral_n, mineral_p, cn_dom, &
         cp_dom, dt, transaction%unsaturated_rates, transaction%saturated_rates)
    call assembleMicrobeMethaneReactionTendencies(parameters, &
         transaction%unsaturated_rates, transaction%unsaturated_tendencies)
    call assembleMicrobeMethaneReactionTendencies(parameters, &
         transaction%saturated_rates, transaction%saturated_tendencies)

    dom_c_tendency = (1._r8 - fraction) * transaction%unsaturated_tendencies%dom_c + &
         fraction * transaction%saturated_tendencies%dom_c
    transaction%nutrient_tendencies%dom_n = dom_c_tendency / cn_dom
    transaction%nutrient_tendencies%dom_p = dom_c_tendency / cp_dom
    transaction%nutrient_tendencies%mineral_n = -transaction%nutrient_tendencies%dom_n
    transaction%nutrient_tendencies%mineral_p = -transaction%nutrient_tendencies%dom_p

    transaction%dom_c = dom_c + dt * dom_c_tendency
    transaction%dom_n = dom_n + dt * transaction%nutrient_tendencies%dom_n
    transaction%dom_p = dom_p + dt * transaction%nutrient_tendencies%dom_p
    transaction%mineral_n = mineral_n + dt * transaction%nutrient_tendencies%mineral_n
    transaction%mineral_p = mineral_p + dt * transaction%nutrient_tendencies%mineral_p

    call applyReactionTendencies(unsaturated_state, &
         transaction%unsaturated_tendencies, transaction%dom_c, dt, &
         transaction%unsaturated_state)
    call applyReactionTendencies(saturated_state, &
         transaction%saturated_tendencies, transaction%dom_c, dt, &
         transaction%saturated_state)

    initial_carbon = max(0._r8, dom_c) + &
         (1._r8 - fraction) * microbeMethaneAdditionalCarbonDensity(unsaturated_state) + &
         fraction * microbeMethaneAdditionalCarbonDensity(saturated_state)
    final_carbon = transaction%dom_c + &
         (1._r8 - fraction) * microbeMethaneAdditionalCarbonDensity(transaction%unsaturated_state) + &
         fraction * microbeMethaneAdditionalCarbonDensity(transaction%saturated_state)
    transaction%carbon_residual = final_carbon - initial_carbon
    initial_nitrogen = max(0._r8, dom_n) + max(0._r8, mineral_n)
    final_nitrogen = transaction%dom_n + transaction%mineral_n
    initial_phosphorus = max(0._r8, dom_p) + max(0._r8, mineral_p)
    final_phosphorus = transaction%dom_p + transaction%mineral_p
    transaction%nitrogen_residual = final_nitrogen - initial_nitrogen
    transaction%phosphorus_residual = final_phosphorus - initial_phosphorus
    transaction%valid = reactionStateIsNonnegative(transaction%unsaturated_state, &
         parameters%mfg_biomass_min) .and. &
         reactionStateIsNonnegative(transaction%saturated_state, parameters%mfg_biomass_min) .and. &
         transaction%dom_c >= -state_tolerance .and. transaction%dom_n >= -state_tolerance .and. &
         transaction%dom_p >= -state_tolerance .and. transaction%mineral_n >= -state_tolerance .and. &
         transaction%mineral_p >= -state_tolerance .and. &
         residualIsClosed(transaction%carbon_residual, initial_carbon, final_carbon) .and. &
         residualIsClosed(transaction%nitrogen_residual, initial_nitrogen, final_nitrogen) .and. &
         residualIsClosed(transaction%phosphorus_residual, initial_phosphorus, final_phosphorus)
  end subroutine advanceMicrobeMethaneReactionLayer

  pure subroutine advanceMicrobeMethaneGasTransport(state, layer_thickness, &
       effective_diffusivity, surface_equilibrium_concentration, surface_conductance, &
       aerenchyma_equilibrium_concentration, aerenchyma_exchange_rate, &
       ch4_ebullition_threshold, ch4_ebullition_activation, dt, updated_state, &
       interface_flux, aerenchyma_flux, ch4_ebullition_loss, surface_diffusive_flux, &
       surface_aerenchyma_flux, surface_ebullition_flux, total_surface_flux, &
       carbon_residual, valid)
    type(microbe_methane_reaction_state_type), intent(in) :: state(:)
    real(r8), intent(in) :: layer_thickness(size(state))
    real(r8), intent(in) :: effective_diffusivity(size(state),microbe_gas_count)
    real(r8), intent(in) :: surface_equilibrium_concentration(microbe_gas_count)
    real(r8), intent(in) :: surface_conductance(microbe_gas_count)
    real(r8), intent(in) :: aerenchyma_equilibrium_concentration(size(state),microbe_gas_count)
    real(r8), intent(in) :: aerenchyma_exchange_rate(size(state),microbe_gas_count)
    real(r8), intent(in) :: ch4_ebullition_threshold(size(state))
    real(r8), intent(in) :: ch4_ebullition_activation(size(state))
    real(r8), intent(in) :: dt
    type(microbe_methane_reaction_state_type), intent(out) :: updated_state(size(state))
    real(r8), intent(out) :: interface_flux(0:size(state),microbe_gas_count)
    real(r8), intent(out) :: aerenchyma_flux(size(state),microbe_gas_count)
    real(r8), intent(out) :: ch4_ebullition_loss(size(state))
    real(r8), intent(out) :: surface_diffusive_flux(microbe_gas_count)
    real(r8), intent(out) :: surface_aerenchyma_flux(microbe_gas_count)
    real(r8), intent(out) :: surface_ebullition_flux
    real(r8), intent(out) :: total_surface_flux(microbe_gas_count)
    real(r8), intent(out) :: carbon_residual
    logical, intent(out) :: valid
    real(r8) :: concentration(size(state),microbe_gas_count)
    real(r8) :: tendency(size(state),microbe_gas_count)
    real(r8) :: initial_carbon, final_carbon
    integer :: j

    do j = 1, size(state)
       call stateToGasVector(state(j), concentration(j,:))
    end do
    call computeMicrobeGasTransport(concentration, layer_thickness, effective_diffusivity, &
         surface_equilibrium_concentration, surface_conductance, &
         aerenchyma_equilibrium_concentration, aerenchyma_exchange_rate, &
         ch4_ebullition_threshold, ch4_ebullition_activation, dt, interface_flux, &
         aerenchyma_flux, ch4_ebullition_loss, tendency, surface_diffusive_flux, &
         surface_aerenchyma_flux, surface_ebullition_flux, total_surface_flux)

    updated_state = state
    do j = 1, size(state)
       updated_state(j)%conc_ch4 = state(j)%conc_ch4 + dt * tendency(j,1)
       updated_state(j)%conc_o2 = state(j)%conc_o2 + dt * tendency(j,2)
       updated_state(j)%conc_co2 = state(j)%conc_co2 + dt * tendency(j,3)
       updated_state(j)%conc_h2 = state(j)%conc_h2 + dt * tendency(j,4)
    end do

    initial_carbon = catomw * sum((concentration(:,microbe_gas_ch4) + &
         concentration(:,microbe_gas_co2)) * layer_thickness)
    final_carbon = catomw * sum(([(updated_state(j)%conc_ch4, j=1,size(state))] + &
         [(updated_state(j)%conc_co2, j=1,size(state))]) * layer_thickness)
    carbon_residual = final_carbon - initial_carbon + dt * &
         microbeMethaneSurfaceCarbonFlux(total_surface_flux)
    valid = all(concentration + dt * tendency >= -state_tolerance) .and. &
         residualIsClosed(carbon_residual, initial_carbon, final_carbon)
  end subroutine advanceMicrobeMethaneGasTransport

  pure subroutine advanceMicrobeMethaneAcetateTransport(state, layer_thickness, &
       effective_diffusivity, dt, updated_state, interface_flux, tendency, residual, valid)
    type(microbe_methane_reaction_state_type), intent(in) :: state(:)
    real(r8), intent(in) :: layer_thickness(size(state))
    real(r8), intent(in) :: effective_diffusivity(size(state))
    real(r8), intent(in) :: dt
    type(microbe_methane_reaction_state_type), intent(out) :: updated_state(size(state))
    real(r8), intent(out) :: interface_flux(0:size(state))
    real(r8), intent(out) :: tendency(size(state))
    real(r8), intent(out) :: residual
    logical, intent(out) :: valid
    real(r8) :: concentration(size(state)), unused_surface_flux
    integer :: j

    do j = 1, size(state)
       concentration(j) = state(j)%acetate_c
    end do
    ! Acetate has closed top and bottom boundaries. DOM retains the Phase 2
    ! decomposition-cascade transport path and is deliberately not handled here.
    call computeMicrobeVerticalDiffusion(concentration, layer_thickness, &
         effective_diffusivity, 0._r8, 0._r8, dt, interface_flux, tendency, &
         unused_surface_flux)
    updated_state = state
    do j = 1, size(state)
       updated_state(j)%acetate_c = state(j)%acetate_c + dt * tendency(j)
    end do
    residual = sum(([(updated_state(j)%acetate_c, j=1,size(state))] - &
         concentration) * layer_thickness)
    valid = all(concentration + dt * tendency >= -state_tolerance) .and. &
         residualIsClosed(residual, sum(concentration * layer_thickness), &
         sum((concentration + dt * tendency) * layer_thickness))
  end subroutine advanceMicrobeMethaneAcetateTransport

  pure real(r8) function microbeMethaneAdditionalCarbonDensity(state) result(carbon)
    type(microbe_methane_reaction_state_type), intent(in) :: state

    ! DOM is excluded because it is already part of ELM's decomp C state.
    carbon = state%acetate_c + state%acetate_methanogen_c + &
         state%h2_methanogen_c + state%aerobic_methanotroph_c + &
         state%anaerobic_methanotroph_c + catomw * &
         (state%conc_ch4 + state%conc_co2)
  end function microbeMethaneAdditionalCarbonDensity

  pure real(r8) function microbeMethaneColumnAdditionalCarbon(saturated_fraction, &
       unsaturated_state, saturated_state, layer_thickness) result(carbon)
    real(r8), intent(in) :: saturated_fraction
    type(microbe_methane_reaction_state_type), intent(in) :: unsaturated_state(:)
    type(microbe_methane_reaction_state_type), intent(in) :: saturated_state(size(unsaturated_state))
    real(r8), intent(in) :: layer_thickness(size(unsaturated_state))
    real(r8) :: fraction
    integer :: j

    fraction = clampUnitInterval(saturated_fraction)
    carbon = 0._r8
    do j = 1, size(unsaturated_state)
       carbon = carbon + layer_thickness(j) * ((1._r8 - fraction) * &
            microbeMethaneAdditionalCarbonDensity(unsaturated_state(j)) + fraction * &
            microbeMethaneAdditionalCarbonDensity(saturated_state(j)))
    end do
  end function microbeMethaneColumnAdditionalCarbon

  pure real(r8) function microbeMethaneSurfaceCarbonFlux(total_surface_flux) result(flux)
    real(r8), intent(in) :: total_surface_flux(microbe_gas_count)
    ! g C m-2 s-1, positive upward. O2 and H2 carry no carbon.
    flux = catomw * (total_surface_flux(microbe_gas_ch4) + &
         total_surface_flux(microbe_gas_co2))
  end function microbeMethaneSurfaceCarbonFlux

  pure real(r8) function microbeMethaneCH4SurfaceFluxKgC(total_surface_flux) result(flux)
    real(r8), intent(in) :: total_surface_flux(microbe_gas_count)
    ! Shared ELM methane exchange contract: kg C m-2 s-1, positive upward.
    flux = catomw * total_surface_flux(microbe_gas_ch4) / 1000._r8
  end function microbeMethaneCH4SurfaceFluxKgC

  pure real(r8) function microbeMethaneCO2Correction(total_surface_flux) result(flux)
    real(r8), intent(in) :: total_surface_flux(microbe_gas_count)
    ! NEE correction: g C m-2 s-1, positive upward.
    flux = catomw * total_surface_flux(microbe_gas_co2)
  end function microbeMethaneCO2Correction

  pure subroutine limitMortalityForDOMNutrients(fraction, mineral_n, mineral_p, &
       cn_dom, cp_dom, dt, unsaturated_rates, saturated_rates)
    real(r8), intent(in) :: fraction, mineral_n, mineral_p, cn_dom, cp_dom, dt
    type(microbe_methane_reaction_rates_type), intent(inout) :: unsaturated_rates
    type(microbe_methane_reaction_rates_type), intent(inout) :: saturated_rates
    real(r8) :: fermentation_loss, mortality_input, permitted_dom_gain, scale

    if (dt <= 0._r8 .or. cn_dom <= 0._r8 .or. cp_dom <= 0._r8) return
    fermentation_loss = 1.5_r8 * catomw * ((1._r8 - fraction) * &
         unsaturated_rates%dom_to_acetate_c + fraction * saturated_rates%dom_to_acetate_c)
    mortality_input = catomw * ((1._r8 - fraction) * mortalityRate(unsaturated_rates) + &
         fraction * mortalityRate(saturated_rates))
    permitted_dom_gain = min(max(0._r8, mineral_n) * cn_dom / dt, &
         max(0._r8, mineral_p) * cp_dom / dt)

    if (mortality_input > fermentation_loss + permitted_dom_gain .and. &
         mortality_input > 0._r8) then
       scale = clampUnitInterval((fermentation_loss + permitted_dom_gain) / mortality_input)
       call scaleMortality(unsaturated_rates, scale)
       call scaleMortality(saturated_rates, scale)
    end if
  end subroutine limitMortalityForDOMNutrients

  pure real(r8) function mortalityRate(rates) result(rate)
    type(microbe_methane_reaction_rates_type), intent(in) :: rates
    rate = rates%acetate_methanogen_mortality_c + rates%h2_methanogen_mortality_c + &
         rates%aerobic_methanotroph_mortality_c + &
         rates%anaerobic_methanotroph_mortality_c
  end function mortalityRate

  pure subroutine scaleMortality(rates, scale)
    type(microbe_methane_reaction_rates_type), intent(inout) :: rates
    real(r8), intent(in) :: scale
    rates%acetate_methanogen_mortality_c = &
         rates%acetate_methanogen_mortality_c * scale
    rates%h2_methanogen_mortality_c = rates%h2_methanogen_mortality_c * scale
    rates%aerobic_methanotroph_mortality_c = &
         rates%aerobic_methanotroph_mortality_c * scale
    rates%anaerobic_methanotroph_mortality_c = &
         rates%anaerobic_methanotroph_mortality_c * scale
  end subroutine scaleMortality

  pure subroutine applyReactionTendencies(state, tendencies, dom_c, dt, updated_state)
    type(microbe_methane_reaction_state_type), intent(in) :: state
    type(microbe_methane_reaction_tendencies_type), intent(in) :: tendencies
    real(r8), intent(in) :: dom_c, dt
    type(microbe_methane_reaction_state_type), intent(out) :: updated_state

    updated_state%dom_c = dom_c
    updated_state%acetate_c = state%acetate_c + dt * tendencies%acetate_c
    updated_state%acetate_methanogen_c = state%acetate_methanogen_c + &
         dt * tendencies%acetate_methanogen_c
    updated_state%h2_methanogen_c = state%h2_methanogen_c + &
         dt * tendencies%h2_methanogen_c
    updated_state%aerobic_methanotroph_c = state%aerobic_methanotroph_c + &
         dt * tendencies%aerobic_methanotroph_c
    updated_state%anaerobic_methanotroph_c = state%anaerobic_methanotroph_c + &
         dt * tendencies%anaerobic_methanotroph_c
    updated_state%conc_ch4 = state%conc_ch4 + dt * tendencies%conc_ch4
    updated_state%conc_o2 = state%conc_o2 + dt * tendencies%conc_o2
    updated_state%conc_co2 = state%conc_co2 + dt * tendencies%conc_co2
    updated_state%conc_h2 = state%conc_h2 + dt * tendencies%conc_h2
  end subroutine applyReactionTendencies

  pure subroutine stateToGasVector(state, concentration)
    type(microbe_methane_reaction_state_type), intent(in) :: state
    real(r8), intent(out) :: concentration(microbe_gas_count)
    concentration(1) = state%conc_ch4
    concentration(2) = state%conc_o2
    concentration(3) = state%conc_co2
    concentration(4) = state%conc_h2
  end subroutine stateToGasVector

  pure logical function reactionStateIsNonnegative(state, biomass_floor) result(valid)
    type(microbe_methane_reaction_state_type), intent(in) :: state
    real(r8), intent(in) :: biomass_floor
    valid = state%dom_c >= -state_tolerance .and. &
         state%acetate_c >= -state_tolerance .and. &
         state%acetate_methanogen_c >= biomass_floor - state_tolerance .and. &
         state%h2_methanogen_c >= biomass_floor - state_tolerance .and. &
         state%aerobic_methanotroph_c >= biomass_floor - state_tolerance .and. &
         state%anaerobic_methanotroph_c >= biomass_floor - state_tolerance .and. &
         state%conc_ch4 >= -state_tolerance .and. state%conc_o2 >= -state_tolerance .and. &
         state%conc_co2 >= -state_tolerance .and. state%conc_h2 >= -state_tolerance
  end function reactionStateIsNonnegative

  pure logical function residualIsClosed(residual, initial_inventory, final_inventory) result(closed)
    real(r8), intent(in) :: residual, initial_inventory, final_inventory
    real(r8) :: scale

    scale = max(1._r8, abs(initial_inventory), abs(final_inventory))
    closed = abs(residual) <= state_tolerance * scale
  end function residualIsClosed

  pure real(r8) function clampUnitInterval(value) result(clamped)
    real(r8), intent(in) :: value
    clamped = min(1._r8, max(0._r8, value))
  end function clampUnitInterval

end module MicrobeMethaneStateUpdateMod
