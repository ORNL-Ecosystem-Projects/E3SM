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
  ! The implicit aqueous solve can be mildly ill-conditioned when porewater
  ! advection and dispersion nearly cancel across very dry interfaces. Keep
  ! its inventory check distinct from the stricter chemistry/state check.
  real(r8), parameter :: aqueous_budget_tolerance = 1.e-10_r8

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
  public :: advanceMicrobeMethaneDOMRelaxation
  public :: advanceMicrobeAqueousTracerTransport
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
    ! Native RD nutrient updates can leave a depleted mineral pool slightly
    ! negative. Preserve that inherited deficit in the transaction ledger;
    ! the mortality limiter treats it as zero available nutrient.
    initial_nitrogen = max(0._r8, dom_n) + mineral_n
    final_nitrogen = transaction%dom_n + transaction%mineral_n
    initial_phosphorus = max(0._r8, dom_p) + mineral_p
    final_phosphorus = transaction%dom_p + transaction%mineral_p
    transaction%nitrogen_residual = final_nitrogen - initial_nitrogen
    transaction%phosphorus_residual = final_phosphorus - initial_phosphorus
    transaction%valid = reactionStateIsNonnegative(transaction%unsaturated_state, &
         parameters%mfg_biomass_min) .and. &
         reactionStateIsNonnegative(transaction%saturated_state, parameters%mfg_biomass_min) .and. &
         transaction%dom_c >= -state_tolerance .and. transaction%dom_n >= -state_tolerance .and. &
         transaction%dom_p >= -state_tolerance .and. &
         transaction%mineral_n >= min(0._r8, mineral_n) - state_tolerance .and. &
         transaction%mineral_p >= min(0._r8, mineral_p) - state_tolerance .and. &
         residualIsClosed(transaction%carbon_residual, initial_carbon, final_carbon) .and. &
         residualIsClosed(transaction%nitrogen_residual, initial_nitrogen, final_nitrogen) .and. &
         residualIsClosed(transaction%phosphorus_residual, initial_phosphorus, final_phosphorus)
  end subroutine advanceMicrobeMethaneReactionLayer

  pure subroutine advanceMicrobeMethaneGasTransport(state, layer_thickness, &
       effective_diffusivity, transport_capacity, surface_equilibrium_concentration, &
       surface_conductance, &
       aerenchyma_equilibrium_concentration, aerenchyma_exchange_rate, &
       aerenchyma_minimum_emission_concentration, aerenchyma_allows_influx, &
       ch4_ebullition_threshold, ch4_ebullition_activation, dt, updated_state, &
       interface_flux, aerenchyma_flux, ch4_ebullition_loss, surface_diffusive_flux, &
       surface_aerenchyma_flux, surface_ebullition_flux, total_surface_flux, &
       carbon_residual, valid)
    type(microbe_methane_reaction_state_type), intent(in) :: state(:)
    real(r8), intent(in) :: layer_thickness(size(state))
    real(r8), intent(in) :: effective_diffusivity(size(state),microbe_gas_count)
    real(r8), intent(in) :: transport_capacity(size(state),microbe_gas_count)
    real(r8), intent(in) :: surface_equilibrium_concentration(microbe_gas_count)
    real(r8), intent(in) :: surface_conductance(microbe_gas_count)
    real(r8), intent(in) :: aerenchyma_equilibrium_concentration(size(state),microbe_gas_count)
    real(r8), intent(in) :: aerenchyma_exchange_rate(size(state),microbe_gas_count)
    real(r8), intent(in) :: aerenchyma_minimum_emission_concentration(size(state),microbe_gas_count)
    logical, intent(in) :: aerenchyma_allows_influx(microbe_gas_count)
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
         transport_capacity, surface_equilibrium_concentration, surface_conductance, &
         aerenchyma_equilibrium_concentration, aerenchyma_exchange_rate, &
         aerenchyma_minimum_emission_concentration, aerenchyma_allows_influx, &
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

  pure subroutine advanceMicrobeAqueousTracerTransport(concentration, &
       layer_thickness, liquid_fraction, diffusion_conductivity, water_flux, &
       mobile_fraction, minimum_liquid_fraction, dt, updated_concentration, &
       advective_flux, diffusive_flux, tendency, boundary_export, residual, valid)
    ! Conservative backward-Euler transport of one bulk-soil solute inventory.
    ! concentration is mass per bulk-soil volume. Fluxes use the mobile
    ! porewater concentration mobile_fraction*concentration/liquid_fraction.
    ! water_flux is positive downward and has units m s-1. The external
    ! concentration is zero at both boundaries, so infiltration is solute-free
    ! and outward flow exports the donor-layer concentration.
    real(r8), intent(in) :: concentration(:)
    real(r8), intent(in) :: layer_thickness(size(concentration))
    real(r8), intent(in) :: liquid_fraction(size(concentration))
    real(r8), intent(in) :: diffusion_conductivity(size(concentration))
    real(r8), intent(in) :: water_flux(0:size(concentration))
    real(r8), intent(in) :: mobile_fraction, minimum_liquid_fraction, dt
    real(r8), intent(out) :: updated_concentration(size(concentration))
    real(r8), intent(out) :: advective_flux(0:size(concentration))
    real(r8), intent(out) :: diffusive_flux(0:size(concentration))
    real(r8), intent(out) :: tendency(size(concentration))
    real(r8), intent(out) :: boundary_export, residual
    logical, intent(out) :: valid
    real(r8) :: porewater_factor(size(concentration))
    real(r8) :: conductance(0:size(concentration))
    real(r8) :: lower(size(concentration)), diagonal(size(concentration))
    real(r8) :: upper(size(concentration)), rhs(size(concentration))
    real(r8) :: cprime(size(concentration)), dprime(size(concentration))
    real(r8) :: left_coefficient, right_coefficient, resistance
    real(r8) :: denominator, initial_inventory, final_inventory
    integer :: j, number_of_layers

    number_of_layers = size(concentration)
    updated_concentration = concentration
    advective_flux = 0._r8
    diffusive_flux = 0._r8
    tendency = 0._r8
    boundary_export = 0._r8
    residual = 0._r8
    valid = .false.
    if (number_of_layers == 0 .or. dt <= 0._r8) return
    if (any(layer_thickness <= 0._r8) .or. any(liquid_fraction < 0._r8) .or. &
         any(diffusion_conductivity < 0._r8) .or. any(concentration < -state_tolerance)) return
    if (mobile_fraction <= 0._r8 .or. mobile_fraction > 1._r8 .or. &
         minimum_liquid_fraction <= 0._r8) return

    porewater_factor = 0._r8
    do j = 1, number_of_layers
       if (liquid_fraction(j) >= minimum_liquid_fraction) then
          porewater_factor(j) = mobile_fraction / liquid_fraction(j)
       end if
    end do

    conductance = 0._r8
    do j = 1, number_of_layers - 1
       if (diffusion_conductivity(j) > 0._r8 .and. &
            diffusion_conductivity(j+1) > 0._r8 .and. &
            porewater_factor(j) > 0._r8 .and. porewater_factor(j+1) > 0._r8) then
          resistance = 0.5_r8 * layer_thickness(j) / diffusion_conductivity(j) + &
               0.5_r8 * layer_thickness(j+1) / diffusion_conductivity(j+1)
          conductance(j) = 1._r8 / resistance
       end if
    end do

    lower = 0._r8
    diagonal = 1._r8
    upper = 0._r8
    rhs = concentration

    ! Internal interfaces. A positive interface flux moves material from j to
    ! j+1; upwinding and the centered diffusive term form an M-matrix.
    do j = 1, number_of_layers - 1
       left_coefficient = (max(water_flux(j), 0._r8) + conductance(j)) * &
            porewater_factor(j)
       right_coefficient = (min(water_flux(j), 0._r8) - conductance(j)) * &
            porewater_factor(j+1)
       diagonal(j) = diagonal(j) + dt * left_coefficient / layer_thickness(j)
       upper(j) = upper(j) + dt * right_coefficient / layer_thickness(j)
       lower(j+1) = lower(j+1) - dt * left_coefficient / layer_thickness(j+1)
       diagonal(j+1) = diagonal(j+1) - dt * right_coefficient / layer_thickness(j+1)
    end do

    ! Zero-concentration external water: downward top inflow and upward bottom
    ! inflow add no solute; upward top flow and downward bottom flow export it.
    if (water_flux(0) < 0._r8) then
       diagonal(1) = diagonal(1) - dt * water_flux(0) * &
            porewater_factor(1) / layer_thickness(1)
    end if
    if (water_flux(number_of_layers) > 0._r8) then
       diagonal(number_of_layers) = diagonal(number_of_layers) + &
            dt * water_flux(number_of_layers) * porewater_factor(number_of_layers) / &
            layer_thickness(number_of_layers)
    end if

    denominator = max(diagonal(1), tiny(1._r8))
    cprime(1) = upper(1) / denominator
    dprime(1) = rhs(1) / denominator
    do j = 2, number_of_layers
       denominator = max(diagonal(j) - lower(j) * cprime(j-1), tiny(1._r8))
       cprime(j) = upper(j) / denominator
       dprime(j) = (rhs(j) - lower(j) * dprime(j-1)) / denominator
    end do
    updated_concentration(number_of_layers) = dprime(number_of_layers)
    do j = number_of_layers - 1, 1, -1
       updated_concentration(j) = dprime(j) - cprime(j) * updated_concentration(j+1)
    end do

    if (water_flux(0) < 0._r8) then
       advective_flux(0) = water_flux(0) * porewater_factor(1) * &
            updated_concentration(1)
    end if
    do j = 1, number_of_layers - 1
       if (water_flux(j) >= 0._r8) then
          advective_flux(j) = water_flux(j) * porewater_factor(j) * &
               updated_concentration(j)
       else
          advective_flux(j) = water_flux(j) * porewater_factor(j+1) * &
               updated_concentration(j+1)
       end if
       diffusive_flux(j) = conductance(j) * &
            (porewater_factor(j) * updated_concentration(j) - &
             porewater_factor(j+1) * updated_concentration(j+1))
    end do
    if (water_flux(number_of_layers) > 0._r8) then
       advective_flux(number_of_layers) = water_flux(number_of_layers) * &
            porewater_factor(number_of_layers) * &
            updated_concentration(number_of_layers)
    end if

    tendency = (updated_concentration - concentration) / dt
    boundary_export = advective_flux(number_of_layers) + &
         diffusive_flux(number_of_layers) - advective_flux(0) - diffusive_flux(0)
    initial_inventory = sum(concentration * layer_thickness)
    final_inventory = sum(updated_concentration * layer_thickness)
    residual = final_inventory - initial_inventory + dt * boundary_export
    valid = all(updated_concentration >= -state_tolerance) .and. &
         boundary_export >= -state_tolerance .and. &
         aqueousResidualIsClosed(residual, initial_inventory, final_inventory)
  end subroutine advanceMicrobeAqueousTracerTransport

  pure subroutine advanceMicrobeMethaneDOMRelaxation(dom_c, dom_n, dom_p, &
       layer_thickness, layer_relaxation_rate, dt, updated_dom_c, updated_dom_n, &
       updated_dom_p, carbon_residual, nitrogen_residual, phosphorus_residual, valid)
    real(r8), intent(in) :: dom_c(:), dom_n(size(dom_c)), dom_p(size(dom_c))
    real(r8), intent(in) :: layer_thickness(size(dom_c))
    real(r8), intent(in) :: layer_relaxation_rate(size(dom_c))
    real(r8), intent(in) :: dt
    real(r8), intent(out) :: updated_dom_c(size(dom_c))
    real(r8), intent(out) :: updated_dom_n(size(dom_c))
    real(r8), intent(out) :: updated_dom_p(size(dom_c))
    real(r8), intent(out) :: carbon_residual, nitrogen_residual, phosphorus_residual
    logical, intent(out) :: valid
    real(r8) :: initial_carbon, initial_nitrogen, initial_phosphorus
    real(r8) :: final_carbon, final_nitrogen, final_phosphorus

    updated_dom_c = dom_c
    updated_dom_n = dom_n
    updated_dom_p = dom_p
    carbon_residual = 0._r8
    nitrogen_residual = 0._r8
    phosphorus_residual = 0._r8
    valid = .false.
    if (size(dom_c) == 0 .or. dt <= 0._r8) return
    if (any(layer_thickness <= 0._r8) .or. any(layer_relaxation_rate < 0._r8)) return
    if (any(dom_c < -state_tolerance) .or. any(dom_n < -state_tolerance) .or. &
         any(dom_p < -state_tolerance)) return

    ! The CLM-Microbe source relaxes each adjacent DOM pair at dom_diffus
    ! [s-1], even though its input metadata labels that value as m2 s-1. This
    ! backward-Euler finite-volume form preserves that equal-layer timescale
    ! while closing the top and bottom boundaries and conserving inventory.
    ! The identical matrix transports C, N, and P without creating or losing
    ! any element; it also preserves a spatially uniform DOM C:N:P ratio.
    call relaxDOMTracer(dom_c, layer_thickness, layer_relaxation_rate, dt, updated_dom_c)
    call relaxDOMTracer(dom_n, layer_thickness, layer_relaxation_rate, dt, updated_dom_n)
    call relaxDOMTracer(dom_p, layer_thickness, layer_relaxation_rate, dt, updated_dom_p)

    initial_carbon = sum(dom_c * layer_thickness)
    initial_nitrogen = sum(dom_n * layer_thickness)
    initial_phosphorus = sum(dom_p * layer_thickness)
    final_carbon = sum(updated_dom_c * layer_thickness)
    final_nitrogen = sum(updated_dom_n * layer_thickness)
    final_phosphorus = sum(updated_dom_p * layer_thickness)
    carbon_residual = final_carbon - initial_carbon
    nitrogen_residual = final_nitrogen - initial_nitrogen
    phosphorus_residual = final_phosphorus - initial_phosphorus
    valid = all(updated_dom_c >= -state_tolerance) .and. &
         all(updated_dom_n >= -state_tolerance) .and. &
         all(updated_dom_p >= -state_tolerance) .and. &
         residualIsClosed(carbon_residual, initial_carbon, final_carbon) .and. &
         residualIsClosed(nitrogen_residual, initial_nitrogen, final_nitrogen) .and. &
         residualIsClosed(phosphorus_residual, initial_phosphorus, final_phosphorus)
  end subroutine advanceMicrobeMethaneDOMRelaxation

  pure subroutine relaxDOMTracer(concentration, layer_thickness, &
       layer_relaxation_rate, dt, updated_concentration)
    real(r8), intent(in) :: concentration(:)
    real(r8), intent(in) :: layer_thickness(size(concentration))
    real(r8), intent(in) :: layer_relaxation_rate(size(concentration))
    real(r8), intent(in) :: dt
    real(r8), intent(out) :: updated_concentration(size(concentration))
    real(r8) :: conductance(0:size(concentration))
    real(r8) :: lower(size(concentration)), diagonal(size(concentration))
    real(r8) :: upper(size(concentration)), rhs(size(concentration))
    real(r8) :: cprime(size(concentration)), dprime(size(concentration))
    real(r8) :: denominator
    integer :: j, number_of_layers

    number_of_layers = size(concentration)
    updated_concentration = concentration
    if (number_of_layers == 0) return

    conductance = 0._r8
    do j = 1, number_of_layers - 1
       ! Using the deeper-layer rate matches the executed CLM loop. The
       ! harmonic storage depth makes conductance symmetric, hence conservative.
       conductance(j) = max(0._r8, layer_relaxation_rate(j+1)) * &
            2._r8 * layer_thickness(j) * layer_thickness(j+1) / &
            (layer_thickness(j) + layer_thickness(j+1))
    end do

    rhs = concentration
    do j = 1, number_of_layers
       lower(j) = -dt * conductance(j-1) / layer_thickness(j)
       upper(j) = -dt * conductance(j) / layer_thickness(j)
       diagonal(j) = 1._r8 - lower(j) - upper(j)
    end do
    lower(1) = 0._r8
    upper(number_of_layers) = 0._r8

    denominator = max(diagonal(1), tiny(1._r8))
    cprime(1) = upper(1) / denominator
    dprime(1) = rhs(1) / denominator
    do j = 2, number_of_layers
       denominator = max(diagonal(j) - lower(j) * cprime(j-1), tiny(1._r8))
       cprime(j) = upper(j) / denominator
       dprime(j) = (rhs(j) - lower(j) * dprime(j-1)) / denominator
    end do
    updated_concentration(number_of_layers) = dprime(number_of_layers)
    do j = number_of_layers - 1, 1, -1
       updated_concentration(j) = dprime(j) - cprime(j) * updated_concentration(j+1)
    end do
  end subroutine relaxDOMTracer

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

  pure logical function aqueousResidualIsClosed(residual, initial_inventory, &
       final_inventory) result(closed)
    real(r8), intent(in) :: residual, initial_inventory, final_inventory
    real(r8) :: scale

    scale = max(1._r8, abs(initial_inventory), abs(final_inventory))
    closed = abs(residual) <= aqueous_budget_tolerance * scale
  end function aqueousResidualIsClosed

  pure real(r8) function clampUnitInterval(value) result(clamped)
    real(r8), intent(in) :: value
    clamped = min(1._r8, max(0._r8, value))
  end function clampUnitInterval

end module MicrobeMethaneStateUpdateMod
