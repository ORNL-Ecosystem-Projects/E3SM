module MicrobeMethaneIsotopeMod

  ! Carbon-isotope routing for the revised microbial methane backend.
  !
  ! The bulk reaction kernel is the only place that selects reaction rates.
  ! This module applies those accepted rates to a passive carbon tracer.  It
  ! therefore cannot choose a different substrate limitation or reaction
  ! extent from bulk carbon, unlike the historical CLM-SPRUCE implementation.

  use shr_kind_mod, only : r8 => shr_kind_r8
  use elm_varcon, only : catomw
  use MicrobeMethaneParamsMod, only : MicrobeMethaneParamsType
  use MicrobeMethaneReactionMod, only : microbe_methane_reaction_state_type
  use MicrobeMethaneReactionMod, only : microbe_methane_reaction_rates_type
  use MicrobeGasTransportMod, only : computeMicrobeImplicitVerticalDiffusion
  use MicrobeGasTransportMod, only : microbe_gas_ch4, microbe_gas_co2
  use MicrobeGasTransportMod, only : microbe_gas_count

  implicit none
  private
  save

  real(r8), parameter :: tracer_tolerance = 1.e-12_r8

  public :: advanceMicrobeMethaneCarbonTracerReactionLayer
  public :: advanceMicrobeMethaneCarbonTracerGasTransport

contains

  pure subroutine advanceMicrobeMethaneCarbonTracerGasTransport( &
       bulk_state, bulk_state_final, tracer_state, layer_thickness, &
       effective_diffusivity, transport_capacity, atmospheric_mobile_bulk, &
       atmospheric_mobile_tracer, &
       surface_conductance, bulk_interface_flux, bulk_aerenchyma_flux, &
       bulk_ch4_ebullition_loss, dt, tracer_state_final, tracer_interface_flux, &
       tracer_aerenchyma_flux, tracer_ch4_ebullition_loss, &
       tracer_surface_flux, residual, valid)
    ! Diffusion is a linear solve on the tracer itself. Plant exchange and
    ! ebullition are nonlinear bulk decisions, so their already accepted bulk
    ! fluxes are multiplied by the appropriate donor isotope ratio.
    type(microbe_methane_reaction_state_type), intent(in) :: bulk_state(:)
    type(microbe_methane_reaction_state_type), intent(in) :: bulk_state_final(size(bulk_state))
    type(microbe_methane_reaction_state_type), intent(in) :: tracer_state(size(bulk_state))
    real(r8), intent(in) :: layer_thickness(size(bulk_state))
    real(r8), intent(in) :: effective_diffusivity(size(bulk_state),microbe_gas_count)
    real(r8), intent(in) :: transport_capacity(size(bulk_state),microbe_gas_count)
    real(r8), intent(in) :: atmospheric_mobile_bulk(microbe_gas_count)
    real(r8), intent(in) :: atmospheric_mobile_tracer(microbe_gas_count)
    real(r8), intent(in) :: surface_conductance(microbe_gas_count)
    real(r8), intent(in) :: bulk_interface_flux(0:size(bulk_state),microbe_gas_count)
    real(r8), intent(in) :: bulk_aerenchyma_flux(size(bulk_state),microbe_gas_count)
    real(r8), intent(in) :: bulk_ch4_ebullition_loss(size(bulk_state))
    real(r8), intent(in) :: dt
    type(microbe_methane_reaction_state_type), intent(out) :: tracer_state_final(size(bulk_state))
    real(r8), intent(out) :: tracer_interface_flux(0:size(bulk_state),microbe_gas_count)
    real(r8), intent(out) :: tracer_aerenchyma_flux(size(bulk_state),microbe_gas_count)
    real(r8), intent(out) :: tracer_ch4_ebullition_loss(size(bulk_state))
    real(r8), intent(out) :: tracer_surface_flux(microbe_gas_count)
    real(r8), intent(out) :: residual
    logical, intent(out) :: valid
    real(r8) :: bulk_concentration(size(bulk_state),microbe_gas_count)
    real(r8) :: bulk_diffused(size(bulk_state),microbe_gas_count)
    real(r8) :: tracer_concentration(size(bulk_state),microbe_gas_count)
    real(r8) :: tracer_diffused(size(bulk_state),microbe_gas_count)
    real(r8) :: tracer_tendency(size(bulk_state))
    real(r8) :: unused_surface_flux, ratio, transfer
    real(r8) :: initial_inventory, final_inventory, inventory_scale
    integer :: gas, j

    tracer_state_final = tracer_state
    tracer_interface_flux = 0._r8
    tracer_aerenchyma_flux = 0._r8
    tracer_ch4_ebullition_loss = 0._r8
    tracer_surface_flux = 0._r8
    residual = 0._r8
    valid = .false.
    if (size(bulk_state) == 0 .or. dt <= 0._r8 .or. &
         any(layer_thickness <= 0._r8)) return

    bulk_concentration = 0._r8
    tracer_concentration = 0._r8
    do j = 1, size(bulk_state)
       bulk_concentration(j,microbe_gas_ch4) = bulk_state(j)%conc_ch4
       bulk_concentration(j,microbe_gas_co2) = bulk_state(j)%conc_co2
       tracer_concentration(j,microbe_gas_ch4) = tracer_state(j)%conc_ch4
       tracer_concentration(j,microbe_gas_co2) = tracer_state(j)%conc_co2
    end do

    bulk_diffused = bulk_concentration
    tracer_diffused = tracer_concentration
    do gas = microbe_gas_ch4, microbe_gas_co2, microbe_gas_co2-microbe_gas_ch4
       do j = 1, size(bulk_state)
          bulk_diffused(j,gas) = bulk_concentration(j,gas) + dt * &
               (bulk_interface_flux(j-1,gas) - bulk_interface_flux(j,gas)) / &
               layer_thickness(j)
       end do
       call computeMicrobeImplicitVerticalDiffusion(tracer_concentration(:,gas), &
            layer_thickness, effective_diffusivity(:,gas), transport_capacity(:,gas), &
            atmospheric_mobile_tracer(gas), surface_conductance(gas), dt, &
            tracer_diffused(:,gas), tracer_interface_flux(:,gas), tracer_tendency, &
            unused_surface_flux)
    end do

    ! Apply accepted plant fluxes using the post-diffusion soil donor for
    ! outward exchange and the atmospheric isotope ratio for inward exchange.
    do gas = microbe_gas_ch4, microbe_gas_co2, microbe_gas_co2-microbe_gas_ch4
       do j = 1, size(bulk_state)
          if (bulk_aerenchyma_flux(j,gas) >= 0._r8) then
             ratio = tracerFraction(tracer_diffused(j,gas), bulk_diffused(j,gas))
          else
             ratio = tracerFraction(atmospheric_mobile_tracer(gas), &
                  atmospheric_mobile_bulk(gas))
          end if
          tracer_aerenchyma_flux(j,gas) = bulk_aerenchyma_flux(j,gas) * ratio
       end do
    end do
    do j = 1, size(bulk_state)
       ratio = tracerFraction(tracer_diffused(j,microbe_gas_ch4), &
            bulk_diffused(j,microbe_gas_ch4))
       tracer_ch4_ebullition_loss(j) = bulk_ch4_ebullition_loss(j) * ratio
       tracer_diffused(j,microbe_gas_ch4) = &
            tracer_diffused(j,microbe_gas_ch4) - dt * &
            (tracer_aerenchyma_flux(j,microbe_gas_ch4) + &
            tracer_ch4_ebullition_loss(j))
       tracer_diffused(j,microbe_gas_co2) = &
            tracer_diffused(j,microbe_gas_co2) - dt * &
            tracer_aerenchyma_flux(j,microbe_gas_co2)
       tracer_state_final(j)%conc_ch4 = tracer_diffused(j,microbe_gas_ch4)
       tracer_state_final(j)%conc_co2 = tracer_diffused(j,microbe_gas_co2)
       tracer_state_final(j)%conc_o2 = 0._r8
       tracer_state_final(j)%conc_h2 = 0._r8
    end do

    do gas = microbe_gas_ch4, microbe_gas_co2, microbe_gas_co2-microbe_gas_ch4
       tracer_surface_flux(gas) = -tracer_interface_flux(0,gas) + &
            sum(tracer_aerenchyma_flux(:,gas) * layer_thickness)
    end do
    tracer_surface_flux(microbe_gas_ch4) = &
         tracer_surface_flux(microbe_gas_ch4) + &
         sum(tracer_ch4_ebullition_loss * layer_thickness)

    initial_inventory = catomw * sum((tracer_concentration(:,microbe_gas_ch4) + &
         tracer_concentration(:,microbe_gas_co2)) * layer_thickness)
    final_inventory = catomw * sum((tracer_diffused(:,microbe_gas_ch4) + &
         tracer_diffused(:,microbe_gas_co2)) * layer_thickness)
    residual = final_inventory - initial_inventory + dt * catomw * &
         (tracer_surface_flux(microbe_gas_ch4) + tracer_surface_flux(microbe_gas_co2))
    inventory_scale = max(1._r8, abs(initial_inventory), abs(final_inventory))
    valid = all(tracer_diffused(:,microbe_gas_ch4) >= -tracer_tolerance) .and. &
         all(tracer_diffused(:,microbe_gas_co2) >= -tracer_tolerance) .and. &
         all(tracer_diffused(:,microbe_gas_ch4) <= &
         [(max(0._r8, bulk_state_final(j)%conc_ch4) + tracer_tolerance, &
         j=1,size(bulk_state))]) .and. &
         all(tracer_diffused(:,microbe_gas_co2) <= &
         [(max(0._r8, bulk_state_final(j)%conc_co2) + tracer_tolerance, &
         j=1,size(bulk_state))]) .and. &
         abs(residual) <= tracer_tolerance * inventory_scale
  end subroutine advanceMicrobeMethaneCarbonTracerGasTransport

  pure subroutine advanceMicrobeMethaneCarbonTracerReactionLayer( &
       bulk_dom_c, bulk_dom_c_final, saturated_fraction, &
       bulk_unsaturated_state, bulk_saturated_state, &
       bulk_unsaturated_state_final, bulk_saturated_state_final, &
       unsaturated_rates, saturated_rates, parameters, dt, &
       tracer_dom_c, tracer_unsaturated_state, tracer_saturated_state, &
       tracer_dom_c_final, tracer_unsaturated_state_final, &
       tracer_saturated_state_final, residual, valid)
    real(r8), intent(in) :: bulk_dom_c, bulk_dom_c_final, saturated_fraction, dt
    type(microbe_methane_reaction_state_type), intent(in) :: bulk_unsaturated_state
    type(microbe_methane_reaction_state_type), intent(in) :: bulk_saturated_state
    type(microbe_methane_reaction_state_type), intent(in) :: bulk_unsaturated_state_final
    type(microbe_methane_reaction_state_type), intent(in) :: bulk_saturated_state_final
    type(microbe_methane_reaction_rates_type), intent(in) :: unsaturated_rates
    type(microbe_methane_reaction_rates_type), intent(in) :: saturated_rates
    type(MicrobeMethaneParamsType), intent(in) :: parameters
    real(r8), intent(in) :: tracer_dom_c
    type(microbe_methane_reaction_state_type), intent(in) :: tracer_unsaturated_state
    type(microbe_methane_reaction_state_type), intent(in) :: tracer_saturated_state
    real(r8), intent(out) :: tracer_dom_c_final
    type(microbe_methane_reaction_state_type), intent(out) :: tracer_unsaturated_state_final
    type(microbe_methane_reaction_state_type), intent(out) :: tracer_saturated_state_final
    real(r8), intent(out) :: residual
    logical, intent(out) :: valid
    real(r8) :: fraction, unsaturated_dom_delta, saturated_dom_delta
    real(r8) :: initial_inventory, final_inventory, inventory_scale
    logical :: unsaturated_valid, saturated_valid

    fraction = clampUnitInterval(saturated_fraction)
    tracer_dom_c_final = tracer_dom_c
    tracer_unsaturated_state_final = tracer_unsaturated_state
    tracer_saturated_state_final = tracer_saturated_state
    residual = 0._r8
    valid = .false.
    if (dt <= 0._r8 .or. bulk_dom_c < -tracer_tolerance .or. &
         tracer_dom_c < -tracer_tolerance) return

    call advancePartition(bulk_dom_c, bulk_unsaturated_state, &
         bulk_unsaturated_state_final, unsaturated_rates, parameters, dt, &
         tracer_dom_c, tracer_unsaturated_state, unsaturated_dom_delta, &
         tracer_unsaturated_state_final, unsaturated_valid)
    call advancePartition(bulk_dom_c, bulk_saturated_state, &
         bulk_saturated_state_final, saturated_rates, parameters, dt, &
         tracer_dom_c, tracer_saturated_state, saturated_dom_delta, &
         tracer_saturated_state_final, saturated_valid)

    tracer_dom_c_final = tracer_dom_c + &
         (1._r8 - fraction) * unsaturated_dom_delta + &
         fraction * saturated_dom_delta

    initial_inventory = tracer_dom_c + &
         (1._r8 - fraction) * tracerCarbonDensity(tracer_unsaturated_state) + &
         fraction * tracerCarbonDensity(tracer_saturated_state)
    final_inventory = tracer_dom_c_final + &
         (1._r8 - fraction) * tracerCarbonDensity(tracer_unsaturated_state_final) + &
         fraction * tracerCarbonDensity(tracer_saturated_state_final)
    residual = final_inventory - initial_inventory
    inventory_scale = max(1._r8, abs(initial_inventory), abs(final_inventory))

    valid = unsaturated_valid .and. saturated_valid .and. &
         tracer_dom_c_final >= -tracer_tolerance .and. &
         tracer_dom_c_final <= max(0._r8, bulk_dom_c_final) + tracer_tolerance .and. &
         abs(residual) <= tracer_tolerance * inventory_scale
  end subroutine advanceMicrobeMethaneCarbonTracerReactionLayer

  pure subroutine advancePartition(bulk_dom_c, bulk_state, bulk_state_final, rates, &
       parameters, dt, tracer_dom_c, tracer_state, tracer_dom_delta, &
       tracer_state_final, valid)
    real(r8), intent(in) :: bulk_dom_c, dt, tracer_dom_c
    type(microbe_methane_reaction_state_type), intent(in) :: bulk_state
    type(microbe_methane_reaction_state_type), intent(in) :: bulk_state_final
    type(microbe_methane_reaction_rates_type), intent(in) :: rates
    type(MicrobeMethaneParamsType), intent(in) :: parameters
    type(microbe_methane_reaction_state_type), intent(in) :: tracer_state
    real(r8), intent(out) :: tracer_dom_delta
    type(microbe_methane_reaction_state_type), intent(out) :: tracer_state_final
    logical, intent(out) :: valid
    type(microbe_methane_reaction_state_type) :: bulk_work, tracer_work
    real(r8) :: extent_dom, extent_acetogenesis, extent_h2_methanogenesis
    real(r8) :: extent_acetoclastic, extent_aerobic_acetate
    real(r8) :: extent_aerobic_ch4, extent_anaerobic_ch4
    real(r8) :: donor, ratio, tracer_transfer, yield, ch4_yield
    real(r8) :: initial_biomass_ratio(4), mortality_extent(4)
    integer :: guild

    bulk_work = bulk_state
    tracer_work = tracer_state
    ! O2 and H2 contain no carbon isotope. Keep them identically zero so this
    ! carbon-only state cannot accidentally be interpreted as a full gas state.
    tracer_work%dom_c = 0._r8
    tracer_work%conc_o2 = 0._r8
    tracer_work%conc_h2 = 0._r8
    tracer_dom_delta = 0._r8

    initial_biomass_ratio(1) = tracerFraction( &
         tracer_state%acetate_methanogen_c, bulk_state%acetate_methanogen_c)
    initial_biomass_ratio(2) = tracerFraction( &
         tracer_state%h2_methanogen_c, bulk_state%h2_methanogen_c)
    initial_biomass_ratio(3) = tracerFraction( &
         tracer_state%aerobic_methanotroph_c, bulk_state%aerobic_methanotroph_c)
    initial_biomass_ratio(4) = tracerFraction( &
         tracer_state%anaerobic_methanotroph_c, bulk_state%anaerobic_methanotroph_c)
    mortality_extent = dt * [rates%acetate_methanogen_mortality_c, &
         rates%h2_methanogen_mortality_c, &
         rates%aerobic_methanotroph_mortality_c, &
         rates%anaerobic_methanotroph_mortality_c]

    ! Stage 1: 1.5 DOM-C -> 1 acetate-C + 0.5 CO2-C.
    extent_dom = dt * rates%dom_to_acetate_c
    donor = 1.5_r8 * catomw * extent_dom
    ratio = tracerFraction(tracer_dom_c, bulk_dom_c)
    tracer_transfer = donor * ratio
    tracer_dom_delta = tracer_dom_delta - tracer_transfer
    bulk_work%acetate_c = bulk_work%acetate_c + catomw * extent_dom
    tracer_work%acetate_c = tracer_work%acetate_c + tracer_transfer / 1.5_r8
    bulk_work%conc_co2 = bulk_work%conc_co2 + 0.5_r8 * extent_dom
    tracer_work%conc_co2 = tracer_work%conc_co2 + tracer_transfer / (3._r8 * catomw)

    ! Stage 2: acetogenesis and hydrogenotrophic methanogenesis withdraw from
    ! one mixed CO2 inventory. Use one donor ratio for both accepted extents.
    extent_acetogenesis = dt * rates%acetogenesis_c
    extent_h2_methanogenesis = dt * rates%hydrogenotrophic_methanogenesis_c
    yield = parameters%h2_methanogen_yield
    donor = extent_acetogenesis + (1._r8 + yield) * extent_h2_methanogenesis
    ratio = tracerFraction(tracer_work%conc_co2, bulk_work%conc_co2)
    tracer_work%conc_co2 = tracer_work%conc_co2 - donor * ratio
    bulk_work%conc_co2 = bulk_work%conc_co2 - donor
    bulk_work%acetate_c = bulk_work%acetate_c + catomw * extent_acetogenesis
    tracer_work%acetate_c = tracer_work%acetate_c + &
         catomw * extent_acetogenesis * ratio
    bulk_work%conc_ch4 = bulk_work%conc_ch4 + extent_h2_methanogenesis
    tracer_work%conc_ch4 = tracer_work%conc_ch4 + extent_h2_methanogenesis * ratio
    bulk_work%h2_methanogen_c = bulk_work%h2_methanogen_c + &
         catomw * yield * extent_h2_methanogenesis
    tracer_work%h2_methanogen_c = tracer_work%h2_methanogen_c + &
         catomw * yield * extent_h2_methanogenesis * ratio

    ! Stage 3: acetoclastic methanogenesis and aerobic acetate oxidation draw
    ! from the same mixed acetate inventory.
    extent_acetoclastic = dt * rates%acetoclastic_methanogenesis_c
    extent_aerobic_acetate = dt * rates%aerobic_acetate_oxidation_c
    donor = catomw * (extent_acetoclastic + extent_aerobic_acetate)
    ratio = tracerFraction(tracer_work%acetate_c, bulk_work%acetate_c)
    tracer_work%acetate_c = tracer_work%acetate_c - donor * ratio
    bulk_work%acetate_c = bulk_work%acetate_c - donor
    yield = parameters%acetate_methanogen_yield
    ch4_yield = parameters%acetoclastic_methanogenesis_ch4_yield
    bulk_work%acetate_methanogen_c = bulk_work%acetate_methanogen_c + &
         catomw * yield * extent_acetoclastic
    tracer_work%acetate_methanogen_c = tracer_work%acetate_methanogen_c + &
         catomw * yield * extent_acetoclastic * ratio
    bulk_work%conc_ch4 = bulk_work%conc_ch4 + &
         ch4_yield * (1._r8 - yield) * extent_acetoclastic
    tracer_work%conc_ch4 = tracer_work%conc_ch4 + &
         ch4_yield * (1._r8 - yield) * extent_acetoclastic * ratio
    bulk_work%conc_co2 = bulk_work%conc_co2 + &
         (1._r8 - ch4_yield) * (1._r8 - yield) * extent_acetoclastic + &
         extent_aerobic_acetate
    tracer_work%conc_co2 = tracer_work%conc_co2 + &
         ((1._r8 - ch4_yield) * (1._r8 - yield) * extent_acetoclastic + &
         extent_aerobic_acetate) * ratio

    ! Stage 4: both methane oxidation pathways withdraw from one mixed CH4
    ! inventory and partition that carbon between biomass and CO2.
    extent_aerobic_ch4 = dt * rates%aerobic_methane_oxidation_c
    extent_anaerobic_ch4 = dt * rates%anaerobic_methane_oxidation_c
    donor = extent_aerobic_ch4 + extent_anaerobic_ch4
    ratio = tracerFraction(tracer_work%conc_ch4, bulk_work%conc_ch4)
    tracer_work%conc_ch4 = tracer_work%conc_ch4 - donor * ratio
    bulk_work%conc_ch4 = bulk_work%conc_ch4 - donor
    yield = parameters%aerobic_methanotroph_yield
    bulk_work%aerobic_methanotroph_c = bulk_work%aerobic_methanotroph_c + &
         catomw * yield * extent_aerobic_ch4
    tracer_work%aerobic_methanotroph_c = tracer_work%aerobic_methanotroph_c + &
         catomw * yield * extent_aerobic_ch4 * ratio
    bulk_work%conc_co2 = bulk_work%conc_co2 + (1._r8 - yield) * extent_aerobic_ch4
    tracer_work%conc_co2 = tracer_work%conc_co2 + &
         (1._r8 - yield) * extent_aerobic_ch4 * ratio
    yield = parameters%anaerobic_methanotroph_yield
    bulk_work%anaerobic_methanotroph_c = bulk_work%anaerobic_methanotroph_c + &
         catomw * yield * extent_anaerobic_ch4
    tracer_work%anaerobic_methanotroph_c = tracer_work%anaerobic_methanotroph_c + &
         catomw * yield * extent_anaerobic_ch4 * ratio
    bulk_work%conc_co2 = bulk_work%conc_co2 + (1._r8 - yield) * extent_anaerobic_ch4
    tracer_work%conc_co2 = tracer_work%conc_co2 + &
         (1._r8 - yield) * extent_anaerobic_ch4 * ratio

    ! Mortality applies only to biomass present at the start of the call; the
    ! bulk limiter explicitly excludes same-step growth. Route it with each
    ! guild's initial tracer ratio and return it to the shared DOM pool.
    do guild = 1, 4
       tracer_transfer = catomw * mortality_extent(guild) * initial_biomass_ratio(guild)
       select case (guild)
       case (1)
          tracer_work%acetate_methanogen_c = &
               tracer_work%acetate_methanogen_c - tracer_transfer
       case (2)
          tracer_work%h2_methanogen_c = &
               tracer_work%h2_methanogen_c - tracer_transfer
       case (3)
          tracer_work%aerobic_methanotroph_c = &
               tracer_work%aerobic_methanotroph_c - tracer_transfer
       case (4)
          tracer_work%anaerobic_methanotroph_c = &
               tracer_work%anaerobic_methanotroph_c - tracer_transfer
       end select
       tracer_dom_delta = tracer_dom_delta + tracer_transfer
    end do

    tracer_state_final = tracer_work
    valid = tracerStateIsValid(tracer_state_final, bulk_state_final)
  end subroutine advancePartition

  pure logical function tracerStateIsValid(tracer_state, bulk_state) result(valid)
    type(microbe_methane_reaction_state_type), intent(in) :: tracer_state, bulk_state
    valid = tracer_state%acetate_c >= -tracer_tolerance .and. &
         tracer_state%acetate_methanogen_c >= -tracer_tolerance .and. &
         tracer_state%h2_methanogen_c >= -tracer_tolerance .and. &
         tracer_state%aerobic_methanotroph_c >= -tracer_tolerance .and. &
         tracer_state%anaerobic_methanotroph_c >= -tracer_tolerance .and. &
         tracer_state%conc_ch4 >= -tracer_tolerance .and. &
         tracer_state%conc_co2 >= -tracer_tolerance .and. &
         tracer_state%acetate_c <= max(0._r8, bulk_state%acetate_c) + tracer_tolerance .and. &
         tracer_state%acetate_methanogen_c <= &
         max(0._r8, bulk_state%acetate_methanogen_c) + tracer_tolerance .and. &
         tracer_state%h2_methanogen_c <= &
         max(0._r8, bulk_state%h2_methanogen_c) + tracer_tolerance .and. &
         tracer_state%aerobic_methanotroph_c <= &
         max(0._r8, bulk_state%aerobic_methanotroph_c) + tracer_tolerance .and. &
         tracer_state%anaerobic_methanotroph_c <= &
         max(0._r8, bulk_state%anaerobic_methanotroph_c) + tracer_tolerance .and. &
         tracer_state%conc_ch4 <= max(0._r8, bulk_state%conc_ch4) + tracer_tolerance .and. &
         tracer_state%conc_co2 <= max(0._r8, bulk_state%conc_co2) + tracer_tolerance
  end function tracerStateIsValid

  pure real(r8) function tracerCarbonDensity(state) result(carbon)
    type(microbe_methane_reaction_state_type), intent(in) :: state
    carbon = state%acetate_c + state%acetate_methanogen_c + &
         state%h2_methanogen_c + state%aerobic_methanotroph_c + &
         state%anaerobic_methanotroph_c + &
         catomw * (state%conc_ch4 + state%conc_co2)
  end function tracerCarbonDensity

  pure real(r8) function tracerFraction(tracer_carbon, bulk_carbon) result(fraction)
    real(r8), intent(in) :: tracer_carbon, bulk_carbon
    if (bulk_carbon <= tiny(1._r8)) then
       fraction = 0._r8
    else
       fraction = clampUnitInterval(max(0._r8, tracer_carbon) / bulk_carbon)
    end if
  end function tracerFraction

  pure real(r8) function clampUnitInterval(value) result(clamped)
    real(r8), intent(in) :: value
    clamped = min(1._r8, max(0._r8, value))
  end function clampUnitInterval

end module MicrobeMethaneIsotopeMod
