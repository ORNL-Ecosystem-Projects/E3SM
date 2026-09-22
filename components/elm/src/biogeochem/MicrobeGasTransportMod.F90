module MicrobeGasTransportMod

  ! Conservative, side-effect-free transport kernels for the revised
  ! microbial methane backend. Gas state is bulk-soil inventory density
  ! (mol m-3 soil); acetate is g C m-3 soil. Flux and tendency units follow
  ! the input unit.

  use shr_kind_mod, only : r8 => shr_kind_r8
  use MicrobeMethaneParamsMod, only : MicrobeMethaneParamsType

  implicit none
  private
  save

  integer, parameter, public :: microbe_gas_ch4 = 1
  integer, parameter, public :: microbe_gas_o2 = 2
  integer, parameter, public :: microbe_gas_co2 = 3
  integer, parameter, public :: microbe_gas_h2 = 4
  integer, parameter, public :: microbe_gas_count = 4

  public :: repartitionMicrobeMethaneScalar
  public :: microbeMethaneBulkConcentration
  public :: computeMicrobeGasTransport
  public :: computeMicrobeVerticalDiffusion
  public :: computeMicrobeImplicitVerticalDiffusion
  public :: computeMicrobeAerenchymaTransport
  public :: computeMicrobeMethaneEbullition
  public :: computeMicrobeTopounitLateralDiffusion
  public :: microbeMethaneEffectiveAqueousDiffusivity
  public :: microbeMethaneTransportResidual

contains

  pure subroutine repartitionMicrobeMethaneScalar(old_saturated_fraction, new_saturated_fraction, &
       unsaturated_concentration, saturated_concentration, repartitioned_unsaturated, &
       repartitioned_saturated)
    real(r8), intent(in) :: old_saturated_fraction, new_saturated_fraction
    real(r8), intent(in) :: unsaturated_concentration, saturated_concentration
    real(r8), intent(out) :: repartitioned_unsaturated, repartitioned_saturated
    real(r8) :: old_fraction, new_fraction, transferred_fraction

    old_fraction = clampUnitInterval(old_saturated_fraction)
    new_fraction = clampUnitInterval(new_saturated_fraction)
    repartitioned_unsaturated = unsaturated_concentration
    repartitioned_saturated = saturated_concentration

    if (new_fraction > old_fraction) then
       ! Newly saturated area carries the unsaturated concentration into the
       ! enlarged saturated partition.
       transferred_fraction = new_fraction - old_fraction
       repartitioned_saturated = (old_fraction * saturated_concentration + &
            transferred_fraction * unsaturated_concentration) / new_fraction
    else if (new_fraction < old_fraction) then
       ! Newly unsaturated area carries the saturated concentration into the
       ! enlarged unsaturated partition.
       transferred_fraction = old_fraction - new_fraction
       repartitioned_unsaturated = ((1._r8 - old_fraction) * unsaturated_concentration + &
            transferred_fraction * saturated_concentration) / (1._r8 - new_fraction)
    end if
  end subroutine repartitionMicrobeMethaneScalar

  pure real(r8) function microbeMethaneBulkConcentration(saturated_fraction, &
       unsaturated_concentration, saturated_concentration) result(bulk_concentration)
    real(r8), intent(in) :: saturated_fraction
    real(r8), intent(in) :: unsaturated_concentration, saturated_concentration
    real(r8) :: fraction

    fraction = clampUnitInterval(saturated_fraction)
    bulk_concentration = (1._r8 - fraction) * unsaturated_concentration + &
         fraction * saturated_concentration
  end function microbeMethaneBulkConcentration

  pure subroutine computeMicrobeVerticalDiffusion(concentration, layer_thickness, &
       effective_diffusivity, atmospheric_equivalent_concentration, surface_conductance, &
       dt, interface_flux, tendency, surface_flux, apply_donor_limit)
    real(r8), intent(in) :: concentration(:)
    real(r8), intent(in) :: layer_thickness(size(concentration))
    real(r8), intent(in) :: effective_diffusivity(size(concentration))
    real(r8), intent(in) :: atmospheric_equivalent_concentration
    real(r8), intent(in) :: surface_conductance
    real(r8), intent(in) :: dt
    ! Interface flux is positive downward. Index 0 is the atmosphere/soil
    ! boundary and index n is the closed lower boundary.
    real(r8), intent(out) :: interface_flux(0:size(concentration))
    real(r8), intent(out) :: tendency(size(concentration))
    ! Surface flux is positive from soil to atmosphere.
    real(r8), intent(out) :: surface_flux
    ! This explicit kernel remains in use for closed-boundary acetate
    ! transport. Gas diffusion uses the implicit capacity-aware kernel below.
    logical, intent(in), optional :: apply_donor_limit
    integer :: j, number_of_layers
    real(r8) :: resistance
    real(r8) :: donor_scale(size(concentration))
    real(r8) :: outgoing_flux
    logical :: limit_inventory

    number_of_layers = size(concentration)
    interface_flux = 0._r8
    tendency = 0._r8
    surface_flux = 0._r8
    if (number_of_layers == 0 .or. dt <= 0._r8) return
    limit_inventory = .true.
    if (present(apply_donor_limit)) limit_inventory = apply_donor_limit

    ! The supplied conductance includes the top half-layer, snow, ponded-water,
    ! and atmospheric-boundary resistances selected by the ELM adapter.
    interface_flux(0) = max(0._r8, surface_conductance) * &
         (atmospheric_equivalent_concentration - concentration(1))

    do j = 1, number_of_layers - 1
       if (effective_diffusivity(j) > 0._r8 .and. &
            effective_diffusivity(j+1) > 0._r8 .and. &
            layer_thickness(j) > 0._r8 .and. layer_thickness(j+1) > 0._r8) then
          resistance = 0.5_r8 * layer_thickness(j) / effective_diffusivity(j) + &
               0.5_r8 * layer_thickness(j+1) / effective_diffusivity(j+1)
          interface_flux(j) = (concentration(j) - concentration(j+1)) / resistance
       end if
    end do
    interface_flux(number_of_layers) = 0._r8

    ! Scale all simultaneous outflows from a donor layer by the same factor.
    ! This retains equal-and-opposite internal fluxes and prevents explicit
    ! transport from making a nonnegative layer negative over dt.
    if (limit_inventory) then
       do j = 1, number_of_layers
          outgoing_flux = max(0._r8, -interface_flux(j-1)) + &
               max(0._r8, interface_flux(j))
          if (outgoing_flux > 0._r8 .and. layer_thickness(j) > 0._r8) then
             donor_scale(j) = clampUnitInterval(max(0._r8, concentration(j)) * &
                  layer_thickness(j) / (dt * outgoing_flux))
          else
             donor_scale(j) = 1._r8
          end if
       end do

       if (interface_flux(0) < 0._r8) interface_flux(0) = interface_flux(0) * donor_scale(1)
       do j = 1, number_of_layers - 1
          if (interface_flux(j) > 0._r8) then
             interface_flux(j) = interface_flux(j) * donor_scale(j)
          else if (interface_flux(j) < 0._r8) then
             interface_flux(j) = interface_flux(j) * donor_scale(j+1)
          end if
       end do
    end if

    do j = 1, number_of_layers
       if (layer_thickness(j) > 0._r8) then
          tendency(j) = (interface_flux(j-1) - interface_flux(j)) / layer_thickness(j)
       end if
    end do
    surface_flux = -interface_flux(0)
  end subroutine computeMicrobeVerticalDiffusion

  pure subroutine computeMicrobeImplicitVerticalDiffusion(concentration, layer_thickness, &
       effective_diffusivity, transport_capacity, atmospheric_mobile_concentration, &
       surface_conductance, dt, updated_concentration, interface_flux, tendency, &
       surface_flux)
    ! Backward-Euler finite-volume diffusion for a bulk-soil inventory.
    ! transport_capacity maps the mobile gas-equivalent concentration to the
    ! stored inventory: C_bulk = capacity * C_mobile. effective_diffusivity is
    ! expressed on the same mobile-concentration basis. This permits gaseous
    ! and aqueous layers to share one conservative solve through Henry-law
    ! capacity and diffusivity conversions supplied by the ELM adapter.
    real(r8), intent(in) :: concentration(:)
    real(r8), intent(in) :: layer_thickness(size(concentration))
    real(r8), intent(in) :: effective_diffusivity(size(concentration))
    real(r8), intent(in) :: transport_capacity(size(concentration))
    real(r8), intent(in) :: atmospheric_mobile_concentration
    real(r8), intent(in) :: surface_conductance
    real(r8), intent(in) :: dt
    real(r8), intent(out) :: updated_concentration(size(concentration))
    real(r8), intent(out) :: interface_flux(0:size(concentration))
    real(r8), intent(out) :: tendency(size(concentration))
    real(r8), intent(out) :: surface_flux
    real(r8) :: capacity(size(concentration))
    real(r8) :: conductance(0:size(concentration))
    real(r8) :: lower(size(concentration)), diagonal(size(concentration))
    real(r8) :: upper(size(concentration)), rhs(size(concentration))
    real(r8) :: cprime(size(concentration)), dprime(size(concentration))
    real(r8) :: mobile_concentration(size(concentration))
    real(r8) :: denominator, resistance
    integer :: j, number_of_layers

    number_of_layers = size(concentration)
    updated_concentration = concentration
    interface_flux = 0._r8
    tendency = 0._r8
    surface_flux = 0._r8
    if (number_of_layers == 0 .or. dt <= 0._r8) return

    capacity = max(transport_capacity, tiny(1._r8))
    conductance = 0._r8
    conductance(0) = max(0._r8, surface_conductance)
    do j = 1, number_of_layers - 1
       if (effective_diffusivity(j) > 0._r8 .and. &
            effective_diffusivity(j+1) > 0._r8 .and. &
            layer_thickness(j) > 0._r8 .and. layer_thickness(j+1) > 0._r8) then
          resistance = 0.5_r8 * layer_thickness(j) / effective_diffusivity(j) + &
               0.5_r8 * layer_thickness(j+1) / effective_diffusivity(j+1)
          conductance(j) = 1._r8 / resistance
       end if
    end do
    conductance(number_of_layers) = 0._r8

    lower = 0._r8
    diagonal = capacity
    upper = 0._r8
    rhs = concentration
    do j = 1, number_of_layers
       if (layer_thickness(j) > 0._r8) then
          lower(j) = -dt * conductance(j-1) / layer_thickness(j)
          upper(j) = -dt * conductance(j) / layer_thickness(j)
          diagonal(j) = capacity(j) - lower(j) - upper(j)
       end if
    end do
    rhs(1) = rhs(1) - lower(1) * max(0._r8, atmospheric_mobile_concentration)
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
    mobile_concentration(number_of_layers) = dprime(number_of_layers)
    do j = number_of_layers - 1, 1, -1
       mobile_concentration(j) = dprime(j) - cprime(j) * mobile_concentration(j+1)
    end do

    updated_concentration = capacity * mobile_concentration
    interface_flux(0) = conductance(0) * &
         (max(0._r8, atmospheric_mobile_concentration) - mobile_concentration(1))
    do j = 1, number_of_layers - 1
       interface_flux(j) = conductance(j) * &
            (mobile_concentration(j) - mobile_concentration(j+1))
    end do
    interface_flux(number_of_layers) = 0._r8
    do j = 1, number_of_layers
       if (layer_thickness(j) > 0._r8) then
          tendency(j) = (interface_flux(j-1) - interface_flux(j)) / &
               layer_thickness(j)
          ! Use the conservative flux divergence as the authoritative update;
          ! this keeps the budget closed to roundoff after the tridiagonal solve.
          updated_concentration(j) = concentration(j) + dt * tendency(j)
       end if
    end do
    surface_flux = -interface_flux(0)
  end subroutine computeMicrobeImplicitVerticalDiffusion

  pure subroutine computeMicrobeAerenchymaTransport(concentration, &
       atmospheric_equivalent_concentration, layer_thickness, layer_exchange_rate, &
       dt, flux_to_atmosphere, tendency, surface_flux, apply_donor_limit, &
       minimum_emission_concentration, allow_influx)
    real(r8), intent(in) :: concentration(:)
    real(r8), intent(in) :: atmospheric_equivalent_concentration(size(concentration))
    real(r8), intent(in) :: layer_thickness(size(concentration))
    ! First-order exchange rate supplied by the ELM root/aerenchyma adapter (s-1).
    real(r8), intent(in) :: layer_exchange_rate(size(concentration))
    real(r8), intent(in) :: dt
    ! Layer flux is positive from soil to atmosphere (input-unit m-3 s-1).
    real(r8), intent(out) :: flux_to_atmosphere(size(concentration))
    real(r8), intent(out) :: tendency(size(concentration))
    real(r8), intent(out) :: surface_flux
    logical, intent(in), optional :: apply_donor_limit
    ! A one-way plant pathway emits only above this concentration. This is
    ! distinct from the atmospheric Henry-law equilibrium used by diffusion.
    real(r8), intent(in), optional :: minimum_emission_concentration(size(concentration))
    logical, intent(in), optional :: allow_influx
    integer :: j
    real(r8) :: equilibrium_flux, exchange_target
    logical :: limit_inventory, permit_influx

    flux_to_atmosphere = 0._r8
    tendency = 0._r8
    surface_flux = 0._r8
    if (dt <= 0._r8) return
    limit_inventory = .true.
    if (present(apply_donor_limit)) limit_inventory = apply_donor_limit
    permit_influx = .true.
    if (present(allow_influx)) permit_influx = allow_influx

    do j = 1, size(concentration)
       exchange_target = atmospheric_equivalent_concentration(j)
       if (.not. permit_influx .and. present(minimum_emission_concentration)) then
          exchange_target = max(exchange_target, &
               max(0._r8, minimum_emission_concentration(j)))
       end if
       flux_to_atmosphere(j) = max(0._r8, layer_exchange_rate(j)) * &
            (concentration(j) - exchange_target)
       if (.not. permit_influx) flux_to_atmosphere(j) = &
            max(0._r8, flux_to_atmosphere(j))
       ! Explicit exchange must not cross its reservoir equilibrium in one
       ! timestep. This bound applies in both directions independently of the
       ! joint soil-inventory limiter used by the combined transport driver.
       equilibrium_flux = (concentration(j) - exchange_target) / dt
       if (flux_to_atmosphere(j) > 0._r8) then
          flux_to_atmosphere(j) = min(flux_to_atmosphere(j), &
               max(0._r8, equilibrium_flux))
       else if (flux_to_atmosphere(j) < 0._r8) then
          flux_to_atmosphere(j) = max(flux_to_atmosphere(j), &
               min(0._r8, equilibrium_flux))
       end if
       if (limit_inventory .and. flux_to_atmosphere(j) > 0._r8) then
          flux_to_atmosphere(j) = min(flux_to_atmosphere(j), &
               max(0._r8, concentration(j)) / dt)
       end if
       tendency(j) = -flux_to_atmosphere(j)
       surface_flux = surface_flux + flux_to_atmosphere(j) * &
            max(0._r8, layer_thickness(j))
    end do
  end subroutine computeMicrobeAerenchymaTransport

  pure subroutine computeMicrobeMethaneEbullition(concentration, threshold, &
       activation_fraction, layer_thickness, dt, loss_rate, tendency, surface_flux)
    real(r8), intent(in) :: concentration(:)
    real(r8), intent(in) :: threshold(size(concentration))
    ! Hydrology, thaw state, and depth attenuation are combined by the ELM
    ! adapter into a bounded activation fraction.
    real(r8), intent(in) :: activation_fraction(size(concentration))
    real(r8), intent(in) :: layer_thickness(size(concentration))
    real(r8), intent(in) :: dt
    real(r8), intent(out) :: loss_rate(size(concentration))
    real(r8), intent(out) :: tendency(size(concentration))
    real(r8), intent(out) :: surface_flux
    integer :: j

    loss_rate = 0._r8
    tendency = 0._r8
    surface_flux = 0._r8
    if (dt <= 0._r8) return

    do j = 1, size(concentration)
       loss_rate(j) = clampUnitInterval(activation_fraction(j)) * &
            max(0._r8, concentration(j) - max(0._r8, threshold(j))) / dt
       tendency(j) = -loss_rate(j)
       surface_flux = surface_flux + loss_rate(j) * max(0._r8, layer_thickness(j))
    end do
  end subroutine computeMicrobeMethaneEbullition

  pure subroutine computeMicrobeGasTransport(concentration, layer_thickness, &
       effective_diffusivity, transport_capacity, surface_equilibrium_concentration, &
       surface_conductance, &
       aerenchyma_equilibrium_concentration, aerenchyma_exchange_rate, &
       aerenchyma_minimum_emission_concentration, aerenchyma_allows_influx, &
       ch4_ebullition_threshold, ch4_ebullition_activation, dt, interface_flux, &
       aerenchyma_flux, ch4_ebullition_loss, tendency, surface_diffusive_flux, &
       surface_aerenchyma_flux, surface_ebullition_flux, total_surface_flux)
    ! Gas dimension order is CH4, O2, CO2, H2. All concentrations use
    ! mol gas m-3 soil; layer tendencies use mol gas m-3 s-1 and surface
    ! fluxes use mol gas m-2 s-1, positive to the atmosphere.
    real(r8), intent(in) :: concentration(:,:)
    real(r8), intent(in) :: layer_thickness(size(concentration,1))
    real(r8), intent(in) :: effective_diffusivity(size(concentration,1),microbe_gas_count)
    real(r8), intent(in) :: transport_capacity(size(concentration,1),microbe_gas_count)
    real(r8), intent(in) :: surface_equilibrium_concentration(microbe_gas_count)
    real(r8), intent(in) :: surface_conductance(microbe_gas_count)
    real(r8), intent(in) :: aerenchyma_equilibrium_concentration(size(concentration,1),microbe_gas_count)
    real(r8), intent(in) :: aerenchyma_exchange_rate(size(concentration,1),microbe_gas_count)
    real(r8), intent(in) :: aerenchyma_minimum_emission_concentration( &
         size(concentration,1),microbe_gas_count)
    logical, intent(in) :: aerenchyma_allows_influx(microbe_gas_count)
    real(r8), intent(in) :: ch4_ebullition_threshold(size(concentration,1))
    real(r8), intent(in) :: ch4_ebullition_activation(size(concentration,1))
    real(r8), intent(in) :: dt
    real(r8), intent(out) :: interface_flux(0:size(concentration,1),microbe_gas_count)
    real(r8), intent(out) :: aerenchyma_flux(size(concentration,1),microbe_gas_count)
    real(r8), intent(out) :: ch4_ebullition_loss(size(concentration,1))
    real(r8), intent(out) :: tendency(size(concentration,1),microbe_gas_count)
    real(r8), intent(out) :: surface_diffusive_flux(microbe_gas_count)
    real(r8), intent(out) :: surface_aerenchyma_flux(microbe_gas_count)
    real(r8), intent(out) :: surface_ebullition_flux
    real(r8), intent(out) :: total_surface_flux(microbe_gas_count)
    real(r8) :: diffusion_tendency(size(concentration,1))
    real(r8) :: diffused_concentration(size(concentration,1),microbe_gas_count)
    real(r8) :: aerenchyma_tendency(size(concentration,1))
    real(r8) :: ebullition_tendency(size(concentration,1))
    real(r8) :: donor_scale(size(concentration,1),microbe_gas_count)
    real(r8) :: outgoing_flux
    integer :: gas, j, number_of_layers

    number_of_layers = size(concentration,1)
    interface_flux = 0._r8
    aerenchyma_flux = 0._r8
    ch4_ebullition_loss = 0._r8
    tendency = 0._r8
    surface_diffusive_flux = 0._r8
    surface_aerenchyma_flux = 0._r8
    surface_ebullition_flux = 0._r8
    total_surface_flux = 0._r8
    if (number_of_layers == 0 .or. dt <= 0._r8) return

    do gas = 1, microbe_gas_count
       call computeMicrobeImplicitVerticalDiffusion(concentration(:,gas), layer_thickness, &
            effective_diffusivity(:,gas), transport_capacity(:,gas), &
            surface_equilibrium_concentration(gas), surface_conductance(gas), dt, &
            diffused_concentration(:,gas), interface_flux(:,gas), diffusion_tendency, &
            surface_diffusive_flux(gas))
       tendency(:,gas) = diffusion_tendency
       call computeMicrobeAerenchymaTransport(diffused_concentration(:,gas), &
            aerenchyma_equilibrium_concentration(:,gas), layer_thickness, &
            aerenchyma_exchange_rate(:,gas), dt, aerenchyma_flux(:,gas), &
            aerenchyma_tendency, surface_aerenchyma_flux(gas), &
            apply_donor_limit=.false., minimum_emission_concentration= &
            aerenchyma_minimum_emission_concentration(:,gas), &
            allow_influx=aerenchyma_allows_influx(gas))
    end do
    call computeMicrobeMethaneEbullition(diffused_concentration(:,microbe_gas_ch4), &
         ch4_ebullition_threshold, ch4_ebullition_activation, layer_thickness, dt, &
         ch4_ebullition_loss, ebullition_tendency, surface_ebullition_flux)

    ! Diffusion is already a positive, conservative implicit transaction.
    ! Jointly limit the remaining simultaneous aerenchyma and ebullition
    ! outflows against the post-diffusion inventory.
    donor_scale = 1._r8
    do gas = 1, microbe_gas_count
       do j = 1, number_of_layers
          outgoing_flux = max(0._r8, aerenchyma_flux(j,gas))
          if (gas == microbe_gas_ch4) outgoing_flux = outgoing_flux + &
               ch4_ebullition_loss(j)
          if (outgoing_flux > 0._r8 .and. layer_thickness(j) > 0._r8) then
             donor_scale(j,gas) = clampUnitInterval(max(0._r8, &
                  diffused_concentration(j,gas)) * &
                  layer_thickness(j) / (dt * outgoing_flux))
          end if
       end do
    end do

    do gas = 1, microbe_gas_count
       do j = 1, number_of_layers
          if (aerenchyma_flux(j,gas) > 0._r8) then
             aerenchyma_flux(j,gas) = aerenchyma_flux(j,gas) * donor_scale(j,gas)
          end if
       end do
    end do
    do j = 1, number_of_layers
       ch4_ebullition_loss(j) = ch4_ebullition_loss(j) * &
            donor_scale(j,microbe_gas_ch4)
    end do

    surface_ebullition_flux = 0._r8
    do gas = 1, microbe_gas_count
       surface_diffusive_flux(gas) = -interface_flux(0,gas)
       surface_aerenchyma_flux(gas) = 0._r8
       do j = 1, number_of_layers
          if (layer_thickness(j) > 0._r8) then
             tendency(j,gas) = tendency(j,gas) - aerenchyma_flux(j,gas)
             surface_aerenchyma_flux(gas) = surface_aerenchyma_flux(gas) + &
                  aerenchyma_flux(j,gas) * layer_thickness(j)
          end if
       end do
    end do
    do j = 1, number_of_layers
       tendency(j,microbe_gas_ch4) = tendency(j,microbe_gas_ch4) - &
            ch4_ebullition_loss(j)
       surface_ebullition_flux = surface_ebullition_flux + &
            ch4_ebullition_loss(j) * max(0._r8, layer_thickness(j))
    end do
    total_surface_flux = surface_diffusive_flux + surface_aerenchyma_flux
    total_surface_flux(microbe_gas_ch4) = total_surface_flux(microbe_gas_ch4) + &
         surface_ebullition_flux
  end subroutine computeMicrobeGasTransport

  pure real(r8) function microbeMethaneEffectiveAqueousDiffusivity(reference_diffusivity, &
       temperature, hydrologic_transport_scalar, parameters) result(effective_diffusivity)
    real(r8), intent(in) :: reference_diffusivity
    real(r8), intent(in) :: temperature
    real(r8), intent(in) :: hydrologic_transport_scalar
    type(MicrobeMethaneParamsType), intent(in) :: parameters
    real(r8) :: temperature_ratio

    temperature_ratio = max(0._r8, temperature) / parameters%aqueous_diffusion_t_ref
    effective_diffusivity = max(0._r8, reference_diffusivity) * &
         parameters%aqueous_gas_diffusion_multiplier * &
         temperature_ratio ** parameters%aqueous_diffusion_temperature_exponent * &
         clampUnitInterval(hydrologic_transport_scalar)
  end function microbeMethaneEffectiveAqueousDiffusivity

  pure subroutine computeMicrobeTopounitLateralDiffusion(concentration, storage_weight, &
       layer_top_elevation, layer_bottom_elevation, effective_diffusivity, &
       edge_source, edge_target, edge_distance, dt, updated_concentration, &
       conservation_residual, valid)
    ! Simultaneous Fickian exchange over an arbitrary graph of topounits.
    ! storage_weight converts a layer concentration to inventory per gridcell
    ! area; the ELM adapter supplies topounit and saturated-area fractions as
    ! well as layer thickness. Absolute elevations allow unequal or vertically
    ! offset soil grids to exchange only over their physical overlap.
    real(r8), intent(in) :: concentration(:,:)
    real(r8), intent(in) :: storage_weight(size(concentration,1),size(concentration,2))
    real(r8), intent(in) :: layer_top_elevation(size(concentration,1),size(concentration,2))
    real(r8), intent(in) :: layer_bottom_elevation(size(concentration,1),size(concentration,2))
    real(r8), intent(in) :: effective_diffusivity(size(concentration,1),size(concentration,2))
    integer, intent(in) :: edge_source(:), edge_target(size(edge_source))
    real(r8), intent(in) :: edge_distance(size(edge_source))
    real(r8), intent(in) :: dt
    real(r8), intent(out) :: updated_concentration(size(concentration,1),size(concentration,2))
    real(r8), intent(out) :: conservation_residual
    logical, intent(out) :: valid
    real(r8) :: outgoing_inventory(size(concentration,1),size(concentration,2))
    real(r8) :: inventory_change(size(concentration,1),size(concentration,2))
    real(r8) :: donor_scale(size(concentration,1),size(concentration,2))
    real(r8) :: overlap, pair_diffusivity, transfer, total_inventory
    real(r8) :: source_area_weight, target_area_weight, edge_area_weight
    real(r8) :: source_top, source_bottom, target_top, target_bottom
    integer :: edge, source, target, source_layer, target_layer

    updated_concentration = concentration
    conservation_residual = 0._r8
    valid = .true.
    if (dt <= 0._r8 .or. size(concentration,1) == 0 .or. &
         size(concentration,2) == 0 .or. size(edge_source) == 0) return

    outgoing_inventory = 0._r8
    inventory_change = 0._r8
    donor_scale = 1._r8

    ! First pass: total every simultaneous outflow from each layer so a node
    ! connected to several neighbors receives one common donor limiter.
    do edge = 1, size(edge_source)
       source = edge_source(edge)
       target = edge_target(edge)
       if (source < 1 .or. source > size(concentration,1) .or. &
            target < 1 .or. target > size(concentration,1) .or. &
            source == target .or. edge_distance(edge) <= 0._r8) cycle
       do source_layer = 1, size(concentration,2)
          if (storage_weight(source,source_layer) <= 0._r8 .or. &
               effective_diffusivity(source,source_layer) <= 0._r8) cycle
          source_top = layer_top_elevation(source,source_layer)
          source_bottom = layer_bottom_elevation(source,source_layer)
          do target_layer = 1, size(concentration,2)
             if (storage_weight(target,target_layer) <= 0._r8 .or. &
                  effective_diffusivity(target,target_layer) <= 0._r8) cycle
             target_top = layer_top_elevation(target,target_layer)
             target_bottom = layer_bottom_elevation(target,target_layer)
             overlap = max(0._r8, min(source_top, target_top) - &
                  max(source_bottom, target_bottom))
             if (overlap <= 0._r8) cycle
             ! Convert the per-interface-area Fickian flux to inventory per
             ! gridcell area using the smaller participating horizontal
             ! footprint. storage_weight is horizontal area fraction times
             ! layer thickness. Without this factor, transfer into a tiny
             ! saturated partition is divided by near-zero storage and can
             ! create an unphysical concentration spike.
             source_area_weight = storage_weight(source,source_layer) / &
                  max(source_top - source_bottom, tiny(1._r8))
             target_area_weight = storage_weight(target,target_layer) / &
                  max(target_top - target_bottom, tiny(1._r8))
             edge_area_weight = min(source_area_weight, target_area_weight)
             pair_diffusivity = 2._r8 / &
                  (1._r8/effective_diffusivity(source,source_layer) + &
                  1._r8/effective_diffusivity(target,target_layer))
             transfer = dt * edge_area_weight * overlap * pair_diffusivity / &
                  edge_distance(edge) * &
                  (concentration(source,source_layer) - &
                  concentration(target,target_layer))
             if (transfer > 0._r8) then
                outgoing_inventory(source,source_layer) = &
                     outgoing_inventory(source,source_layer) + transfer
             else if (transfer < 0._r8) then
                outgoing_inventory(target,target_layer) = &
                     outgoing_inventory(target,target_layer) - transfer
             end if
          end do
       end do
    end do

    do source = 1, size(concentration,1)
       do source_layer = 1, size(concentration,2)
          if (outgoing_inventory(source,source_layer) > 0._r8) then
             donor_scale(source,source_layer) = clampUnitInterval( &
                  max(0._r8, concentration(source,source_layer)) * &
                  storage_weight(source,source_layer) / &
                  outgoing_inventory(source,source_layer))
          end if
       end do
    end do

    ! Second pass: apply each pair transfer with its donor's common scale.
    do edge = 1, size(edge_source)
       source = edge_source(edge)
       target = edge_target(edge)
       if (source < 1 .or. source > size(concentration,1) .or. &
            target < 1 .or. target > size(concentration,1) .or. &
            source == target .or. edge_distance(edge) <= 0._r8) cycle
       do source_layer = 1, size(concentration,2)
          if (storage_weight(source,source_layer) <= 0._r8 .or. &
               effective_diffusivity(source,source_layer) <= 0._r8) cycle
          source_top = layer_top_elevation(source,source_layer)
          source_bottom = layer_bottom_elevation(source,source_layer)
          do target_layer = 1, size(concentration,2)
             if (storage_weight(target,target_layer) <= 0._r8 .or. &
                  effective_diffusivity(target,target_layer) <= 0._r8) cycle
             target_top = layer_top_elevation(target,target_layer)
             target_bottom = layer_bottom_elevation(target,target_layer)
             overlap = max(0._r8, min(source_top, target_top) - &
                  max(source_bottom, target_bottom))
             if (overlap <= 0._r8) cycle
             source_area_weight = storage_weight(source,source_layer) / &
                  max(source_top - source_bottom, tiny(1._r8))
             target_area_weight = storage_weight(target,target_layer) / &
                  max(target_top - target_bottom, tiny(1._r8))
             edge_area_weight = min(source_area_weight, target_area_weight)
             pair_diffusivity = 2._r8 / &
                  (1._r8/effective_diffusivity(source,source_layer) + &
                  1._r8/effective_diffusivity(target,target_layer))
             transfer = dt * edge_area_weight * overlap * pair_diffusivity / &
                  edge_distance(edge) * &
                  (concentration(source,source_layer) - &
                  concentration(target,target_layer))
             if (transfer > 0._r8) then
                transfer = transfer * donor_scale(source,source_layer)
             else if (transfer < 0._r8) then
                transfer = transfer * donor_scale(target,target_layer)
             end if
             inventory_change(source,source_layer) = &
                  inventory_change(source,source_layer) - transfer
             inventory_change(target,target_layer) = &
                  inventory_change(target,target_layer) + transfer
          end do
       end do
    end do

    total_inventory = 0._r8
    do source = 1, size(concentration,1)
       do source_layer = 1, size(concentration,2)
          if (storage_weight(source,source_layer) > 0._r8) then
             updated_concentration(source,source_layer) = max(0._r8, &
                  (concentration(source,source_layer) * &
                  storage_weight(source,source_layer) + &
                  inventory_change(source,source_layer)) / &
                  storage_weight(source,source_layer))
             total_inventory = total_inventory + abs(concentration(source,source_layer)) * &
                  storage_weight(source,source_layer)
          end if
       end do
    end do
    conservation_residual = sum((updated_concentration - concentration) * storage_weight)
    valid = all(updated_concentration >= 0._r8) .and. &
         abs(conservation_residual) <= 1.e-11_r8 * max(1._r8, total_inventory)
  end subroutine computeMicrobeTopounitLateralDiffusion

  pure real(r8) function microbeMethaneTransportResidual(tendency, layer_thickness, &
       surface_flux) result(residual)
    real(r8), intent(in) :: tendency(:)
    real(r8), intent(in) :: layer_thickness(size(tendency))
    real(r8), intent(in) :: surface_flux

    residual = sum(tendency * layer_thickness) + surface_flux
  end function microbeMethaneTransportResidual

  pure real(r8) function clampUnitInterval(value) result(clamped)
    real(r8), intent(in) :: value
    clamped = min(1._r8, max(0._r8, value))
  end function clampUnitInterval

end module MicrobeGasTransportMod
