module MicrobeGasTransportMod

  ! Conservative, side-effect-free transport kernels for the revised
  ! microbial methane backend. Concentrations may be gas (mol m-3 soil) or
  ! acetate (g C m-3 soil); flux and tendency units follow the input unit.

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
  public :: computeMicrobeAerenchymaTransport
  public :: computeMicrobeMethaneEbullition
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
    ! The combined gas driver disables this per-path limiter and applies one
    ! joint limiter across diffusion, aerenchyma, and ebullition instead.
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

  pure subroutine computeMicrobeAerenchymaTransport(concentration, &
       atmospheric_equivalent_concentration, layer_thickness, layer_exchange_rate, &
       dt, flux_to_atmosphere, tendency, surface_flux, apply_donor_limit)
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
    integer :: j
    real(r8) :: equilibrium_flux
    logical :: limit_inventory

    flux_to_atmosphere = 0._r8
    tendency = 0._r8
    surface_flux = 0._r8
    if (dt <= 0._r8) return
    limit_inventory = .true.
    if (present(apply_donor_limit)) limit_inventory = apply_donor_limit

    do j = 1, size(concentration)
       flux_to_atmosphere(j) = max(0._r8, layer_exchange_rate(j)) * &
            (concentration(j) - atmospheric_equivalent_concentration(j))
       ! Explicit exchange must not cross its reservoir equilibrium in one
       ! timestep. This bound applies in both directions independently of the
       ! joint soil-inventory limiter used by the combined transport driver.
       equilibrium_flux = (concentration(j) - &
            atmospheric_equivalent_concentration(j)) / dt
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
       effective_diffusivity, surface_equilibrium_concentration, surface_conductance, &
       aerenchyma_equilibrium_concentration, aerenchyma_exchange_rate, &
       ch4_ebullition_threshold, ch4_ebullition_activation, dt, interface_flux, &
       aerenchyma_flux, ch4_ebullition_loss, tendency, surface_diffusive_flux, &
       surface_aerenchyma_flux, surface_ebullition_flux, total_surface_flux)
    ! Gas dimension order is CH4, O2, CO2, H2. All concentrations use
    ! mol gas m-3 soil; layer tendencies use mol gas m-3 s-1 and surface
    ! fluxes use mol gas m-2 s-1, positive to the atmosphere.
    real(r8), intent(in) :: concentration(:,:)
    real(r8), intent(in) :: layer_thickness(size(concentration,1))
    real(r8), intent(in) :: effective_diffusivity(size(concentration,1),microbe_gas_count)
    real(r8), intent(in) :: surface_equilibrium_concentration(microbe_gas_count)
    real(r8), intent(in) :: surface_conductance(microbe_gas_count)
    real(r8), intent(in) :: aerenchyma_equilibrium_concentration(size(concentration,1),microbe_gas_count)
    real(r8), intent(in) :: aerenchyma_exchange_rate(size(concentration,1),microbe_gas_count)
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
       call computeMicrobeVerticalDiffusion(concentration(:,gas), layer_thickness, &
            effective_diffusivity(:,gas), surface_equilibrium_concentration(gas), &
            surface_conductance(gas), dt, interface_flux(:,gas), diffusion_tendency, &
            surface_diffusive_flux(gas), apply_donor_limit=.false.)
       call computeMicrobeAerenchymaTransport(concentration(:,gas), &
            aerenchyma_equilibrium_concentration(:,gas), layer_thickness, &
            aerenchyma_exchange_rate(:,gas), dt, aerenchyma_flux(:,gas), &
            aerenchyma_tendency, surface_aerenchyma_flux(gas), &
            apply_donor_limit=.false.)
    end do
    call computeMicrobeMethaneEbullition(concentration(:,microbe_gas_ch4), &
         ch4_ebullition_threshold, ch4_ebullition_activation, layer_thickness, dt, &
         ch4_ebullition_loss, ebullition_tendency, surface_ebullition_flux)

    ! Jointly limit raw diffusion, aerenchyma loss, and ebullition with one
    ! donor scale. This preserves their relative rates while remaining safe
    ! without using simultaneous incoming flux as immediately available.
    donor_scale = 1._r8
    do gas = 1, microbe_gas_count
       do j = 1, number_of_layers
          outgoing_flux = max(0._r8, -interface_flux(j-1,gas)) + &
               max(0._r8, interface_flux(j,gas)) + &
               max(0._r8, aerenchyma_flux(j,gas))
          if (gas == microbe_gas_ch4) outgoing_flux = outgoing_flux + &
               ch4_ebullition_loss(j)
          if (outgoing_flux > 0._r8 .and. layer_thickness(j) > 0._r8) then
             donor_scale(j,gas) = clampUnitInterval(max(0._r8, concentration(j,gas)) * &
                  layer_thickness(j) / (dt * outgoing_flux))
          end if
       end do
    end do

    do gas = 1, microbe_gas_count
       if (interface_flux(0,gas) < 0._r8) then
          interface_flux(0,gas) = interface_flux(0,gas) * donor_scale(1,gas)
       end if
       do j = 1, number_of_layers - 1
          if (interface_flux(j,gas) > 0._r8) then
             interface_flux(j,gas) = interface_flux(j,gas) * donor_scale(j,gas)
          else if (interface_flux(j,gas) < 0._r8) then
             interface_flux(j,gas) = interface_flux(j,gas) * donor_scale(j+1,gas)
          end if
       end do
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
             tendency(j,gas) = (interface_flux(j-1,gas) - &
                  interface_flux(j,gas)) / layer_thickness(j) - &
                  aerenchyma_flux(j,gas)
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
