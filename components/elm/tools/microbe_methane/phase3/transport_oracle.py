#!/usr/bin/env python3
"""Independent scalar/profile oracle for Phase 3 Step 3 transport."""

from __future__ import annotations

from typing import Mapping, Sequence


GAS_CH4 = 0
GAS_COUNT = 4


def _clamp(value: float) -> float:
    return min(1.0, max(0.0, value))


def repartition(
    old_saturated_fraction: float,
    new_saturated_fraction: float,
    unsaturated_concentration: float,
    saturated_concentration: float,
) -> tuple[float, float]:
    old_fraction = _clamp(old_saturated_fraction)
    new_fraction = _clamp(new_saturated_fraction)
    new_unsaturated = unsaturated_concentration
    new_saturated = saturated_concentration
    if new_fraction > old_fraction:
        transferred = new_fraction - old_fraction
        new_saturated = (
            old_fraction * saturated_concentration
            + transferred * unsaturated_concentration
        ) / new_fraction
    elif new_fraction < old_fraction:
        transferred = old_fraction - new_fraction
        new_unsaturated = (
            (1.0 - old_fraction) * unsaturated_concentration
            + transferred * saturated_concentration
        ) / (1.0 - new_fraction)
    return new_unsaturated, new_saturated


def bulk_concentration(
    saturated_fraction: float,
    unsaturated_concentration: float,
    saturated_concentration: float,
) -> float:
    fraction = _clamp(saturated_fraction)
    return (
        (1.0 - fraction) * unsaturated_concentration
        + fraction * saturated_concentration
    )


def vertical_diffusion(
    concentration: Sequence[float],
    layer_thickness: Sequence[float],
    effective_diffusivity: Sequence[float],
    atmospheric_equivalent_concentration: float,
    surface_conductance: float,
    dt: float,
    apply_donor_limit: bool = True,
) -> tuple[list[float], list[float], float]:
    number_of_layers = len(concentration)
    interface_flux = [0.0] * (number_of_layers + 1)
    tendency = [0.0] * number_of_layers
    if number_of_layers == 0 or dt <= 0.0:
        return interface_flux, tendency, 0.0

    interface_flux[0] = max(0.0, surface_conductance) * (
        atmospheric_equivalent_concentration - concentration[0]
    )
    for layer in range(number_of_layers - 1):
        if (
            effective_diffusivity[layer] > 0.0
            and effective_diffusivity[layer + 1] > 0.0
            and layer_thickness[layer] > 0.0
            and layer_thickness[layer + 1] > 0.0
        ):
            resistance = (
                0.5 * layer_thickness[layer] / effective_diffusivity[layer]
                + 0.5
                * layer_thickness[layer + 1]
                / effective_diffusivity[layer + 1]
            )
            interface_flux[layer + 1] = (
                concentration[layer] - concentration[layer + 1]
            ) / resistance

    if apply_donor_limit:
        donor_scale = [1.0] * number_of_layers
        for layer in range(number_of_layers):
            outgoing_flux = max(0.0, -interface_flux[layer]) + max(
                0.0, interface_flux[layer + 1]
            )
            if outgoing_flux > 0.0 and layer_thickness[layer] > 0.0:
                donor_scale[layer] = _clamp(
                    max(0.0, concentration[layer])
                    * layer_thickness[layer]
                    / (dt * outgoing_flux)
                )

        if interface_flux[0] < 0.0:
            interface_flux[0] *= donor_scale[0]
        for interface in range(1, number_of_layers):
            if interface_flux[interface] > 0.0:
                interface_flux[interface] *= donor_scale[interface - 1]
            elif interface_flux[interface] < 0.0:
                interface_flux[interface] *= donor_scale[interface]

    for layer in range(number_of_layers):
        if layer_thickness[layer] > 0.0:
            tendency[layer] = (
                interface_flux[layer] - interface_flux[layer + 1]
            ) / layer_thickness[layer]
    surface_flux = -interface_flux[0]
    return interface_flux, tendency, surface_flux


def aerenchyma_transport(
    concentration: Sequence[float],
    atmospheric_equivalent_concentration: Sequence[float],
    layer_thickness: Sequence[float],
    layer_exchange_rate: Sequence[float],
    dt: float,
    apply_donor_limit: bool = True,
) -> tuple[list[float], list[float], float]:
    flux_to_atmosphere = [0.0] * len(concentration)
    tendency = [0.0] * len(concentration)
    if dt <= 0.0:
        return flux_to_atmosphere, tendency, 0.0
    for layer, value in enumerate(concentration):
        flux = max(0.0, layer_exchange_rate[layer]) * (
            value - atmospheric_equivalent_concentration[layer]
        )
        if apply_donor_limit and flux > 0.0:
            flux = min(flux, max(0.0, value) / dt)
        flux_to_atmosphere[layer] = flux
        tendency[layer] = -flux
    surface_flux = sum(
        flux * max(0.0, thickness)
        for flux, thickness in zip(flux_to_atmosphere, layer_thickness)
    )
    return flux_to_atmosphere, tendency, surface_flux


def methane_ebullition(
    concentration: Sequence[float],
    threshold: Sequence[float],
    activation_fraction: Sequence[float],
    layer_thickness: Sequence[float],
    dt: float,
) -> tuple[list[float], list[float], float]:
    loss_rate = [0.0] * len(concentration)
    tendency = [0.0] * len(concentration)
    if dt <= 0.0:
        return loss_rate, tendency, 0.0
    for layer, value in enumerate(concentration):
        loss = (
            _clamp(activation_fraction[layer])
            * max(0.0, value - max(0.0, threshold[layer]))
            / dt
        )
        loss_rate[layer] = loss
        tendency[layer] = -loss
    surface_flux = sum(
        loss * max(0.0, thickness)
        for loss, thickness in zip(loss_rate, layer_thickness)
    )
    return loss_rate, tendency, surface_flux


def gas_transport(
    concentration: Sequence[Sequence[float]],
    layer_thickness: Sequence[float],
    effective_diffusivity: Sequence[Sequence[float]],
    surface_equilibrium_concentration: Sequence[float],
    surface_conductance: Sequence[float],
    aerenchyma_equilibrium_concentration: Sequence[Sequence[float]],
    aerenchyma_exchange_rate: Sequence[Sequence[float]],
    ch4_ebullition_threshold: Sequence[float],
    ch4_ebullition_activation: Sequence[float],
    dt: float,
) -> tuple[
    list[list[float]],
    list[list[float]],
    list[float],
    list[list[float]],
    list[float],
]:
    """Run the combined four-gas transport interface with joint donor limiting."""
    number_of_layers = len(concentration)
    interfaces = [[0.0] * GAS_COUNT for _ in range(number_of_layers + 1)]
    aerenchyma_flux = [[0.0] * GAS_COUNT for _ in range(number_of_layers)]
    ch4_ebullition_loss = [0.0] * number_of_layers
    tendency = [[0.0] * GAS_COUNT for _ in range(number_of_layers)]
    total_surface_flux = [0.0] * GAS_COUNT
    if number_of_layers == 0 or dt <= 0.0:
        return (
            interfaces,
            aerenchyma_flux,
            ch4_ebullition_loss,
            tendency,
            total_surface_flux,
        )

    for gas in range(GAS_COUNT):
        gas_concentration = [layer[gas] for layer in concentration]
        gas_diffusivity = [layer[gas] for layer in effective_diffusivity]
        gas_aerenchyma_equilibrium = [
            layer[gas] for layer in aerenchyma_equilibrium_concentration
        ]
        gas_aerenchyma_rate = [layer[gas] for layer in aerenchyma_exchange_rate]
        gas_interfaces, _diffusion_tendency, _surface_diffusion = vertical_diffusion(
            gas_concentration,
            layer_thickness,
            gas_diffusivity,
            surface_equilibrium_concentration[gas],
            surface_conductance[gas],
            dt,
            apply_donor_limit=False,
        )
        gas_aerenchyma_flux, _aerenchyma_tendency, _surface_aerenchyma = (
            aerenchyma_transport(
                gas_concentration,
                gas_aerenchyma_equilibrium,
                layer_thickness,
                gas_aerenchyma_rate,
                dt,
                apply_donor_limit=False,
            )
        )
        for interface, value in enumerate(gas_interfaces):
            interfaces[interface][gas] = value
        for layer, value in enumerate(gas_aerenchyma_flux):
            aerenchyma_flux[layer][gas] = value

    ch4_concentration = [layer[GAS_CH4] for layer in concentration]
    ch4_ebullition_loss, _ebullition_tendency, _surface_ebullition = methane_ebullition(
        ch4_concentration,
        ch4_ebullition_threshold,
        ch4_ebullition_activation,
        layer_thickness,
        dt,
    )

    donor_scale = [[1.0] * GAS_COUNT for _ in range(number_of_layers)]
    for gas in range(GAS_COUNT):
        for layer in range(number_of_layers):
            outgoing = max(0.0, -interfaces[layer][gas]) + max(
                0.0, interfaces[layer + 1][gas]
            )
            outgoing += max(0.0, aerenchyma_flux[layer][gas])
            if gas == GAS_CH4:
                outgoing += ch4_ebullition_loss[layer]
            if outgoing > 0.0 and layer_thickness[layer] > 0.0:
                donor_scale[layer][gas] = _clamp(
                    max(0.0, concentration[layer][gas])
                    * layer_thickness[layer]
                    / (dt * outgoing)
                )

    for gas in range(GAS_COUNT):
        if interfaces[0][gas] < 0.0:
            interfaces[0][gas] *= donor_scale[0][gas]
        for interface in range(1, number_of_layers):
            if interfaces[interface][gas] > 0.0:
                interfaces[interface][gas] *= donor_scale[interface - 1][gas]
            elif interfaces[interface][gas] < 0.0:
                interfaces[interface][gas] *= donor_scale[interface][gas]
        for layer in range(number_of_layers):
            if aerenchyma_flux[layer][gas] > 0.0:
                aerenchyma_flux[layer][gas] *= donor_scale[layer][gas]
    for layer in range(number_of_layers):
        ch4_ebullition_loss[layer] *= donor_scale[layer][GAS_CH4]

    for gas in range(GAS_COUNT):
        surface_diffusion = -interfaces[0][gas]
        surface_aerenchyma = 0.0
        for layer in range(number_of_layers):
            if layer_thickness[layer] > 0.0:
                tendency[layer][gas] = (
                    interfaces[layer][gas] - interfaces[layer + 1][gas]
                ) / layer_thickness[layer] - aerenchyma_flux[layer][gas]
                surface_aerenchyma += (
                    aerenchyma_flux[layer][gas] * layer_thickness[layer]
                )
        total_surface_flux[gas] = surface_diffusion + surface_aerenchyma
    surface_ebullition = 0.0
    for layer in range(number_of_layers):
        tendency[layer][GAS_CH4] -= ch4_ebullition_loss[layer]
        surface_ebullition += ch4_ebullition_loss[layer] * max(
            0.0, layer_thickness[layer]
        )
    total_surface_flux[GAS_CH4] += surface_ebullition

    return (
        interfaces,
        aerenchyma_flux,
        ch4_ebullition_loss,
        tendency,
        total_surface_flux,
    )


def effective_aqueous_diffusivity(
    reference_diffusivity: float,
    temperature: float,
    hydrologic_transport_scalar: float,
    parameters: Mapping[str, float],
) -> float:
    reference_temperature = parameters["microbe_methane_aqueous_diffusion_t_ref"]
    exponent = parameters["microbe_methane_aqueous_diffusion_temperature_exponent"]
    multiplier = parameters["microbe_methane_aqueous_gas_diffusion_multiplier"]
    return (
        max(0.0, reference_diffusivity)
        * multiplier
        * (max(0.0, temperature) / reference_temperature) ** exponent
        * _clamp(hydrologic_transport_scalar)
    )


def transport_residual(
    tendency: Sequence[float], layer_thickness: Sequence[float], surface_flux: float
) -> float:
    return sum(
        rate * thickness for rate, thickness in zip(tendency, layer_thickness)
    ) + surface_flux
