#!/usr/bin/env python3
"""Independent scalar/profile oracle for Phase 3 Step 3 transport."""

from __future__ import annotations

import math
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


def observed_dom_profile_restoration(
    dom_c: Sequence[float],
    dom_n: Sequence[float],
    dom_p: Sequence[float],
    target_dom_c: Sequence[float],
    cn_dom: float,
    cp_dom: float,
    relaxation_timescale_days: float,
    dt: float,
) -> tuple[list[float], list[float], list[float], list[float]]:
    fraction = 1.0 - math.exp(-dt / (86400.0 * relaxation_timescale_days))
    updated_c = [old + fraction * (target - old)
                 for old, target in zip(dom_c, target_dom_c)]
    updated_n = [old + fraction * (target / cn_dom - old)
                 for old, target in zip(dom_n, target_dom_c)]
    updated_p = [old + fraction * (target / cp_dom - old)
                 for old, target in zip(dom_p, target_dom_c)]
    source = [(new - old) / dt for new, old in zip(updated_c, dom_c)]
    return updated_c, updated_n, updated_p, source


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


def implicit_vertical_diffusion(
    concentration: Sequence[float],
    layer_thickness: Sequence[float],
    effective_diffusivity: Sequence[float],
    transport_capacity: Sequence[float],
    atmospheric_mobile_concentration: float,
    surface_conductance: float,
    dt: float,
) -> tuple[list[float], list[float], list[float], float]:
    """Backward-Euler diffusion of a bulk inventory on a mobile-gas basis."""
    number_of_layers = len(concentration)
    updated_concentration = list(concentration)
    interface_flux = [0.0] * (number_of_layers + 1)
    tendency = [0.0] * number_of_layers
    if number_of_layers == 0 or dt <= 0.0:
        return updated_concentration, interface_flux, tendency, 0.0

    tiny = 2.2250738585072014e-308
    capacity = [max(value, tiny) for value in transport_capacity]
    conductance = [0.0] * (number_of_layers + 1)
    conductance[0] = max(0.0, surface_conductance)
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
            conductance[layer + 1] = 1.0 / resistance

    lower = [0.0] * number_of_layers
    diagonal = list(capacity)
    upper = [0.0] * number_of_layers
    rhs = list(concentration)
    for layer in range(number_of_layers):
        if layer_thickness[layer] > 0.0:
            lower[layer] = -dt * conductance[layer] / layer_thickness[layer]
            upper[layer] = -dt * conductance[layer + 1] / layer_thickness[layer]
            diagonal[layer] = capacity[layer] - lower[layer] - upper[layer]
    rhs[0] -= lower[0] * max(0.0, atmospheric_mobile_concentration)
    lower[0] = 0.0
    upper[-1] = 0.0

    cprime = [0.0] * number_of_layers
    dprime = [0.0] * number_of_layers
    denominator = max(diagonal[0], tiny)
    cprime[0] = upper[0] / denominator
    dprime[0] = rhs[0] / denominator
    for layer in range(1, number_of_layers):
        denominator = max(
            diagonal[layer] - lower[layer] * cprime[layer - 1], tiny
        )
        cprime[layer] = upper[layer] / denominator
        dprime[layer] = (
            rhs[layer] - lower[layer] * dprime[layer - 1]
        ) / denominator

    mobile_concentration = [0.0] * number_of_layers
    mobile_concentration[-1] = dprime[-1]
    for layer in range(number_of_layers - 2, -1, -1):
        mobile_concentration[layer] = (
            dprime[layer] - cprime[layer] * mobile_concentration[layer + 1]
        )

    interface_flux[0] = conductance[0] * (
        max(0.0, atmospheric_mobile_concentration) - mobile_concentration[0]
    )
    for layer in range(number_of_layers - 1):
        interface_flux[layer + 1] = conductance[layer + 1] * (
            mobile_concentration[layer] - mobile_concentration[layer + 1]
        )
    for layer in range(number_of_layers):
        if layer_thickness[layer] > 0.0:
            tendency[layer] = (
                interface_flux[layer] - interface_flux[layer + 1]
            ) / layer_thickness[layer]
            updated_concentration[layer] = (
                concentration[layer] + dt * tendency[layer]
            )
    surface_flux = -interface_flux[0]
    return updated_concentration, interface_flux, tendency, surface_flux


def dom_cnp_relaxation(
    dom_c: Sequence[float],
    dom_n: Sequence[float],
    dom_p: Sequence[float],
    layer_thickness: Sequence[float],
    layer_relaxation_rate: Sequence[float],
    dt: float,
) -> tuple[list[float], list[float], list[float], tuple[float, float, float]]:
    """Conservative implicit DOM C/N/P relaxation used by the parity test."""

    def relax(concentration: Sequence[float]) -> list[float]:
        number_of_layers = len(concentration)
        if number_of_layers == 0 or dt <= 0.0:
            return list(concentration)
        conductance = [0.0] * (number_of_layers + 1)
        for layer in range(number_of_layers - 1):
            upper_depth = layer_thickness[layer]
            lower_depth = layer_thickness[layer + 1]
            harmonic_storage_depth = (
                2.0 * upper_depth * lower_depth / (upper_depth + lower_depth)
            )
            conductance[layer + 1] = (
                max(0.0, layer_relaxation_rate[layer + 1])
                * harmonic_storage_depth
            )

        lower = [0.0] * number_of_layers
        diagonal = [0.0] * number_of_layers
        upper = [0.0] * number_of_layers
        for layer in range(number_of_layers):
            lower[layer] = -dt * conductance[layer] / layer_thickness[layer]
            upper[layer] = -dt * conductance[layer + 1] / layer_thickness[layer]
            diagonal[layer] = 1.0 - lower[layer] - upper[layer]
        lower[0] = 0.0
        upper[-1] = 0.0

        cprime = [0.0] * number_of_layers
        dprime = [0.0] * number_of_layers
        cprime[0] = upper[0] / diagonal[0]
        dprime[0] = concentration[0] / diagonal[0]
        for layer in range(1, number_of_layers):
            denominator = diagonal[layer] - lower[layer] * cprime[layer - 1]
            cprime[layer] = upper[layer] / denominator
            dprime[layer] = (
                concentration[layer] - lower[layer] * dprime[layer - 1]
            ) / denominator
        updated = [0.0] * number_of_layers
        updated[-1] = dprime[-1]
        for layer in range(number_of_layers - 2, -1, -1):
            updated[layer] = dprime[layer] - cprime[layer] * updated[layer + 1]
        return updated

    updated_c = relax(dom_c)
    updated_n = relax(dom_n)
    updated_p = relax(dom_p)

    def inventory(values: Sequence[float]) -> float:
        return sum(value * depth for value, depth in zip(values, layer_thickness))

    residuals = (
        inventory(updated_c) - inventory(dom_c),
        inventory(updated_n) - inventory(dom_n),
        inventory(updated_p) - inventory(dom_p),
    )
    return updated_c, updated_n, updated_p, residuals


def nh4_dissolved_fraction(
    liquid_fraction: float,
    dry_bulk_density: float,
    partition_coefficient: float,
) -> float:
    """Equilibrium fraction of total NH4-N present in soil water."""
    theta = max(0.0, liquid_fraction)
    sorption_capacity = max(0.0, dry_bulk_density) * max(0.0, partition_coefficient)
    denominator = theta + sorption_capacity
    return theta / denominator if denominator > 0.0 else 0.0


def aqueous_tracer_transport(
    concentration: Sequence[float],
    layer_thickness: Sequence[float],
    liquid_fraction: Sequence[float],
    diffusion_conductivity: Sequence[float],
    water_flux: Sequence[float],
    mobile_fraction: float | Sequence[float],
    minimum_liquid_fraction: float,
    dt: float,
) -> tuple[list[float], list[float], list[float], list[float], float, float]:
    """Implicit porewater advection-diffusion of a bulk-soil inventory."""
    n = len(concentration)
    if isinstance(mobile_fraction, (int, float)):
        layer_mobile_fraction = [float(mobile_fraction)] * n
    else:
        layer_mobile_fraction = list(mobile_fraction)
        if len(layer_mobile_fraction) != n:
            raise ValueError("mobile_fraction must be scalar or have one value per layer")
    porewater_factor = [
        mobile / theta if theta >= minimum_liquid_fraction else 0.0
        for mobile, theta in zip(layer_mobile_fraction, liquid_fraction)
    ]
    conductance = [0.0] * (n + 1)
    for layer in range(n - 1):
        if (
            diffusion_conductivity[layer] > 0.0
            and diffusion_conductivity[layer + 1] > 0.0
            and porewater_factor[layer] > 0.0
            and porewater_factor[layer + 1] > 0.0
        ):
            resistance = (
                0.5 * layer_thickness[layer] / diffusion_conductivity[layer]
                + 0.5
                * layer_thickness[layer + 1]
                / diffusion_conductivity[layer + 1]
            )
            conductance[layer + 1] = 1.0 / resistance

    lower = [0.0] * n
    diagonal = [1.0] * n
    upper = [0.0] * n
    for interface in range(1, n):
        left = interface - 1
        right = interface
        q = water_flux[interface]
        left_coefficient = (max(q, 0.0) + conductance[interface]) * porewater_factor[left]
        right_coefficient = (min(q, 0.0) - conductance[interface]) * porewater_factor[right]
        diagonal[left] += dt * left_coefficient / layer_thickness[left]
        upper[left] += dt * right_coefficient / layer_thickness[left]
        lower[right] -= dt * left_coefficient / layer_thickness[right]
        diagonal[right] -= dt * right_coefficient / layer_thickness[right]
    if water_flux[0] < 0.0:
        diagonal[0] -= dt * water_flux[0] * porewater_factor[0] / layer_thickness[0]
    if water_flux[n] > 0.0:
        diagonal[-1] += dt * water_flux[n] * porewater_factor[-1] / layer_thickness[-1]

    cprime = [0.0] * n
    dprime = [0.0] * n
    cprime[0] = upper[0] / diagonal[0]
    dprime[0] = concentration[0] / diagonal[0]
    for layer in range(1, n):
        denominator = diagonal[layer] - lower[layer] * cprime[layer - 1]
        cprime[layer] = upper[layer] / denominator
        dprime[layer] = (
            concentration[layer] - lower[layer] * dprime[layer - 1]
        ) / denominator
    updated = [0.0] * n
    updated[-1] = dprime[-1]
    for layer in range(n - 2, -1, -1):
        updated[layer] = dprime[layer] - cprime[layer] * updated[layer + 1]

    advective = [0.0] * (n + 1)
    diffusive = [0.0] * (n + 1)
    if water_flux[0] < 0.0:
        advective[0] = water_flux[0] * porewater_factor[0] * updated[0]
    for interface in range(1, n):
        left = interface - 1
        right = interface
        donor = left if water_flux[interface] >= 0.0 else right
        advective[interface] = (
            water_flux[interface] * porewater_factor[donor] * updated[donor]
        )
        diffusive[interface] = conductance[interface] * (
            porewater_factor[left] * updated[left]
            - porewater_factor[right] * updated[right]
        )
    if water_flux[n] > 0.0:
        advective[n] = water_flux[n] * porewater_factor[-1] * updated[-1]
    tendency = [(new - old) / dt for new, old in zip(updated, concentration)]
    boundary_export = advective[-1] + diffusive[-1] - advective[0] - diffusive[0]
    initial_inventory = sum(c * dz for c, dz in zip(concentration, layer_thickness))
    final_inventory = sum(c * dz for c, dz in zip(updated, layer_thickness))
    residual = final_inventory - initial_inventory + dt * boundary_export
    return updated, advective, diffusive, tendency, boundary_export, residual


def topounit_lateral_diffusion(
    concentration: Sequence[Sequence[float]],
    storage_weight: Sequence[Sequence[float]],
    layer_top_elevation: Sequence[Sequence[float]],
    layer_bottom_elevation: Sequence[Sequence[float]],
    effective_diffusivity: Sequence[Sequence[float]],
    edges: Sequence[tuple[int, int, float]],
    dt: float,
) -> tuple[list[list[float]], float]:
    """Conservative simultaneous diffusion over a zero-based topounit graph."""
    node_count = len(concentration)
    layer_count = len(concentration[0]) if node_count else 0
    updated = [list(profile) for profile in concentration]
    if node_count == 0 or layer_count == 0 or dt <= 0.0 or not edges:
        return updated, 0.0

    outgoing = [[0.0] * layer_count for _ in range(node_count)]
    change = [[0.0] * layer_count for _ in range(node_count)]

    def transfers():
        for source, target, distance in edges:
            if (
                source < 0
                or source >= node_count
                or target < 0
                or target >= node_count
                or source == target
                or distance <= 0.0
            ):
                continue
            for source_layer in range(layer_count):
                if (
                    storage_weight[source][source_layer] <= 0.0
                    or effective_diffusivity[source][source_layer] <= 0.0
                ):
                    continue
                for target_layer in range(layer_count):
                    if (
                        storage_weight[target][target_layer] <= 0.0
                        or effective_diffusivity[target][target_layer] <= 0.0
                    ):
                        continue
                    overlap = max(
                        0.0,
                        min(
                            layer_top_elevation[source][source_layer],
                            layer_top_elevation[target][target_layer],
                        )
                        - max(
                            layer_bottom_elevation[source][source_layer],
                            layer_bottom_elevation[target][target_layer],
                        ),
                    )
                    if overlap <= 0.0:
                        continue
                    source_diffusivity = effective_diffusivity[source][source_layer]
                    target_diffusivity = effective_diffusivity[target][target_layer]
                    pair_diffusivity = 2.0 / (
                        1.0 / source_diffusivity + 1.0 / target_diffusivity
                    )
                    source_thickness = (
                        layer_top_elevation[source][source_layer]
                        - layer_bottom_elevation[source][source_layer]
                    )
                    target_thickness = (
                        layer_top_elevation[target][target_layer]
                        - layer_bottom_elevation[target][target_layer]
                    )
                    edge_area_weight = min(
                        storage_weight[source][source_layer] / source_thickness,
                        storage_weight[target][target_layer] / target_thickness,
                    )
                    transfer = (
                        dt
                        * edge_area_weight
                        * overlap
                        * pair_diffusivity
                        / distance
                        * (
                            concentration[source][source_layer]
                            - concentration[target][target_layer]
                        )
                    )
                    yield source, source_layer, target, target_layer, transfer

    for source, source_layer, target, target_layer, transfer in transfers():
        if transfer > 0.0:
            outgoing[source][source_layer] += transfer
        elif transfer < 0.0:
            outgoing[target][target_layer] -= transfer

    donor_scale = [[1.0] * layer_count for _ in range(node_count)]
    for node in range(node_count):
        for layer in range(layer_count):
            if outgoing[node][layer] > 0.0:
                donor_scale[node][layer] = _clamp(
                    max(0.0, concentration[node][layer])
                    * storage_weight[node][layer]
                    / outgoing[node][layer]
                )

    for source, source_layer, target, target_layer, transfer in transfers():
        if transfer > 0.0:
            transfer *= donor_scale[source][source_layer]
        elif transfer < 0.0:
            transfer *= donor_scale[target][target_layer]
        change[source][source_layer] -= transfer
        change[target][target_layer] += transfer

    for node in range(node_count):
        for layer in range(layer_count):
            weight = storage_weight[node][layer]
            if weight > 0.0:
                updated[node][layer] = max(
                    0.0,
                    (concentration[node][layer] * weight + change[node][layer])
                    / weight,
                )
    residual = sum(
        (updated[node][layer] - concentration[node][layer])
        * storage_weight[node][layer]
        for node in range(node_count)
        for layer in range(layer_count)
    )
    return updated, residual


def aerenchyma_transport(
    concentration: Sequence[float],
    atmospheric_equivalent_concentration: Sequence[float],
    layer_thickness: Sequence[float],
    layer_exchange_rate: Sequence[float],
    dt: float,
    apply_donor_limit: bool = True,
    minimum_emission_concentration: Sequence[float] | None = None,
    allow_influx: bool = True,
) -> tuple[list[float], list[float], float]:
    flux_to_atmosphere = [0.0] * len(concentration)
    tendency = [0.0] * len(concentration)
    if dt <= 0.0:
        return flux_to_atmosphere, tendency, 0.0
    for layer, value in enumerate(concentration):
        exchange_target = atmospheric_equivalent_concentration[layer]
        if not allow_influx and minimum_emission_concentration is not None:
            exchange_target = max(
                exchange_target, max(0.0, minimum_emission_concentration[layer])
            )
        flux = max(0.0, layer_exchange_rate[layer]) * (
            value - exchange_target
        )
        if not allow_influx:
            flux = max(0.0, flux)
        equilibrium_flux = (value - exchange_target) / dt
        if flux > 0.0:
            flux = min(flux, max(0.0, equilibrium_flux))
        elif flux < 0.0:
            flux = max(flux, min(0.0, equilibrium_flux))
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
    transport_capacity: Sequence[Sequence[float]],
    surface_equilibrium_concentration: Sequence[float],
    surface_conductance: Sequence[float],
    aerenchyma_equilibrium_concentration: Sequence[Sequence[float]],
    aerenchyma_exchange_rate: Sequence[Sequence[float]],
    aerenchyma_minimum_emission_concentration: Sequence[Sequence[float]],
    aerenchyma_allows_influx: Sequence[bool],
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
    """Run implicit diffusion, then jointly limit plant and bubble losses."""
    number_of_layers = len(concentration)
    interfaces = [[0.0] * GAS_COUNT for _ in range(number_of_layers + 1)]
    aerenchyma_flux = [[0.0] * GAS_COUNT for _ in range(number_of_layers)]
    ch4_ebullition_loss = [0.0] * number_of_layers
    tendency = [[0.0] * GAS_COUNT for _ in range(number_of_layers)]
    diffused_concentration = [[0.0] * GAS_COUNT for _ in range(number_of_layers)]
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
        gas_capacity = [layer[gas] for layer in transport_capacity]
        gas_aerenchyma_equilibrium = [
            layer[gas] for layer in aerenchyma_equilibrium_concentration
        ]
        gas_aerenchyma_rate = [layer[gas] for layer in aerenchyma_exchange_rate]
        gas_aerenchyma_minimum_emission = [
            layer[gas] for layer in aerenchyma_minimum_emission_concentration
        ]
        (
            gas_diffused_concentration,
            gas_interfaces,
            diffusion_tendency,
            _surface_diffusion,
        ) = implicit_vertical_diffusion(
            gas_concentration,
            layer_thickness,
            gas_diffusivity,
            gas_capacity,
            surface_equilibrium_concentration[gas],
            surface_conductance[gas],
            dt,
        )
        gas_aerenchyma_flux, _aerenchyma_tendency, _surface_aerenchyma = (
            aerenchyma_transport(
                gas_diffused_concentration,
                gas_aerenchyma_equilibrium,
                layer_thickness,
                gas_aerenchyma_rate,
                dt,
                apply_donor_limit=False,
                minimum_emission_concentration=gas_aerenchyma_minimum_emission,
                allow_influx=aerenchyma_allows_influx[gas],
            )
        )
        for interface, value in enumerate(gas_interfaces):
            interfaces[interface][gas] = value
        for layer, value in enumerate(gas_aerenchyma_flux):
            aerenchyma_flux[layer][gas] = value
        for layer, value in enumerate(gas_diffused_concentration):
            diffused_concentration[layer][gas] = value
            tendency[layer][gas] = diffusion_tendency[layer]

    ch4_concentration = [layer[GAS_CH4] for layer in diffused_concentration]
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
            outgoing = max(0.0, aerenchyma_flux[layer][gas])
            if gas == GAS_CH4:
                outgoing += ch4_ebullition_loss[layer]
            if outgoing > 0.0 and layer_thickness[layer] > 0.0:
                donor_scale[layer][gas] = _clamp(
                    max(0.0, diffused_concentration[layer][gas])
                    * layer_thickness[layer]
                    / (dt * outgoing)
                )

    for gas in range(GAS_COUNT):
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
                tendency[layer][gas] -= aerenchyma_flux[layer][gas]
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


def saturated_dom_macrodispersion_conductivity(
    liquid_fraction: float,
    relative_liquid_saturation: float,
    thawed_fraction: float,
    macrodispersion: float,
    saturation_threshold: float,
) -> float:
    """Layer theta*D contribution for thawed, nearly saturated DOM mixing."""
    saturation = _clamp(relative_liquid_saturation)
    threshold = _clamp(saturation_threshold)
    if threshold < 1.0:
        saturation_scalar = _clamp((saturation - threshold) / (1.0 - threshold))
    else:
        saturation_scalar = 1.0 if saturation >= 1.0 else 0.0
    return (
        max(0.0, liquid_fraction)
        * _clamp(thawed_fraction)
        * max(0.0, macrodispersion)
        * saturation_scalar
    )


def zwt_saturated_layer_fraction(
    layer_top: float,
    layer_bottom: float,
    water_table_depth: float,
) -> float:
    """Fraction of a soil layer below the connected water table."""
    if layer_bottom <= layer_top:
        raise ValueError("layer_bottom must be deeper than layer_top")
    water_table = max(0.0, water_table_depth)
    return _clamp(
        (layer_bottom - max(layer_top, water_table))
        / (layer_bottom - layer_top)
    )


def legacy_fickian_gas_diffusivity(
    gas_index: int,
    temperature: float,
    parameters: Mapping[str, float],
) -> float:
    """Effective coefficient in the active CLM-Microbe vertical equation."""
    source_fick_d_w = (1.49e-5, 2.10e-5, 1.92e-5, 4.50e-5)
    if gas_index < 0 or gas_index >= len(source_fick_d_w):
        return 0.0
    reference_temperature = parameters["microbe_methane_aqueous_diffusion_t_ref"]
    exponent = parameters["microbe_methane_aqueous_diffusion_temperature_exponent"]
    multiplier = parameters["microbe_methane_aqueous_gas_diffusion_multiplier"]
    return (
        source_fick_d_w[gas_index]
        * 1.0e-3
        * multiplier
        * (max(0.0, temperature) / reference_temperature) ** exponent
    )


def transport_residual(
    tendency: Sequence[float], layer_thickness: Sequence[float], surface_flux: float
) -> float:
    return sum(
        rate * thickness for rate, thickness in zip(tendency, layer_thickness)
    ) + surface_flux
