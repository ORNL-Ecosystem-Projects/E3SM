#!/usr/bin/env python3
"""Independent transaction oracle for Phase 3 Step 4 state updates."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Mapping

from reaction_oracle import CATOMW, Environment, Rates, State, Tendencies, assemble, compute


@dataclass
class NutrientTendencies:
    dom_n: float = 0.0
    dom_p: float = 0.0
    mineral_n: float = 0.0
    mineral_p: float = 0.0


@dataclass
class ReactionTransaction:
    unsaturated_state: State
    saturated_state: State
    unsaturated_rates: Rates
    saturated_rates: Rates
    unsaturated_tendencies: Tendencies
    saturated_tendencies: Tendencies
    nutrient_tendencies: NutrientTendencies
    dom_c: float
    dom_n: float
    dom_p: float
    mineral_n: float
    mineral_p: float
    carbon_residual: float
    nitrogen_residual: float
    phosphorus_residual: float
    valid: bool


def additional_carbon_density(state: State) -> float:
    """Carbon outside the authoritative ELM DOM pool, in g C m-3."""
    return (
        state.acetate_c
        + state.acetate_methanogen_c
        + state.h2_methanogen_c
        + state.aerobic_methanotroph_c
        + state.anaerobic_methanotroph_c
        + CATOMW * (state.conc_ch4 + state.conc_co2)
    )


def column_additional_carbon(
    saturated_fraction: float,
    unsaturated_state: list[State],
    saturated_state: list[State],
    layer_thickness: list[float],
) -> float:
    fraction = _clamp(saturated_fraction)
    return sum(
        dz
        * (
            (1.0 - fraction) * additional_carbon_density(unsaturated)
            + fraction * additional_carbon_density(saturated)
        )
        for unsaturated, saturated, dz in zip(
            unsaturated_state, saturated_state, layer_thickness
        )
    )


def surface_carbon_flux(total_surface_flux: list[float]) -> float:
    return CATOMW * (total_surface_flux[0] + total_surface_flux[2])


def ch4_surface_flux_kgc(total_surface_flux: list[float]) -> float:
    return CATOMW * total_surface_flux[0] / 1000.0


def co2_correction(total_surface_flux: list[float]) -> float:
    return CATOMW * total_surface_flux[2]


def advance_reaction_layer(
    dom_c: float,
    dom_n: float,
    dom_p: float,
    mineral_n: float,
    mineral_p: float,
    saturated_fraction: float,
    unsaturated_state: State,
    saturated_state: State,
    unsaturated_environment: Environment,
    saturated_environment: Environment,
    parameters: Mapping[str, float],
    cn_dom: float,
    cp_dom: float,
    dt: float,
) -> ReactionTransaction:
    fraction = _clamp(saturated_fraction)
    available_dom_c = min(
        max(0.0, dom_c), max(0.0, dom_n) * cn_dom, max(0.0, dom_p) * cp_dom
    )
    unsaturated_work = replace(unsaturated_state, dom_c=available_dom_c)
    saturated_work = replace(saturated_state, dom_c=available_dom_c)
    unsaturated_rates, unsaturated_tendencies = compute(
        unsaturated_work, unsaturated_environment, parameters, dt
    )
    saturated_rates, saturated_tendencies = compute(
        saturated_work, saturated_environment, parameters, dt
    )
    _limit_mortality_for_nutrients(
        fraction,
        mineral_n,
        mineral_p,
        cn_dom,
        cp_dom,
        dt,
        unsaturated_rates,
        saturated_rates,
    )
    unsaturated_tendencies = assemble(parameters, unsaturated_rates)
    saturated_tendencies = assemble(parameters, saturated_rates)
    dom_c_tendency = (
        (1.0 - fraction) * unsaturated_tendencies.dom_c
        + fraction * saturated_tendencies.dom_c
    )
    nutrients = NutrientTendencies(
        dom_n=dom_c_tendency / cn_dom,
        dom_p=dom_c_tendency / cp_dom,
        mineral_n=-dom_c_tendency / cn_dom,
        mineral_p=-dom_c_tendency / cp_dom,
    )
    new_dom_c = dom_c + dt * dom_c_tendency
    new_dom_n = dom_n + dt * nutrients.dom_n
    new_dom_p = dom_p + dt * nutrients.dom_p
    new_mineral_n = mineral_n + dt * nutrients.mineral_n
    new_mineral_p = mineral_p + dt * nutrients.mineral_p
    new_unsaturated = _apply(unsaturated_state, unsaturated_tendencies, new_dom_c, dt)
    new_saturated = _apply(saturated_state, saturated_tendencies, new_dom_c, dt)
    initial_carbon = dom_c + (1.0 - fraction) * additional_carbon_density(
        unsaturated_state
    ) + fraction * additional_carbon_density(saturated_state)
    final_carbon = new_dom_c + (1.0 - fraction) * additional_carbon_density(
        new_unsaturated
    ) + fraction * additional_carbon_density(new_saturated)
    carbon_residual = final_carbon - initial_carbon
    initial_nitrogen = max(0.0, dom_n) + mineral_n
    final_nitrogen = new_dom_n + new_mineral_n
    nitrogen_residual = final_nitrogen - initial_nitrogen
    initial_phosphorus = max(0.0, dom_p) + mineral_p
    final_phosphorus = new_dom_p + new_mineral_p
    phosphorus_residual = final_phosphorus - initial_phosphorus
    floor = parameters["microbe_methane_mfg_biomass_min"]
    return ReactionTransaction(
        unsaturated_state=new_unsaturated,
        saturated_state=new_saturated,
        unsaturated_rates=unsaturated_rates,
        saturated_rates=saturated_rates,
        unsaturated_tendencies=unsaturated_tendencies,
        saturated_tendencies=saturated_tendencies,
        nutrient_tendencies=nutrients,
        dom_c=new_dom_c,
        dom_n=new_dom_n,
        dom_p=new_dom_p,
        mineral_n=new_mineral_n,
        mineral_p=new_mineral_p,
        carbon_residual=carbon_residual,
        nitrogen_residual=nitrogen_residual,
        phosphorus_residual=phosphorus_residual,
        valid=(
            _state_valid(new_unsaturated, floor)
            and _state_valid(new_saturated, floor)
            and min(new_dom_c, new_dom_n, new_dom_p) >= -1.0e-12
            and new_mineral_n >= min(0.0, mineral_n) - 1.0e-12
            and new_mineral_p >= min(0.0, mineral_p) - 1.0e-12
            and _residual_closed(carbon_residual, initial_carbon, final_carbon)
            and _residual_closed(
                nitrogen_residual, initial_nitrogen, final_nitrogen
            )
            and _residual_closed(
                phosphorus_residual, initial_phosphorus, final_phosphorus
            )
        ),
    )


def _limit_mortality_for_nutrients(
    fraction: float,
    mineral_n: float,
    mineral_p: float,
    cn_dom: float,
    cp_dom: float,
    dt: float,
    unsaturated_rates: Rates,
    saturated_rates: Rates,
) -> None:
    fermentation_loss = 1.5 * CATOMW * (
        (1.0 - fraction) * unsaturated_rates.dom_to_acetate_c
        + fraction * saturated_rates.dom_to_acetate_c
    )
    mortality_input = CATOMW * (
        (1.0 - fraction) * _mortality(unsaturated_rates)
        + fraction * _mortality(saturated_rates)
    )
    permitted_dom_gain = min(
        max(0.0, mineral_n) * cn_dom / dt,
        max(0.0, mineral_p) * cp_dom / dt,
    )
    if mortality_input > fermentation_loss + permitted_dom_gain:
        scale = _clamp((fermentation_loss + permitted_dom_gain) / mortality_input)
        for rates in (unsaturated_rates, saturated_rates):
            for name in (
                "acetate_methanogen_mortality_c",
                "h2_methanogen_mortality_c",
                "aerobic_methanotroph_mortality_c",
                "anaerobic_methanotroph_mortality_c",
            ):
                setattr(rates, name, getattr(rates, name) * scale)


def _mortality(rates: Rates) -> float:
    return sum(
        getattr(rates, name)
        for name in (
            "acetate_methanogen_mortality_c",
            "h2_methanogen_mortality_c",
            "aerobic_methanotroph_mortality_c",
            "anaerobic_methanotroph_mortality_c",
        )
    )


def _apply(state: State, tendency: Tendencies, dom_c: float, dt: float) -> State:
    return State(
        **{
            name: dom_c
            if name == "dom_c"
            else getattr(state, name) + dt * getattr(tendency, name)
            for name in State.__dataclass_fields__
        }
    )


def _state_valid(state: State, floor: float) -> bool:
    return (
        min(
            state.dom_c,
            state.acetate_c,
            state.conc_ch4,
            state.conc_o2,
            state.conc_co2,
            state.conc_h2,
        )
        >= -1.0e-12
        and min(
            state.acetate_methanogen_c,
            state.h2_methanogen_c,
            state.aerobic_methanotroph_c,
            state.anaerobic_methanotroph_c,
        )
        >= floor - 1.0e-12
    )


def _clamp(value: float) -> float:
    return min(1.0, max(0.0, value))


def _residual_closed(residual: float, initial: float, final: float) -> bool:
    return abs(residual) <= 1.0e-12 * max(1.0, abs(initial), abs(final))
