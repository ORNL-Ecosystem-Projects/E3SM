#!/usr/bin/env python3
"""Independent scalar oracle for the Phase 3 Step 2 reaction kernel."""

from __future__ import annotations

from dataclasses import dataclass, fields
from math import log10
from typing import Mapping


CATOMW = 12.011
SECSPDAY = 86400.0


@dataclass
class State:
    dom_c: float = 0.0
    acetate_c: float = 0.0
    acetate_methanogen_c: float = 0.0
    h2_methanogen_c: float = 0.0
    aerobic_methanotroph_c: float = 0.0
    anaerobic_methanotroph_c: float = 0.0
    conc_ch4: float = 0.0
    conc_o2: float = 0.0
    conc_co2: float = 0.0
    conc_h2: float = 0.0


@dataclass
class Environment:
    soil_temperature: float = 273.15
    soil_ph: float = 7.0
    dom_fermentation_scalar: float = 0.0
    aerobic_acetate_oxidation_scalar: float = 0.0


@dataclass
class Rates:
    dom_to_acetate_c: float = 0.0
    acetogenesis_c: float = 0.0
    acetoclastic_methanogenesis_c: float = 0.0
    hydrogenotrophic_methanogenesis_c: float = 0.0
    aerobic_acetate_oxidation_c: float = 0.0
    aerobic_methane_oxidation_c: float = 0.0
    anaerobic_methane_oxidation_c: float = 0.0
    acetate_methanogen_mortality_c: float = 0.0
    h2_methanogen_mortality_c: float = 0.0
    aerobic_methanotroph_mortality_c: float = 0.0
    anaerobic_methanotroph_mortality_c: float = 0.0
    effective_soil_ph: float = 0.0
    ph_response: float = 0.0


@dataclass
class Tendencies(State):
    pass


def _p(parameters: Mapping[str, float], short_name: str) -> float:
    return parameters[f"microbe_methane_{short_name}"]


def _clamp(value: float) -> float:
    return min(1.0, max(0.0, value))


def ph_response(ph: float, ph_min: float, ph_opt: float, ph_max: float) -> float:
    if ph <= ph_min or ph >= ph_max:
        return 0.0
    numerator = (ph - ph_min) * (ph - ph_max)
    denominator = numerator - (ph - ph_opt) ** 2
    return 0.0 if denominator == 0.0 else _clamp(numerator / denominator)


def q10_response(q10: float, temperature: float, reference_temperature: float) -> float:
    return q10 ** ((temperature - reference_temperature) / 10.0)


def _monod(substrate: float, half_saturation: float) -> float:
    substrate = max(0.0, substrate)
    return substrate / (substrate + half_saturation)


def _supply_scale(available_rate: float, demand_rate: float) -> float:
    return 1.0 if demand_rate <= 0.0 else _clamp(max(0.0, available_rate) / demand_rate)


def potential_rates(
    state: State, environment: Environment, parameters: Mapping[str, float]
) -> Rates:
    rates = Rates()
    dom_mmol_c = max(0.0, state.dom_c) / CATOMW * 1000.0
    acetate_mmol_c = max(0.0, state.acetate_c) / CATOMW * 1000.0
    ch4_mmol = max(0.0, state.conc_ch4) * 1000.0
    o2_mmol = max(0.0, state.conc_o2) * 1000.0
    co2_mmol = max(0.0, state.conc_co2) * 1000.0
    h2_mmol = max(0.0, state.conc_h2) * 1000.0

    effective_ph = environment.soil_ph
    if effective_ph > _p(parameters, "acetate_ph_trigger") and acetate_mmol_c > 0.0:
        effective_ph = -log10(
            10.0 ** (-effective_ph)
            + _p(parameters, "acidification_coefficient") * acetate_mmol_c
        )
    rates.effective_soil_ph = effective_ph
    rates.ph_response = ph_response(
        effective_ph,
        _p(parameters, "ph_min"),
        _p(parameters, "ph_opt"),
        _p(parameters, "ph_max"),
    )
    acetate_feedback = _p(parameters, "acetate_feedback_half_saturation") / (
        acetate_mmol_c + _p(parameters, "acetate_feedback_half_saturation")
    )
    co2_inhibition = max(
        0.0,
        1.0 - min(1.0, co2_mmol / _p(parameters, "h2_methanogenesis_co2_inhibition_scale")),
    )
    o2_inhibition = max(
        0.0, 1.0 - min(1.0, o2_mmol / _p(parameters, "aom_o2_inhibition_scale"))
    )

    rates.dom_to_acetate_c = (
        2.0
        / 3.0
        * _p(parameters, "acetate_prod_max")
        * _monod(dom_mmol_c, _p(parameters, "k_acetate"))
        * q10_response(
            _p(parameters, "dom_to_acetate_q10"),
            environment.soil_temperature,
            _p(parameters, "reaction_t_ref"),
        )
        * rates.ph_response
        * acetate_feedback
        * _clamp(environment.dom_fermentation_scalar)
        * 1.0e-3
        / SECSPDAY
    )
    rates.acetogenesis_c = (
        _p(parameters, "acetogenesis_max")
        * _monod(h2_mmol, _p(parameters, "k_acetogenesis_h2"))
        * _monod(co2_mmol, _p(parameters, "k_acetogenesis_co2"))
        * q10_response(
            _p(parameters, "acetogenesis_q10"),
            environment.soil_temperature,
            _p(parameters, "reaction_t_ref"),
        )
        * rates.ph_response
        * 1.0e-3
        / SECSPDAY
    )
    rates.hydrogenotrophic_methanogenesis_c = (
        _p(parameters, "h2_methanogen_growth_rate")
        / _p(parameters, "h2_methanogen_yield")
        / SECSPDAY
        * max(0.0, state.h2_methanogen_c)
        / CATOMW
        * _monod(h2_mmol, _p(parameters, "k_h2_methanogenesis_h2"))
        * _monod(co2_mmol, _p(parameters, "k_h2_methanogenesis_co2"))
        * q10_response(
            _p(parameters, "h2_methanogenesis_q10"),
            environment.soil_temperature,
            _p(parameters, "reaction_t_ref"),
        )
        * rates.ph_response
        * co2_inhibition
    )
    rates.acetoclastic_methanogenesis_c = (
        _p(parameters, "acetate_methanogen_growth_rate")
        / _p(parameters, "acetate_methanogen_yield")
        / SECSPDAY
        * max(0.0, state.acetate_methanogen_c)
        / CATOMW
        * _monod(acetate_mmol_c, _p(parameters, "k_acetoclastic_methanogenesis_acetate"))
        * q10_response(
            _p(parameters, "acetoclastic_methanogenesis_q10"),
            environment.soil_temperature,
            _p(parameters, "reaction_t_ref"),
        )
        * rates.ph_response
    )
    rates.aerobic_acetate_oxidation_c = (
        _p(parameters, "aerobic_acetate_oxidation_rate")
        / SECSPDAY
        * max(0.0, state.acetate_c)
        / CATOMW
        * _monod(o2_mmol, _p(parameters, "k_acetate_prod_o2"))
        * _clamp(environment.aerobic_acetate_oxidation_scalar)
    )
    rates.aerobic_methane_oxidation_c = (
        _p(parameters, "aerobic_methanotroph_growth_rate")
        / _p(parameters, "aerobic_methanotroph_yield")
        / SECSPDAY
        * max(0.0, state.aerobic_methanotroph_c)
        / CATOMW
        * _monod(ch4_mmol, _p(parameters, "k_aerobic_oxidation_ch4"))
        * _monod(o2_mmol, _p(parameters, "k_aerobic_oxidation_o2"))
        * q10_response(
            _p(parameters, "aerobic_oxidation_q10"),
            environment.soil_temperature,
            _p(parameters, "reaction_t_ref"),
        )
        * rates.ph_response
    )
    rates.anaerobic_methane_oxidation_c = (
        _p(parameters, "anaerobic_methanotroph_growth_rate")
        / _p(parameters, "anaerobic_methanotroph_yield")
        / SECSPDAY
        * max(0.0, state.anaerobic_methanotroph_c)
        / CATOMW
        * _monod(ch4_mmol, _p(parameters, "k_anaerobic_oxidation_ch4"))
        * q10_response(
            _p(parameters, "anaerobic_oxidation_q10"),
            environment.soil_temperature,
            _p(parameters, "aom_t_ref"),
        )
        * rates.ph_response
        * o2_inhibition
    )
    for rate_name, parameter_name, state_name in (
        ("acetate_methanogen_mortality_c", "acetate_methanogen_death_rate", "acetate_methanogen_c"),
        ("h2_methanogen_mortality_c", "h2_methanogen_death_rate", "h2_methanogen_c"),
        ("aerobic_methanotroph_mortality_c", "aerobic_methanotroph_death_rate", "aerobic_methanotroph_c"),
        ("anaerobic_methanotroph_mortality_c", "anaerobic_methanotroph_death_rate", "anaerobic_methanotroph_c"),
    ):
        setattr(
            rates,
            rate_name,
            _p(parameters, parameter_name)
            / SECSPDAY
            * max(0.0, getattr(state, state_name))
            / CATOMW
            * rates.ph_response,
        )
    return rates


def limit_rates(state: State, parameters: Mapping[str, float], dt: float, rates: Rates) -> Rates:
    if dt <= 0.0:
        return Rates()
    rates.dom_to_acetate_c = min(
        rates.dom_to_acetate_c, max(0.0, state.dom_c) / CATOMW / (1.5 * dt)
    )
    available_h2 = max(0.0, state.conc_h2) / dt + rates.dom_to_acetate_c / 6.0
    available_co2 = max(0.0, state.conc_co2) / dt + 0.5 * rates.dom_to_acetate_c
    h2_demand = 2.0 * rates.acetogenesis_c + 4.0 * rates.hydrogenotrophic_methanogenesis_c
    co2_demand = rates.acetogenesis_c + (
        1.0 + _p(parameters, "h2_methanogen_yield")
    ) * rates.hydrogenotrophic_methanogenesis_c
    scale = min(_supply_scale(available_h2, h2_demand), _supply_scale(available_co2, co2_demand))
    rates.acetogenesis_c *= scale
    rates.hydrogenotrophic_methanogenesis_c *= scale

    available_acetate = (
        max(0.0, state.acetate_c) / CATOMW / dt
        + rates.dom_to_acetate_c
        + rates.acetogenesis_c
    )
    scale = _supply_scale(
        available_acetate,
        rates.acetoclastic_methanogenesis_c + rates.aerobic_acetate_oxidation_c,
    )
    rates.acetoclastic_methanogenesis_c *= scale
    rates.aerobic_acetate_oxidation_c *= scale

    available_ch4 = (
        max(0.0, state.conc_ch4) / dt
        + _p(parameters, "acetoclastic_methanogenesis_ch4_yield")
        * (1.0 - _p(parameters, "acetate_methanogen_yield"))
        * rates.acetoclastic_methanogenesis_c
        + rates.hydrogenotrophic_methanogenesis_c
    )
    scale = _supply_scale(
        available_ch4,
        rates.aerobic_methane_oxidation_c + rates.anaerobic_methane_oxidation_c,
    )
    rates.aerobic_methane_oxidation_c *= scale
    rates.anaerobic_methane_oxidation_c *= scale

    available_o2 = max(0.0, state.conc_o2) / dt
    o2_demand = (
        _p(parameters, "aerobic_decomp_o2_c_ratio") * rates.aerobic_acetate_oxidation_c
        + _p(parameters, "aerobic_oxidation_o2_ch4_ratio")
        * rates.aerobic_methane_oxidation_c
    )
    scale = _supply_scale(available_o2, o2_demand)
    rates.aerobic_acetate_oxidation_c *= scale
    rates.aerobic_methane_oxidation_c *= scale

    for rate_name, state_name in (
        ("acetate_methanogen_mortality_c", "acetate_methanogen_c"),
        ("h2_methanogen_mortality_c", "h2_methanogen_c"),
        ("aerobic_methanotroph_mortality_c", "aerobic_methanotroph_c"),
        ("anaerobic_methanotroph_mortality_c", "anaerobic_methanotroph_c"),
    ):
        setattr(
            rates,
            rate_name,
            min(getattr(rates, rate_name), max(0.0, getattr(state, state_name)) / CATOMW / dt),
        )
    return rates


def assemble(parameters: Mapping[str, float], rates: Rates) -> Tendencies:
    y_acetate = _p(parameters, "acetate_methanogen_yield")
    y_h2 = _p(parameters, "h2_methanogen_yield")
    y_aerobic = _p(parameters, "aerobic_methanotroph_yield")
    y_anaerobic = _p(parameters, "anaerobic_methanotroph_yield")
    acetate_nonbiomass = (1.0 - y_acetate) * rates.acetoclastic_methanogenesis_c
    acetate_ch4 = _p(parameters, "acetoclastic_methanogenesis_ch4_yield") * acetate_nonbiomass
    acetate_co2 = (
        1.0 - _p(parameters, "acetoclastic_methanogenesis_ch4_yield")
    ) * acetate_nonbiomass
    return Tendencies(
        dom_c=CATOMW
        * (
            -1.5 * rates.dom_to_acetate_c
            + rates.acetate_methanogen_mortality_c
            + rates.h2_methanogen_mortality_c
            + rates.aerobic_methanotroph_mortality_c
            + rates.anaerobic_methanotroph_mortality_c
        ),
        acetate_c=CATOMW
        * (
            rates.dom_to_acetate_c
            + rates.acetogenesis_c
            - rates.acetoclastic_methanogenesis_c
            - rates.aerobic_acetate_oxidation_c
        ),
        acetate_methanogen_c=CATOMW
        * (y_acetate * rates.acetoclastic_methanogenesis_c - rates.acetate_methanogen_mortality_c),
        h2_methanogen_c=CATOMW
        * (y_h2 * rates.hydrogenotrophic_methanogenesis_c - rates.h2_methanogen_mortality_c),
        aerobic_methanotroph_c=CATOMW
        * (y_aerobic * rates.aerobic_methane_oxidation_c - rates.aerobic_methanotroph_mortality_c),
        anaerobic_methanotroph_c=CATOMW
        * (y_anaerobic * rates.anaerobic_methane_oxidation_c - rates.anaerobic_methanotroph_mortality_c),
        conc_ch4=acetate_ch4
        + rates.hydrogenotrophic_methanogenesis_c
        - rates.aerobic_methane_oxidation_c
        - rates.anaerobic_methane_oxidation_c,
        conc_o2=-_p(parameters, "aerobic_decomp_o2_c_ratio")
        * rates.aerobic_acetate_oxidation_c
        - _p(parameters, "aerobic_oxidation_o2_ch4_ratio")
        * rates.aerobic_methane_oxidation_c,
        conc_co2=0.5 * rates.dom_to_acetate_c
        - rates.acetogenesis_c
        - (1.0 + y_h2) * rates.hydrogenotrophic_methanogenesis_c
        + acetate_co2
        + rates.aerobic_acetate_oxidation_c
        + (1.0 - y_aerobic) * rates.aerobic_methane_oxidation_c
        + (1.0 - y_anaerobic) * rates.anaerobic_methane_oxidation_c,
        conc_h2=rates.dom_to_acetate_c / 6.0
        - 2.0 * rates.acetogenesis_c
        - 4.0 * rates.hydrogenotrophic_methanogenesis_c,
    )


def compute(
    state: State, environment: Environment, parameters: Mapping[str, float], dt: float
) -> tuple[Rates, Tendencies]:
    rates = limit_rates(state, parameters, dt, potential_rates(state, environment, parameters))
    return rates, assemble(parameters, rates)


def carbon_residual(tendencies: Tendencies) -> float:
    organic = sum(
        getattr(tendencies, item.name)
        for item in fields(Tendencies)
        if item.name.endswith("_c")
    )
    return organic / CATOMW + tendencies.conc_ch4 + tendencies.conc_co2
