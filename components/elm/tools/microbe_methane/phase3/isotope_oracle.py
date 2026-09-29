#!/usr/bin/env python3
"""Passive carbon-tracer routing through accepted microbial methane rates."""

from __future__ import annotations

from dataclasses import replace
from typing import Mapping

from reaction_oracle import CATOMW, Rates, State
from state_update_oracle import ReactionTransaction, additional_carbon_density


def _fraction(tracer: float, bulk: float) -> float:
    return 0.0 if bulk <= 0.0 else min(1.0, max(0.0, tracer) / bulk)


def _advance_partition(
    bulk_dom: float,
    bulk: State,
    tracer_dom: float,
    tracer: State,
    rates: Rates,
    parameters: Mapping[str, float],
    dt: float,
) -> tuple[float, State]:
    b = replace(bulk)
    t = replace(tracer, dom_c=0.0, conc_o2=0.0, conc_h2=0.0)
    dom_delta = 0.0
    initial_biomass_ratios = [
        _fraction(tracer.acetate_methanogen_c, bulk.acetate_methanogen_c),
        _fraction(tracer.h2_methanogen_c, bulk.h2_methanogen_c),
        _fraction(tracer.aerobic_methanotroph_c, bulk.aerobic_methanotroph_c),
        _fraction(tracer.anaerobic_methanotroph_c, bulk.anaerobic_methanotroph_c),
    ]

    # DOM -> acetate + CO2.
    extent = dt * rates.dom_to_acetate_c
    transfer = 1.5 * CATOMW * extent * _fraction(tracer_dom, bulk_dom)
    dom_delta -= transfer
    b.acetate_c += CATOMW * extent
    t.acetate_c += transfer / 1.5
    b.conc_co2 += 0.5 * extent
    t.conc_co2 += transfer / (3.0 * CATOMW)

    # CO2 -> acetate, CH4, and hydrogenotrophic biomass.
    acetogenesis = dt * rates.acetogenesis_c
    h2_meth = dt * rates.hydrogenotrophic_methanogenesis_c
    y = parameters["microbe_methane_h2_methanogen_yield"]
    ratio = _fraction(t.conc_co2, b.conc_co2)
    donor = acetogenesis + (1.0 + y) * h2_meth
    b.conc_co2 -= donor
    t.conc_co2 -= donor * ratio
    b.acetate_c += CATOMW * acetogenesis
    t.acetate_c += CATOMW * acetogenesis * ratio
    b.conc_ch4 += h2_meth
    t.conc_ch4 += h2_meth * ratio
    b.h2_methanogen_c += CATOMW * y * h2_meth
    t.h2_methanogen_c += CATOMW * y * h2_meth * ratio

    # Acetate -> biomass, CH4, and CO2; or aerobic CO2.
    acetoclastic = dt * rates.acetoclastic_methanogenesis_c
    aerobic_acetate = dt * rates.aerobic_acetate_oxidation_c
    ratio = _fraction(t.acetate_c, b.acetate_c)
    donor = CATOMW * (acetoclastic + aerobic_acetate)
    b.acetate_c -= donor
    t.acetate_c -= donor * ratio
    y = parameters["microbe_methane_acetate_methanogen_yield"]
    fch4 = parameters["microbe_methane_acetoclastic_methanogenesis_ch4_yield"]
    b.acetate_methanogen_c += CATOMW * y * acetoclastic
    t.acetate_methanogen_c += CATOMW * y * acetoclastic * ratio
    b.conc_ch4 += fch4 * (1.0 - y) * acetoclastic
    t.conc_ch4 += fch4 * (1.0 - y) * acetoclastic * ratio
    co2 = (1.0 - fch4) * (1.0 - y) * acetoclastic + aerobic_acetate
    b.conc_co2 += co2
    t.conc_co2 += co2 * ratio

    # CH4 -> methanotroph biomass + CO2.
    aerobic_ch4 = dt * rates.aerobic_methane_oxidation_c
    anaerobic_ch4 = dt * rates.anaerobic_methane_oxidation_c
    ratio = _fraction(t.conc_ch4, b.conc_ch4)
    donor = aerobic_ch4 + anaerobic_ch4
    b.conc_ch4 -= donor
    t.conc_ch4 -= donor * ratio
    for extent, yield_name, biomass_name in (
        (aerobic_ch4, "aerobic_methanotroph_yield", "aerobic_methanotroph_c"),
        (anaerobic_ch4, "anaerobic_methanotroph_yield", "anaerobic_methanotroph_c"),
    ):
        y = parameters[f"microbe_methane_{yield_name}"]
        setattr(b, biomass_name, getattr(b, biomass_name) + CATOMW * y * extent)
        setattr(t, biomass_name, getattr(t, biomass_name) + CATOMW * y * extent * ratio)
        b.conc_co2 += (1.0 - y) * extent
        t.conc_co2 += (1.0 - y) * extent * ratio

    # Mortality uses the initial biomass ratio because bulk mortality cannot
    # consume biomass produced during this call.
    for ratio, rate_name, biomass_name in zip(
        initial_biomass_ratios,
        (
            "acetate_methanogen_mortality_c",
            "h2_methanogen_mortality_c",
            "aerobic_methanotroph_mortality_c",
            "anaerobic_methanotroph_mortality_c",
        ),
        (
            "acetate_methanogen_c",
            "h2_methanogen_c",
            "aerobic_methanotroph_c",
            "anaerobic_methanotroph_c",
        ),
    ):
        transfer = CATOMW * dt * getattr(rates, rate_name) * ratio
        setattr(t, biomass_name, getattr(t, biomass_name) - transfer)
        dom_delta += transfer
    return dom_delta, t


def advance_carbon_tracer_layer(
    bulk_dom: float,
    tracer_dom: float,
    saturated_fraction: float,
    bulk_unsaturated: State,
    bulk_saturated: State,
    tracer_unsaturated: State,
    tracer_saturated: State,
    transaction: ReactionTransaction,
    parameters: Mapping[str, float],
    dt: float,
) -> tuple[float, State, State, float]:
    fraction = min(1.0, max(0.0, saturated_fraction))
    du, tu = _advance_partition(
        bulk_dom, bulk_unsaturated, tracer_dom, tracer_unsaturated,
        transaction.unsaturated_rates, parameters, dt,
    )
    ds, ts = _advance_partition(
        bulk_dom, bulk_saturated, tracer_dom, tracer_saturated,
        transaction.saturated_rates, parameters, dt,
    )
    new_dom = tracer_dom + (1.0 - fraction) * du + fraction * ds
    initial = tracer_dom + (1.0 - fraction) * additional_carbon_density(
        tracer_unsaturated
    ) + fraction * additional_carbon_density(tracer_saturated)
    final = new_dom + (1.0 - fraction) * additional_carbon_density(tu) + fraction * additional_carbon_density(ts)
    return new_dom, tu, ts, final - initial
