#!/usr/bin/env python3
"""Unit and closed-box tests for the Phase 3 Step 2 reaction kernel."""

from __future__ import annotations

import math
import re
import unittest
from dataclasses import fields
from pathlib import Path

import phase3
from reaction_oracle import (
    CATOMW,
    Environment,
    Rates,
    State,
    assemble,
    carbon_residual,
    compute,
    limit_rates,
    ph_response,
    potential_rates,
    q10_response,
)


PHASE3_DIR = Path(__file__).resolve().parent
ELM_DIR = PHASE3_DIR.parents[2]
MANIFEST = PHASE3_DIR / "phase3_reference_parameters.json"
KERNEL = ELM_DIR / "src" / "biogeochem" / "MicrobeMethaneReactionMod.F90"
DRIVER = ELM_DIR / "src" / "main" / "elm_driver.F90"
DESIGN = ELM_DIR / "docs" / "dev-guide" / "spruce-microbe-methane-design.md"


class Phase3ReactionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        document = phase3.load_manifest(MANIFEST)
        cls.parameters = {
            name: float(details[1])
            for name, details in phase3.parameter_map(document).items()
        }
        cls.source = KERNEL.read_text(encoding="utf-8")

    def test_kernel_is_side_effect_free_and_not_dispatched(self) -> None:
        self.assertIn("pure subroutine computeMicrobeMethaneReactionTendencies", self.source)
        self.assertIn("intent(in) :: state", self.source)
        self.assertNotIn("use MicrobeMethaneMod", self.source)
        self.assertNotRegex(self.source, r"state%[a-z0-9_]+\s*=")
        self.assertNotIn("computeMicrobeMethaneReactionTendencies", DRIVER.read_text(encoding="utf-8"))
        for routine in (
            "computeMicrobeMethanePotentialRates",
            "limitMicrobeMethaneReactionRates",
            "assembleMicrobeMethaneReactionTendencies",
        ):
            self.assertIn(f"public :: {routine}", self.source)

    def test_temperature_and_ph_response_functions(self) -> None:
        self.assertAlmostEqual(q10_response(2.0, 296.65, 286.65), 2.0)
        self.assertAlmostEqual(q10_response(2.0, 276.65, 286.65), 0.5)
        self.assertAlmostEqual(ph_response(7.0, 4.0, 7.0, 10.0), 1.0)
        self.assertEqual(ph_response(4.0, 4.0, 7.0, 10.0), 0.0)
        self.assertEqual(ph_response(10.0, 4.0, 7.0, 10.0), 0.0)
        self.assertGreater(ph_response(6.0, 4.0, 7.0, 10.0), 0.0)

    def test_zero_state_has_zero_tendencies(self) -> None:
        rates, tendencies = compute(State(), Environment(), self.parameters, 1800.0)
        for item in fields(Rates):
            if item.name not in {"effective_soil_ph", "ph_response", "oxygen_stress"}:
                self.assertEqual(getattr(rates, item.name), 0.0)
        self.assertEqual(rates.oxygen_stress, 1.0)
        for item in fields(State):
            self.assertEqual(getattr(tendencies, item.name), 0.0)

    def test_every_reaction_and_mortality_path_conserves_carbon(self) -> None:
        rate_names = [
            item.name
            for item in fields(Rates)
            if item.name not in {"effective_soil_ph", "ph_response"}
        ]
        for name in rate_names:
            with self.subTest(rate=name):
                rates = Rates()
                setattr(rates, name, 1.25e-7)
                tendencies = assemble(self.parameters, rates)
                self.assertAlmostEqual(carbon_residual(tendencies), 0.0, places=20)

    def test_combined_closed_box_conserves_carbon(self) -> None:
        rates = Rates(
            dom_to_acetate_c=1.1e-8,
            acetogenesis_c=2.2e-8,
            acetoclastic_methanogenesis_c=3.3e-8,
            hydrogenotrophic_methanogenesis_c=4.4e-8,
            aerobic_acetate_oxidation_c=5.5e-8,
            aerobic_methane_oxidation_c=6.6e-8,
            anaerobic_methane_oxidation_c=7.7e-8,
            acetate_methanogen_mortality_c=8.8e-9,
            h2_methanogen_mortality_c=9.9e-9,
            aerobic_methanotroph_mortality_c=1.2e-8,
            anaerobic_methanotroph_mortality_c=1.3e-8,
        )
        self.assertAlmostEqual(carbon_residual(assemble(self.parameters, rates)), 0.0, places=20)

    def test_competing_substrate_limiter_prevents_negative_state(self) -> None:
        state = State(
            dom_c=1.0e-8,
            acetate_c=1.0e-8,
            acetate_methanogen_c=2.0,
            h2_methanogen_c=2.0,
            aerobic_methanotroph_c=2.0,
            anaerobic_methanotroph_c=2.0,
            conc_ch4=1.0e-10,
            conc_o2=1.0e-10,
            conc_co2=1.0e-10,
            conc_h2=1.0e-10,
        )
        environment = Environment(
            soil_temperature=306.65,
            soil_ph=7.0,
            dom_fermentation_scalar=1.0,
            aerobic_acetate_oxidation_scalar=1.0,
        )
        dt = 86400.0
        _rates, tendencies = compute(state, environment, self.parameters, dt)
        for item in fields(State):
            updated = getattr(state, item.name) + dt * getattr(tendencies, item.name)
            self.assertGreaterEqual(updated, -1.0e-14, item.name)
        self.assertAlmostEqual(carbon_residual(tendencies), 0.0, places=20)

    def test_elm_and_microbial_processes_share_finite_oxygen(self) -> None:
        dt = 1.0
        state = State(acetate_c=1.0, conc_ch4=1.0, conc_o2=3.5e-7)
        environment = Environment(elm_aerobic_o2_demand=3.0e-7)
        rates = Rates(aerobic_methane_oxidation_c=2.0e-7)
        limited = limit_rates(state, environment, self.parameters, dt, rates)
        tendencies = assemble(self.parameters, limited)

        # Total potential demand is 3e-7 from ELM plus 4e-7 from methane
        # oxidation. A 3.5e-7 supply therefore applies a common 0.5 stress.
        self.assertAlmostEqual(limited.oxygen_stress, 0.5)
        self.assertAlmostEqual(limited.elm_aerobic_o2_consumption, 1.5e-7)
        self.assertAlmostEqual(limited.aerobic_methane_oxidation_pre_o2_c, 2.0e-7)
        self.assertAlmostEqual(limited.aerobic_methane_oxidation_c, 1.0e-7)
        self.assertAlmostEqual(tendencies.conc_o2, -state.conc_o2 / dt)

    def test_one_layer_legacy_executable_rate_parity(self) -> None:
        """Unit-converted ELM parameters preserve comparable legacy rates."""
        molar_state = {
            "dom": 1.5,
            "acetate": 0.2,
            "acetate_methanogen": 0.01,
            "h2_methanogen": 0.08,
            "aerobic_methanotroph": 0.015,
            "anaerobic_methanotroph": 0.012,
            "ch4": 0.3,
            "o2": 0.5,
            "co2": 0.08,
            "h2": 0.004,
        }
        state = State(
            dom_c=molar_state["dom"] * CATOMW,
            acetate_c=molar_state["acetate"] * CATOMW,
            acetate_methanogen_c=molar_state["acetate_methanogen"] * CATOMW,
            h2_methanogen_c=molar_state["h2_methanogen"] * CATOMW,
            aerobic_methanotroph_c=molar_state["aerobic_methanotroph"] * CATOMW,
            anaerobic_methanotroph_c=molar_state["anaerobic_methanotroph"] * CATOMW,
            conc_ch4=molar_state["ch4"],
            conc_o2=molar_state["o2"],
            conc_co2=molar_state["co2"],
            conc_h2=molar_state["h2"],
        )
        environment = Environment(
            soil_temperature=286.65,
            soil_ph=6.5,
            dom_fermentation_scalar=0.65,
            aerobic_acetate_oxidation_scalar=1.0,
        )
        elm_rates = potential_rates(state, environment, self.parameters)

        def monod(substrate: float, half_saturation: float) -> float:
            return substrate / (substrate + half_saturation)

        effective_ph = -math.log10(
            10.0 ** (-environment.soil_ph)
            + 4.2e-9 * molar_state["acetate"]
        )
        ph_scalar = ph_response(effective_ph, 4.0, 7.0, 10.0)
        acetate_feedback = 0.1 / (molar_state["acetate"] + 0.1)
        co2_inhibition = 1.0 - min(1.0, molar_state["co2"] / 9.2)
        o2_inhibition = 1.0 - min(1.0, molar_state["o2"] / 4.6)

        # These are the active CLM-SPRUCE expressions evaluated after its
        # gC/12 conversion, so concentrations are mol m-3 and rates are per s.
        legacy_rates = {
            "dom_to_acetate_c": (
                (2.0 / 3.0) * 2.4e-6
                * monod(molar_state["dom"], 16.0)
                * ph_scalar * acetate_feedback
                * environment.dom_fermentation_scalar
            ),
            "acetogenesis_c": (
                5.0e-8
                * monod(molar_state["h2"], 0.00165)
                * monod(molar_state["h2_methanogen"], 0.0825)
                * ph_scalar
            ),
            "hydrogenotrophic_methanogenesis_c": (
                0.01 / 0.015 * molar_state["h2_methanogen"]
                * monod(molar_state["h2"], 7.75e-5)
                * monod(molar_state["co2"], 3.1e-4)
                * ph_scalar * co2_inhibition
            ),
            "acetoclastic_methanogenesis_c": (
                0.008 / 0.2 * molar_state["acetate_methanogen"]
                * monod(molar_state["acetate"], 0.05)
                * ph_scalar
            ),
            "aerobic_acetate_oxidation_c": (
                0.05 * min(molar_state["acetate"], molar_state["o2"])
                * monod(molar_state["o2"], 0.004)
            ),
            "aerobic_methane_oxidation_c": (
                0.008 / 0.4 * molar_state["aerobic_methanotroph"]
                * monod(molar_state["ch4"], 1.0)
                * monod(molar_state["o2"], 4.0)
                * ph_scalar
            ),
            "acetate_methanogen_mortality_c": (
                0.002 * molar_state["acetate_methanogen"] * ph_scalar
            ),
            "h2_methanogen_mortality_c": (
                0.001 * molar_state["h2_methanogen"] * ph_scalar
            ),
            "aerobic_methanotroph_mortality_c": (
                0.002 * molar_state["aerobic_methanotroph"] * ph_scalar
            ),
            "anaerobic_methanotroph_mortality_c": (
                0.002 * molar_state["anaerobic_methanotroph"] * ph_scalar
            ),
        }
        for rate_name, legacy_rate in legacy_rates.items():
            with self.subTest(rate=rate_name):
                self.assertTrue(
                    math.isclose(
                        getattr(elm_rates, rate_name),
                        legacy_rate,
                        rel_tol=2.0e-14,
                        abs_tol=0.0,
                    ),
                    f"ELM={getattr(elm_rates, rate_name)!r}, CLM={legacy_rate!r}",
                )

        # Legacy AOM subtracts 13.5 (C) directly from a Kelvin temperature.
        # ELM intentionally uses 286.65 K, so this is a classified source
        # repair rather than a parameter-conversion parity failure.
        legacy_aom = (
            0.004 / 0.15 * molar_state["anaerobic_methanotroph"]
            * monod(molar_state["ch4"], 1.5)
            * 1.2 ** ((environment.soil_temperature - 13.5) / 10.0)
            * ph_scalar * o2_inhibition
        )
        self.assertFalse(math.isclose(elm_rates.anaerobic_methane_oxidation_c, legacy_aom))
        self.assertAlmostEqual(
            legacy_aom / elm_rates.anaerobic_methane_oxidation_c,
            1.2 ** ((286.65 - 13.5) / 10.0),
        )

    def test_named_parameters_and_conservation_repairs_are_explicit(self) -> None:
        used = set(re.findall(r"parameters%([a-z0-9_]+)", self.source))
        required = {
            "acetate_prod_max",
            "acetogenesis_max",
            "h2_methanogen_yield",
            "acetate_methanogen_yield",
            "aerobic_methanotroph_yield",
            "anaerobic_methanotroph_yield",
            "aerobic_acetate_oxidation_rate",
            "aerobic_oxidation_o2_ch4_ratio",
        }
        self.assertTrue(required <= used)
        design = DESIGN.read_text(encoding="utf-8")
        self.assertIn("mortality carbon returns to DOM", design)
        self.assertIn("additional factor of four", design)
        self.assertRegex(
            design,
            r"dissolved CO2\s+rather than methanogen biomass",
        )
        self.assertIn("does not reproduce those carbon leaks", design)


if __name__ == "__main__":
    unittest.main()
