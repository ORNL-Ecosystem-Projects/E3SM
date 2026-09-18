#!/usr/bin/env python3
"""Unit and closed-box tests for the Phase 3 Step 2 reaction kernel."""

from __future__ import annotations

import re
import unittest
from dataclasses import fields
from pathlib import Path

import phase3
from reaction_oracle import (
    Environment,
    Rates,
    State,
    assemble,
    carbon_residual,
    compute,
    ph_response,
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
            if item.name not in {"effective_soil_ph", "ph_response"}:
                self.assertEqual(getattr(rates, item.name), 0.0)
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
