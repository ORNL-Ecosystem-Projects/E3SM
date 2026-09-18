#!/usr/bin/env python3
"""Conservation tests for Phase 3 Step 4 state transactions."""

from __future__ import annotations

import unittest
from pathlib import Path

import phase3
from reaction_oracle import CATOMW, Environment, State, compute
from state_update_oracle import (
    additional_carbon_density,
    advance_reaction_layer,
    ch4_surface_flux_kgc,
    co2_correction,
    column_additional_carbon,
    surface_carbon_flux,
)


PHASE3_DIR = Path(__file__).resolve().parent
ELM_DIR = PHASE3_DIR.parents[2]
MANIFEST = PHASE3_DIR / "phase3_reference_parameters.json"
STATE_UPDATE = ELM_DIR / "src" / "biogeochem" / "MicrobeMethaneStateUpdateMod.F90"
REACTIONS = ELM_DIR / "src" / "biogeochem" / "MicrobeMethaneReactionMod.F90"
DRIVER = ELM_DIR / "src" / "main" / "elm_driver.F90"


class Phase3StateUpdateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        document = phase3.load_manifest(MANIFEST)
        cls.parameters = {
            name: float(details[1])
            for name, details in phase3.parameter_map(document).items()
        }
        cls.source = STATE_UPDATE.read_text(encoding="utf-8")

    def test_functional_biomass_floor_is_reserved_before_commit(self) -> None:
        floor = self.parameters["microbe_methane_mfg_biomass_min"]
        state = State(
            acetate_methanogen_c=floor,
            h2_methanogen_c=floor,
            aerobic_methanotroph_c=floor,
            anaerobic_methanotroph_c=floor,
        )
        rates, tendencies = compute(
            state,
            Environment(soil_temperature=286.65, soil_ph=7.0),
            self.parameters,
            86400.0,
        )
        for name in (
            "acetate_methanogen_mortality_c",
            "h2_methanogen_mortality_c",
            "aerobic_methanotroph_mortality_c",
            "anaerobic_methanotroph_mortality_c",
        ):
            self.assertEqual(getattr(rates, name), 0.0)
        self.assertEqual(tendencies.dom_c, 0.0)
        self.assertIn("state%acetate_methanogen_c - parameters%mfg_biomass_min", REACTIONS.read_text())

    def test_partitioned_reaction_transaction_closes_c_n_and_p(self) -> None:
        state = State(
            acetate_c=0.8,
            acetate_methanogen_c=0.2,
            h2_methanogen_c=0.2,
            aerobic_methanotroph_c=0.2,
            anaerobic_methanotroph_c=0.2,
            conc_ch4=0.003,
            conc_o2=0.015,
            conc_co2=0.01,
            conc_h2=0.001,
        )
        transaction = advance_reaction_layer(
            dom_c=12.0,
            dom_n=1.2,
            dom_p=0.12,
            mineral_n=0.5,
            mineral_p=0.05,
            saturated_fraction=0.35,
            unsaturated_state=state,
            saturated_state=state,
            unsaturated_environment=Environment(286.65, 7.0, 0.25, 0.75),
            saturated_environment=Environment(286.65, 7.0, 1.0, 0.0),
            parameters=self.parameters,
            cn_dom=10.0,
            cp_dom=100.0,
            dt=1800.0,
        )
        self.assertTrue(transaction.valid)
        self.assertAlmostEqual(transaction.carbon_residual, 0.0, places=12)
        self.assertAlmostEqual(transaction.nitrogen_residual, 0.0, places=14)
        self.assertAlmostEqual(transaction.phosphorus_residual, 0.0, places=14)

    def test_dom_reaction_is_limited_by_associated_n_and_p(self) -> None:
        state = State(
            acetate_methanogen_c=1.0,
            h2_methanogen_c=1.0,
            aerobic_methanotroph_c=1.0,
            anaerobic_methanotroph_c=1.0,
        )
        transaction = advance_reaction_layer(
            dom_c=100.0,
            dom_n=1.0e-8,
            dom_p=1.0e-9,
            mineral_n=0.0,
            mineral_p=0.0,
            saturated_fraction=0.5,
            unsaturated_state=state,
            saturated_state=state,
            unsaturated_environment=Environment(306.65, 7.0, 1.0, 0.0),
            saturated_environment=Environment(306.65, 7.0, 1.0, 0.0),
            parameters=self.parameters,
            cn_dom=10.0,
            cp_dom=100.0,
            dt=86400.0,
        )
        self.assertTrue(transaction.valid)
        self.assertGreaterEqual(transaction.dom_n, -1.0e-14)
        self.assertGreaterEqual(transaction.dom_p, -1.0e-14)

    def test_mortality_never_immobilizes_unavailable_mineral_nutrients(self) -> None:
        state = State(
            acetate_methanogen_c=10.0,
            h2_methanogen_c=10.0,
            aerobic_methanotroph_c=10.0,
            anaerobic_methanotroph_c=10.0,
        )
        transaction = advance_reaction_layer(
            dom_c=1.0,
            dom_n=0.1,
            dom_p=0.01,
            mineral_n=0.0,
            mineral_p=0.0,
            saturated_fraction=0.5,
            unsaturated_state=state,
            saturated_state=state,
            unsaturated_environment=Environment(286.65, 7.0, 0.0, 0.0),
            saturated_environment=Environment(286.65, 7.0, 0.0, 0.0),
            parameters=self.parameters,
            cn_dom=10.0,
            cp_dom=100.0,
            dt=86400.0,
        )
        self.assertTrue(transaction.valid)
        self.assertAlmostEqual(transaction.dom_c, 1.0)
        self.assertAlmostEqual(transaction.mineral_n, 0.0)
        self.assertAlmostEqual(transaction.mineral_p, 0.0)

    def test_budget_conversions_have_explicit_units_and_sign(self) -> None:
        fluxes = [2.0e-7, -1.0e-6, 3.0e-7, 4.0e-8]
        self.assertAlmostEqual(surface_carbon_flux(fluxes), CATOMW * 5.0e-7)
        self.assertAlmostEqual(ch4_surface_flux_kgc(fluxes), CATOMW * 2.0e-10)
        self.assertAlmostEqual(co2_correction(fluxes), CATOMW * 3.0e-7)
        state = State(acetate_c=2.0, conc_ch4=0.1, conc_co2=0.2)
        self.assertAlmostEqual(
            additional_carbon_density(state), 2.0 + CATOMW * 0.3
        )
        self.assertAlmostEqual(
            column_additional_carbon(0.25, [state], [State()], [0.4]),
            0.75 * additional_carbon_density(state) * 0.4,
        )

    def test_state_update_is_pure_and_still_not_dispatched(self) -> None:
        for routine in (
            "advanceMicrobeMethaneReactionLayer",
            "advanceMicrobeMethaneGasTransport",
            "advanceMicrobeMethaneAcetateTransport",
        ):
            self.assertIn(f"pure subroutine {routine}", self.source)
            self.assertNotIn(routine, DRIVER.read_text(encoding="utf-8"))
        self.assertIn("DOM is excluded", self.source)
        self.assertIn("positive upward", self.source)


if __name__ == "__main__":
    unittest.main()
