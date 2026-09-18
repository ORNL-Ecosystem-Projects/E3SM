#!/usr/bin/env python3
"""Structural integration tests for the Phase 3 Step 4 ELM adapter."""

from __future__ import annotations

import unittest
from pathlib import Path


PHASE3_DIR = Path(__file__).resolve().parent
ELM_DIR = PHASE3_DIR.parents[2]
ADAPTER = ELM_DIR / "src" / "biogeochem" / "MicrobeMethaneMod.F90"
BALANCE = ELM_DIR / "src" / "biogeochem" / "EcosystemBalanceCheckMod.F90"
BUDGET = ELM_DIR / "src" / "biogeochem" / "CNPBudgetMod.F90"
DRIVER = ELM_DIR / "src" / "main" / "elm_driver.F90"


class Phase3AdapterTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.adapter = ADAPTER.read_text(encoding="utf-8")
        cls.balance = BALANCE.read_text(encoding="utf-8")
        cls.budget = BUDGET.read_text(encoding="utf-8")
        cls.driver = DRIVER.read_text(encoding="utf-8")

    def test_adapter_composes_all_transaction_kernels(self) -> None:
        self.assertIn("procedure, public  :: Advance", self.adapter)
        for routine in (
            "advanceMicrobeMethaneReactionLayer",
            "advanceMicrobeMethaneAcetateTransport",
            "advanceMicrobeMethaneGasTransport",
        ):
            self.assertIn(f"call {routine}", self.adapter)

        validation = self.adapter.index("if (.not. column_valid)")
        commit = self.adapter.index("call commitColumnState")
        self.assertLess(validation, commit)

    def test_authoritative_elm_pools_are_gathered_and_committed(self) -> None:
        for statement in (
            "dom_c(j) = col_cs%decomp_cpools_vr(c,j,i_dom)",
            "dom_n(j) = col_ns%decomp_npools_vr(c,j,i_dom)",
            "dom_p(j) = col_ps%decomp_ppools_vr(c,j,i_dom)",
            "mineral_n(j) = col_ns%smin_nh4_vr(c,j)",
            "mineral_p(j) = col_ps%solutionp_vr(c,j)",
            "col_cs%decomp_cpools_vr(c,j,i_dom) = dom_c(j)",
            "col_ns%decomp_npools_vr(c,j,i_dom) = dom_n(j)",
            "col_ps%decomp_ppools_vr(c,j,i_dom) = dom_p(j)",
            "col_ns%smin_nh4_vr(c,j) = mineral_n(j)",
            "col_ns%sminn_vr(c,j) = mineral_n(j) + col_ns%smin_no3_vr(c,j)",
            "col_ps%solutionp_vr(c,j) = mineral_p(j)",
        ):
            self.assertIn(statement, self.adapter)

    def test_transport_uses_elm_hydrology_and_atmospheric_boundary(self) -> None:
        for source in (
            "soilhydrology_vars%fsat_col(c)",
            "col_ws%frac_h2osfc(c)",
            "col_ws%h2osoi_liq(column,layer)",
            "col_ws%h2osoi_ice(column,layer)",
            "soilstate_vars%rootfr_col(c,j)",
            "atm2lnd_vars%forc_pbot_downscaled_col(c)",
            "ground_conductance_patch(bounds%begp:bounds%endp)",
        ):
            self.assertIn(source, self.adapter)
        self.assertIn("division by a", self.adapter)
        self.assertIn("converts it to the s-1 rate", self.adapter)

    def test_storage_flux_diagnostics_and_balance_interfaces_are_explicit(self) -> None:
        for name, units in (
            ("MM_ADDITIONAL_C", "gC/m^2"),
            ("MM_SURFACE_C_FLUX", "gC/m^2/s"),
            ("MM_SURFACE_CH4_FLUX", "kgC/m^2/s"),
            ("MM_SURFACE_CO2_FLUX", "gC/m^2/s"),
        ):
            self.assertIn(f"fname='{name}', units='{units}'", self.adapter)

        self.assertIn("call this%UpdateAdditionalCarbon(bounds)", self.adapter)
        self.assertIn("additional_carbon_col", self.balance)
        self.assertIn("surface_carbon_flux_col", self.balance)
        self.assertIn("surface_carbon_flux_col", self.budget)

    def test_legacy_dispatch_is_unchanged_until_step_5(self) -> None:
        self.assertIn("call CH4", self.driver)
        self.assertNotIn("microbe_methane_vars%Advance", self.driver)


if __name__ == "__main__":
    unittest.main()
