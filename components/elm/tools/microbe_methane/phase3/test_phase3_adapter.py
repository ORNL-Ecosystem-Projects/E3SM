#!/usr/bin/env python3
"""Structural integration tests for the Phase 3 ELM adapter and dispatch."""

from __future__ import annotations

import unittest
from pathlib import Path


PHASE3_DIR = Path(__file__).resolve().parent
ELM_DIR = PHASE3_DIR.parents[2]
ADAPTER = ELM_DIR / "src" / "biogeochem" / "MicrobeMethaneMod.F90"
BALANCE = ELM_DIR / "src" / "biogeochem" / "EcosystemBalanceCheckMod.F90"
BUDGET = ELM_DIR / "src" / "biogeochem" / "CNPBudgetMod.F90"
NITRIFICATION = ELM_DIR / "src" / "biogeochem" / "NitrifDenitrifMod.F90"
LEGACY_CH4 = ELM_DIR / "src" / "biogeochem" / "CH4Mod.F90"
DRIVER = ELM_DIR / "src" / "main" / "elm_driver.F90"
LND2ATM = ELM_DIR / "src" / "main" / "lnd2atmMod.F90"
SOIL_TRANSPORT = ELM_DIR / "src" / "biogeochem" / "SoilLittVertTranspMod.F90"
VARCTL = ELM_DIR / "src" / "main" / "elm_varctl.F90"
NAMELIST_DEFINITION = ELM_DIR / "bld" / "namelist_files" / "namelist_definition.xml"
CASCADE = ELM_DIR / "src" / "biogeochem" / "DecompCascadeCNMod.F90"


class Phase3AdapterTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.adapter = ADAPTER.read_text(encoding="utf-8")
        cls.balance = BALANCE.read_text(encoding="utf-8")
        cls.budget = BUDGET.read_text(encoding="utf-8")
        cls.nitrification = NITRIFICATION.read_text(encoding="utf-8")
        cls.legacy_ch4 = LEGACY_CH4.read_text(encoding="utf-8")
        cls.driver = DRIVER.read_text(encoding="utf-8")
        cls.lnd2atm = LND2ATM.read_text(encoding="utf-8")
        cls.soil_transport = SOIL_TRANSPORT.read_text(encoding="utf-8")
        cls.varctl = VARCTL.read_text(encoding="utf-8")
        cls.namelist_definition = NAMELIST_DEFINITION.read_text(encoding="utf-8")
        cls.cascade = CASCADE.read_text(encoding="utf-8")

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

    def test_clm_saturated_biomass_reset_is_not_ported(self) -> None:
        self.assertNotIn("1.e-5_r8 * catomw", self.adapter)
        self.assertNotIn("initial_biomass_carbon_flux_col", self.adapter)
        self.assertNotIn("legacy_biomass_reset", self.adapter)

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
            "soilstate_vars%rootfr_patch(bounds%begp:bounds%endp,1:nlevdecomp)",
            "validRootFraction(root_fraction_col(c,j))",
            "atm2lnd_vars%forc_pbot_downscaled_col(c)",
            "ground_conductance_patch(bounds%begp:bounds%endp)",
        ):
            self.assertIn(source, self.adapter)
        self.assertIn("division by a", self.adapter)
        self.assertIn("converts it to the s-1 rate", self.adapter)
        self.assertIn("call p2c(bounds, nlevdecomp", self.adapter)
        self.assertIn("value_in <= 1._r8", self.adapter)

    def test_physical_aqueous_transport_has_one_authoritative_path(self) -> None:
        self.assertIn("1.e-3_r8 * col_wf%qflx_adv", self.adapter)
        self.assertIn("water_flux(0) = max(0._r8, water_flux(0))", self.adapter)
        self.assertGreaterEqual(
            self.adapter.count("call advanceMicrobeAqueousTracerTransport"), 5
        )
        self.assertIn(
            ".not. (use_microbe_aqueous_transport .and. is_dissolved(s))",
            self.soil_transport,
        )

    def test_nonbog_lateral_gas_edge_has_a_default_on_diagnostic_gate(self) -> None:
        self.assertIn("use_microbe_nonbog_lateral_gas_transport", self.adapter)
        self.assertIn(
            "if (.not. use_microbe_nonbog_lateral_gas_transport) then",
            self.adapter,
        )
        self.assertIn(
            ".not. top_pp%is_bog(t) .or. .not. top_pp%is_bog(target_t)",
            self.adapter,
        )

    def test_aqueous_exports_are_in_balance_and_budget_interfaces(self) -> None:
        for field in (
            "aqueous_carbon_export_col",
            "aqueous_nitrogen_export_col",
            "aqueous_phosphorus_export_col",
        ):
            self.assertIn(field, self.driver)
            self.assertIn(field, self.balance)
        self.assertIn("aqueous_carbon_export_col", self.budget)
        # The reconstructed branch has no separate DOC-runoff budget slot;
        # book physical DOM export in its existing SOM-C leaching category.
        self.assertIn("nf = f_som_c_leached", self.budget)

    def test_native_elm_uses_explicit_soil_ph_fallback(self) -> None:
        self.assertIn(
            "unsaturated_environment%soil_ph = MicrobeMethaneParamsInst%soil_ph_fallback",
            self.adapter,
        )
        self.assertIn(
            "saturated_environment%soil_ph = MicrobeMethaneParamsInst%soil_ph_fallback",
            self.adapter,
        )
        self.assertNotIn(
            "unsaturated_environment%soil_ph = chemstate_vars%soil_pH(c,j)",
            self.adapter,
        )

    def test_storage_flux_diagnostics_and_balance_interfaces_are_explicit(self) -> None:
        for name, units in (
            ("MM_ADDITIONAL_C", "gC/m^2"),
            ("MM_SURFACE_C_FLUX", "gC/m^2/s"),
            ("MM_SURFACE_CH4_FLUX", "kgC/m^2/s"),
            ("MM_SURFACE_CO2_FLUX", "gC/m^2/s"),
            ("MM_CH4_PROD", "gC/m^2/s"),
            ("MM_CH4_OXID", "gC/m^2/s"),
            ("MM_CH4_PROD_UNSAT", "gC/m^2/s"),
            ("MM_CH4_PROD_SAT", "gC/m^2/s"),
            ("MM_CH4_OXID_UNSAT", "gC/m^2/s"),
            ("MM_CH4_OXID_SAT", "gC/m^2/s"),
            ("MM_CH4_OXID_AER", "gC/m^2/s"),
            ("MM_CH4_OXID_AER_UNSAT", "gC/m^2/s"),
            ("MM_CH4_OXID_AER_SAT", "gC/m^2/s"),
            ("MM_CH4_OXID_AER_PREO2", "gC/m^2/s"),
            ("MM_CH4_OXID_AER_PREO2_UNSAT", "gC/m^2/s"),
            ("MM_CH4_OXID_AER_PREO2_SAT", "gC/m^2/s"),
            ("MM_CH4_OXID_AOM", "gC/m^2/s"),
            ("MM_CH4_OXID_AOM_UNSAT", "gC/m^2/s"),
            ("MM_CH4_OXID_AOM_SAT", "gC/m^2/s"),
            ("MM_O2_STRESS_UNSAT", "1"),
            ("MM_O2_STRESS_SAT", "1"),
            ("MM_SAT_FRACTION", "1"),
        ):
            self.assertIn(f"fname='{name}', units='{units}'", self.adapter)

        self.assertIn("call this%UpdateAdditionalCarbon(bounds)", self.adapter)
        self.assertIn("additional_carbon_col", self.balance)
        self.assertIn("surface_carbon_flux_col", self.balance)
        self.assertIn("lateral_carbon_flux_col", self.balance)
        self.assertIn("surface_carbon_flux_col", self.budget)

        for rate in (
            "acetoclastic_methanogenesis_c",
            "hydrogenotrophic_methanogenesis_c",
            "aerobic_methane_oxidation_c",
            "anaerobic_methane_oxidation_c",
        ):
            self.assertIn(rate, self.adapter)

        production_start = self.adapter.index("this%ch4_production_unsat_col(c) =")
        oxidation_start = self.adapter.index("this%ch4_oxidation_unsat_col(c) =", production_start)
        production_diagnostic = self.adapter[production_start:oxidation_start]
        self.assertEqual(
            production_diagnostic.count("acetoclastic_methanogenesis_ch4_yield"),
            2,
        )
        self.assertIn(
            "this%ch4_production_col(c) = this%ch4_production_unsat_col(c)",
            self.adapter,
        )
        self.assertIn(
            "this%ch4_oxidation_col(c) = this%ch4_oxidation_unsat_col(c)",
            self.adapter,
        )
        self.assertIn(
            "this%ch4_aerobic_oxidation_col(c) = this%ch4_aerobic_oxidation_unsat_col(c)",
            self.adapter,
        )
        self.assertIn(
            "this%ch4_aerobic_oxidation_pre_o2_col(c) =",
            self.adapter,
        )
        self.assertIn(
            "this%ch4_anaerobic_oxidation_col(c) = this%ch4_anaerobic_oxidation_unsat_col(c)",
            self.adapter,
        )
        self.assertIn(
            "this%ch4_oxidation_unsat_col(c) = this%ch4_aerobic_oxidation_unsat_col(c)",
            self.adapter,
        )
        self.assertEqual(
            production_diagnostic.count("1._r8 - MicrobeMethaneParamsInst%acetate_methanogen_yield"),
            2,
        )

    def test_step_5_dispatch_selects_exactly_one_backend(self) -> None:
        revised = self.driver.index("microbe_methane_vars%Advance")
        legacy = self.driver.index("call CH4", revised)
        dispatch = self.driver.rfind("if (use_lch4", 0, revised)
        fallback = self.driver.index("else if (use_lch4", revised, legacy)

        self.assertLess(dispatch, revised)
        self.assertLess(revised, fallback)
        self.assertLess(fallback, legacy)

    def test_phase2_legacy_methane_bridge_is_an_explicit_default_off_option(self) -> None:
        self.assertIn(
            "logical, public :: use_legacy_ch4_with_microbe = .false.",
            self.varctl,
        )
        self.assertIn('id="use_legacy_ch4_with_microbe"', self.namelist_definition)
        self.assertIn('value=".false."', self.namelist_definition)

        revised = self.driver.index("microbe_methane_vars%Advance")
        dispatch = self.driver.rfind("if (use_lch4", 0, revised)
        self.assertIn(".not. use_legacy_ch4_with_microbe", self.driver[dispatch:revised])
        self.assertIn("if (use_microbe_methane) then", self.cascade)
        self.assertIn(
            ".not. use_microbe_methane .or. use_legacy_ch4_with_microbe",
            self.legacy_ch4,
        )
        self.assertGreaterEqual(
            self.lnd2atm.count(
                "use_microbe_methane .and. .not. use_legacy_ch4_with_microbe"
            ),
            2,
        )

    def test_revised_accounting_is_option_gated(self) -> None:
        for field in (
            "additional_carbon_col=microbe_methane_vars%additional_carbon_col",
            "surface_carbon_flux_col=microbe_methane_vars%surface_carbon_flux_col",
            "lateral_carbon_flux_col=microbe_methane_vars%lateral_carbon_flux_col",
        ):
            self.assertIn(field, self.driver)

        self.assertIn(
            "if (use_microbe_methane .and. .not. use_legacy_ch4_with_microbe) then",
            self.driver,
        )
        self.assertIn("call CH4", self.driver)

    def test_revised_fluxes_feed_atmosphere_without_changing_offline_policy(self) -> None:
        self.assertIn("microbe_methane_vars%surface_ch4_flux_col", self.lnd2atm)
        self.assertIn("microbe_methane_vars%surface_co2_flux_col", self.lnd2atm)
        self.assertIn("if (.not. ch4offline) then", self.lnd2atm)

    def test_revised_adapter_publishes_legacy_oxygen_interface(self) -> None:
        self.assertIn("call syncLegacyOxygenBridge", self.adapter)
        for field in (
            "ch4_vars%finundated_col",
            "ch4_vars%conc_o2_unsat_col",
            "ch4_vars%conc_o2_sat_col",
            "ch4_vars%o2stress_unsat_col",
            "ch4_vars%o2stress_sat_col",
            "ch4_vars%o2_decomp_depth_unsat_col",
            "ch4_vars%o2_decomp_depth_sat_col",
        ):
            self.assertIn(field, self.adapter)

        self.assertIn(
            ".not. use_microbe_methane .or. use_legacy_ch4_with_microbe",
            self.legacy_ch4,
        )
        self.assertIn(
            "max(0._r8, unsaturated(layer)%conc_o2)",
            self.adapter,
        )
        self.assertIn(
            "max(0._r8, saturated(layer)%conc_o2)",
            self.adapter,
        )

    def test_all_aerobic_processes_share_revised_oxygen(self) -> None:
        self.assertNotIn("bulk_o2_available_rate", self.nitrification)
        self.assertNotIn("call consumeELMAerobicOxygen", self.adapter)
        self.assertIn("elm_aerobic_o2_demand", self.adapter)
        self.assertIn("elmAerobicOxygenDemand", self.adapter)
        self.assertIn("col_cf%phr_vr(column,layer)", self.adapter)
        self.assertNotIn("col_cf%rr_vr(column,layer)", self.adapter)
        self.assertIn("root_respiration_col_vr", self.adapter)
        self.assertIn("veg_cf%rr(p) * soilstate_vars%rootfr_patch(p,j)", self.adapter)
        self.assertIn("filter(nc)%num_soilp, filter(nc)%soilp", self.driver)
        self.assertIn("col_nf%pot_f_nit_vr(column,layer)", self.adapter)
        self.assertIn("2._r8 / 14._r8", self.adapter)
        self.assertIn("%unsaturated_rates%oxygen_stress", self.adapter)
        self.assertIn("%saturated_rates%oxygen_stress", self.adapter)


if __name__ == "__main__":
    unittest.main()
