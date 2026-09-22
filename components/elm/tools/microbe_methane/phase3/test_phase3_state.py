#!/usr/bin/env python3
"""Focused structural tests for Phase 3 Step 1."""

from __future__ import annotations

import re
import sys
import tempfile
import unittest
from pathlib import Path


PHASE3_DIR = Path(__file__).resolve().parent
ELM_DIR = PHASE3_DIR.parents[2]
sys.path.insert(0, str(PHASE3_DIR))

import phase3  # noqa: E402


MANIFEST = PHASE3_DIR / "phase3_reference_parameters.json"
PARAMS_FORTRAN = ELM_DIR / "src" / "biogeochem" / "MicrobeMethaneParamsMod.F90"
STATE_FORTRAN = ELM_DIR / "src" / "biogeochem" / "MicrobeMethaneMod.F90"
READ_PARAMS_FORTRAN = ELM_DIR / "src" / "main" / "readParamsMod.F90"
RESTART_FORTRAN = ELM_DIR / "src" / "main" / "restFileMod.F90"
DRIVER_FORTRAN = ELM_DIR / "src" / "main" / "elm_driver.F90"
DESIGN_DOCUMENT = ELM_DIR / "docs" / "dev-guide" / "spruce-microbe-methane-design.md"

STATE_HISTORY_NAMES = {
    "MM_ACETATE_C_UNSAT", "MM_ACETATE_C_SAT",
    "MM_ACET_METH_C_UNSAT", "MM_ACET_METH_C_SAT",
    "MM_H2_METH_C_UNSAT", "MM_H2_METH_C_SAT",
    "MM_AER_METHANOTROPH_C_UNSAT", "MM_AER_METHANOTROPH_C_SAT",
    "MM_ANAER_METHANOTROPH_C_UNSAT", "MM_ANAER_METHANOTROPH_C_SAT",
    "MM_CONC_CH4_UNSAT", "MM_CONC_CH4_SAT",
    "MM_CONC_O2_UNSAT", "MM_CONC_O2_SAT",
    "MM_CONC_CO2_UNSAT", "MM_CONC_CO2_SAT",
    "MM_CONC_H2_UNSAT", "MM_CONC_H2_SAT",
}

CLM_MICROBE_RUNTIME_VALUES = {
    "microbe_methane_k_acetate": 16.0,
    "microbe_methane_acetate_prod_max": 2.4e-6,
    "microbe_methane_k_acetate_prod_o2": 0.004,
    "microbe_methane_acetogenesis_max": 5.0e-8,
    "microbe_methane_k_acetogenesis_h2": 0.00165,
    "microbe_methane_k_acetogenesis_co2": 0.0825,
    "microbe_methane_acetogenesis_q10": 2.0,
    "microbe_methane_h2_methanogen_growth_rate": 0.01,
    "microbe_methane_h2_methanogen_death_rate": 0.001,
    "microbe_methane_h2_methanogen_yield": 0.015,
    "microbe_methane_acetate_methanogen_growth_rate": 0.008,
    "microbe_methane_acetate_methanogen_death_rate": 0.002,
    "microbe_methane_acetate_methanogen_yield": 0.2,
    "microbe_methane_aerobic_methanotroph_growth_rate": 0.008,
    "microbe_methane_aerobic_methanotroph_death_rate": 0.002,
    "microbe_methane_aerobic_methanotroph_yield": 0.4,
    "microbe_methane_anaerobic_methanotroph_growth_rate": 0.004,
    "microbe_methane_anaerobic_methanotroph_death_rate": 0.002,
    "microbe_methane_anaerobic_methanotroph_yield": 0.15,
    "microbe_methane_dom_to_acetate_q10": 3.0,
    "microbe_methane_k_h2_methanogenesis_h2": 7.75e-5,
    "microbe_methane_k_h2_methanogenesis_co2": 3.1e-4,
    "microbe_methane_h2_methanogenesis_q10": 2.0,
    "microbe_methane_k_acetoclastic_methanogenesis_acetate": 0.05,
    "microbe_methane_acetoclastic_methanogenesis_q10": 2.0,
    "microbe_methane_acetoclastic_methanogenesis_ch4_yield": 0.5,
    "microbe_methane_k_aerobic_oxidation_ch4": 1.0,
    "microbe_methane_k_aerobic_oxidation_o2": 4.0,
    "microbe_methane_aerobic_oxidation_q10": 1.2,
    "microbe_methane_aerobic_oxidation_o2_ch4_ratio": 2.0,
    "microbe_methane_k_anaerobic_oxidation_ch4": 1.5,
    "microbe_methane_anaerobic_oxidation_q10": 1.2,
    "microbe_methane_aerobic_decomp_o2_c_ratio": 0.002,
    "microbe_methane_ch4_transport_threshold": 0.05,
    "microbe_methane_dom_diffusivity": 1.8e-7,
    "microbe_methane_aqueous_gas_diffusion_multiplier": 2.0,
    "microbe_methane_plant_transport_coefficient": 0.007,
    "microbe_methane_atmospheric_ch4_mixing_ratio": 1.7e-6,
}

LEGACY_MOLAR_CONCENTRATION_PARAMETERS = {
    "microbe_methane_k_acetate",
    "microbe_methane_k_acetate_prod_o2",
    "microbe_methane_k_acetogenesis_h2",
    "microbe_methane_k_acetogenesis_co2",
    "microbe_methane_k_h2_methanogenesis_h2",
    "microbe_methane_k_h2_methanogenesis_co2",
    "microbe_methane_k_acetoclastic_methanogenesis_acetate",
    "microbe_methane_k_aerobic_oxidation_ch4",
    "microbe_methane_k_aerobic_oxidation_o2",
    "microbe_methane_k_anaerobic_oxidation_ch4",
    "microbe_methane_ch4_transport_threshold",
}

LEGACY_VOLUMETRIC_RATE_PARAMETERS = {
    "microbe_methane_acetate_prod_max",
    "microbe_methane_acetogenesis_max",
}

LEGACY_SPECIFIC_RATE_PARAMETERS = {
    "microbe_methane_h2_methanogen_growth_rate",
    "microbe_methane_h2_methanogen_death_rate",
    "microbe_methane_acetate_methanogen_growth_rate",
    "microbe_methane_acetate_methanogen_death_rate",
    "microbe_methane_aerobic_methanotroph_growth_rate",
    "microbe_methane_aerobic_methanotroph_death_rate",
    "microbe_methane_anaerobic_methanotroph_growth_rate",
    "microbe_methane_anaerobic_methanotroph_death_rate",
}


def executable_equivalent_value(name: str, raw_value: float) -> float:
    """Convert an input-file number to the explicit units used by ELM."""
    if name in LEGACY_MOLAR_CONCENTRATION_PARAMETERS:
        return raw_value * 1000.0
    if name in LEGACY_VOLUMETRIC_RATE_PARAMETERS:
        return raw_value * 1000.0 * 86400.0
    if name in LEGACY_SPECIFIC_RATE_PARAMETERS:
        return raw_value * 86400.0
    return raw_value


class Phase3StateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.document = phase3.load_manifest(MANIFEST)
        cls.parameters = phase3.parameter_map(cls.document)
        cls.parameter_source = PARAMS_FORTRAN.read_text(encoding="utf-8")
        cls.state_source = STATE_FORTRAN.read_text(encoding="utf-8")

    def test_parameter_manifest_exactly_matches_fortran_reader(self) -> None:
        manifest_names = set(self.parameters)
        read_names = set(re.findall(
            r"call read_scalar\(ncid, '([^']+)'", self.parameter_source
        ))
        self.assertEqual(len(manifest_names), 73)
        self.assertEqual(
            self.document["schema"],
            "elm_microbe_methane_phase3_aqueous_transport_v1",
        )
        self.assertEqual(read_names, manifest_names)
        for name, (source_id, _value, units, provenance) in self.parameters.items():
            self.assertTrue(name.startswith("microbe_methane_"))
            self.assertTrue(source_id)
            self.assertTrue(units)
            self.assertTrue(provenance)
        self.assertIn("REFERENCE TEST VALUES ONLY", self.document["warning"])

        design_names = set(re.findall(
            r"`(microbe_methane_[a-z0-9_]+)`",
            DESIGN_DOCUMENT.read_text(encoding="utf-8"),
        ))
        self.assertTrue(manifest_names <= design_names)

    def test_temperature_unit_defects_are_explicitly_corrected(self) -> None:
        aom = self.parameters["microbe_methane_aom_t_ref"]
        thaw = self.parameters["microbe_methane_transport_thaw_threshold"]
        self.assertAlmostEqual(aom[1], 286.65)
        self.assertAlmostEqual(thaw[1], 273.05)
        self.assertEqual(aom[3], "clm_microbe_9c2e0a_source_bug_corrected_kelvin")
        self.assertEqual(thaw[3], "clm_microbe_9c2e0a_source_bug_corrected_kelvin")

    def test_runtime_file_values_are_converted_to_elm_declared_units(self) -> None:
        self.assertEqual(
            self.document["source_revision"],
            "CLM-Microbe 9c2e0a048bb3799669d32b91e6cc76efb36d4b75",
        )
        self.assertEqual(
            self.document["source_parameter_file"],
            "inputdata/lnd/clm2/paramdata/microbepar_in",
        )
        self.assertEqual(len(CLM_MICROBE_RUNTIME_VALUES), 38)
        for name, raw_value in CLM_MICROBE_RUNTIME_VALUES.items():
            _source_id, value, _units, provenance = self.parameters[name]
            expected = executable_equivalent_value(name, raw_value)
            self.assertAlmostEqual(value, expected)
            self.assertTrue(provenance.startswith("clm_microbe_9c2e0a_runtime_file"))

    def test_source_defaults_and_literals_use_elm_declared_units(self) -> None:
        expected = {
            "microbe_methane_mfg_biomass_min": 1.0e-15 * 12.011,
            "microbe_methane_h2_plant_transport_threshold": 0.000473 * 1000.0,
            "microbe_methane_acidification_coefficient": 4.2e-9 / 1000.0,
            "microbe_methane_acetate_feedback_half_saturation": 0.1 * 1000.0,
            "microbe_methane_h2_methanogenesis_co2_inhibition_scale": 9.2 * 1000.0,
            "microbe_methane_aom_o2_inhibition_scale": 4.6 * 1000.0,
            "microbe_methane_aerobic_acetate_oxidation_rate": 0.05 * 86400.0,
        }
        for name, converted_value in expected.items():
            _source_id, value, _units, provenance = self.parameters[name]
            self.assertAlmostEqual(value, converted_value)
            self.assertIn("legacy_", provenance)

    def test_state_is_partitioned_without_duplicate_decomp_pools(self) -> None:
        pointers = set(re.findall(
            r"real\(r8\), pointer :: ([a-z0-9_]+)\(:,:\)", self.state_source
        ))
        diagnostic_pointers = {
            "o2_stress_unsat_col", "o2_stress_sat_col",
            "dom_advective_flux_col", "dom_diffusive_flux_col",
        }
        self.assertTrue(diagnostic_pointers <= pointers)
        prognostic_pointers = pointers - diagnostic_pointers
        self.assertEqual(len(prognostic_pointers), 18)
        for forbidden in ("dom", "bacteria", "fungi"):
            self.assertFalse(any(forbidden in name for name in prognostic_pointers))
        for state in ("acetate", "acetate_methanogen", "h2_methanogen",
                      "aerobic_methanotroph", "anaerobic_methanotroph",
                      "conc_ch4", "conc_o2", "conc_co2", "conc_h2"):
            self.assertIn(f"{state}_c_unsat_col" if not state.startswith("conc_")
                          else f"{state}_unsat_col", prognostic_pointers)
            self.assertIn(f"{state}_c_sat_col" if not state.startswith("conc_")
                          else f"{state}_sat_col", prognostic_pointers)
        self.assertIn("sat_fraction_previous_col(:)", self.state_source)

    def test_cold_start_uses_named_seed_and_zero_substrates_and_gases(self) -> None:
        self.assertIn(
            "biomass_seed = MicrobeMethaneParamsInst%mfg_biomass_min",
            self.state_source,
        )
        guild_assignments = re.findall(
            r"this%([a-z0-9_]*(?:methanogen|methanotroph)[a-z0-9_]*)"
            r"\(c,1:nlevdecomp\) = biomass_seed",
            self.state_source,
        )
        self.assertEqual(len(guild_assignments), 8)
        self.assertIn("this%sat_fraction_previous_col = 0._r8", self.state_source)
        self.assertNotIn("1.e-5_r8 * catomw", self.state_source)
        for name in ("acetate_c_unsat_col", "acetate_c_sat_col",
                     "conc_ch4_unsat_col", "conc_ch4_sat_col",
                     "conc_o2_unsat_col", "conc_o2_sat_col",
                     "conc_co2_unsat_col", "conc_co2_sat_col",
                     "conc_h2_unsat_col", "conc_h2_sat_col"):
            self.assertIn(f"this%{name}(c,1:nlevdecomp) = 0._r8", self.state_source)

    def test_history_and_restart_cover_every_prognostic_field(self) -> None:
        history_block = self.state_source.split("subroutine InitHistory", 1)[1].split(
            "end subroutine InitHistory", 1
        )[0]
        restart_block = self.state_source.split("subroutine Restart", 1)[1].split(
            "end subroutine Restart", 1
        )[0]
        history_names = set(re.findall(r"call add_state\('([^']+)'", history_block))
        restart_names = set(re.findall(r"call restart_state\('([^']+)'", restart_block))
        self.assertEqual(history_names, STATE_HISTORY_NAMES)
        self.assertEqual(restart_names, STATE_HISTORY_NAMES)
        self.assertIn("MM_SAT_FRACTION_PREVIOUS", restart_block)
        self.assertIn("nsrest /= nsrStartup", restart_block)
        self.assertIn("A legacy restart cannot be continued or branched", restart_block)
        self.assertIn("default='inactive'", history_block)

    def test_lifecycle_is_conditional_and_legacy_backend_is_retained(self) -> None:
        read_source = READ_PARAMS_FORTRAN.read_text(encoding="utf-8")
        restart_source = RESTART_FORTRAN.read_text(encoding="utf-8")
        driver_source = DRIVER_FORTRAN.read_text(encoding="utf-8")
        self.assertIn("call readMicrobeMethaneParams(ncid)", read_source)
        self.assertEqual(restart_source.count("microbe_methane_vars%Restart"), 3)
        self.assertGreaterEqual(restart_source.count("ch4_vars%restart"), 3)
        self.assertIn("call CH4", driver_source)
        self.assertNotIn("microbe_methane_vars%Step", driver_source)
        self.assertIn(
            "if (.not. use_microbe_methane .or. use_legacy_ch4_with_microbe) return",
            self.state_source,
        )

    def test_parameter_file_injection_round_trip(self) -> None:
        try:
            import netCDF4  # type: ignore
        except ImportError:
            self.skipTest("netCDF4 is unavailable; Docker runs this test")

        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            input_path = temporary / "phase2.nc"
            output_path = temporary / "phase3.nc"
            with netCDF4.Dataset(input_path, "w") as dataset:
                dataset.createDimension("pft", 4)
                dataset.createDimension("allpfts", 1)
                sentinel = dataset.createVariable("k_dom", "f8", ("pft",))
                sentinel[:] = [0.1, 0.2, 0.3, 0.4]

            phase3.inject_parameters(input_path, output_path, MANIFEST)
            phase3.validate_parameter_file(output_path, MANIFEST)

            with netCDF4.Dataset(output_path) as dataset:
                variable = dataset.variables["microbe_methane_mfg_biomass_min"]
                self.assertEqual(variable.dimensions, ("allpfts",))
                self.assertEqual(variable.legacy_name, "MFGbiomin")
                self.assertEqual(variable.science_status, "phase3_reference_test_only")
                self.assertEqual(list(dataset.variables["k_dom"][:]), [0.1, 0.2, 0.3, 0.4])

            with netCDF4.Dataset(output_path, "a") as dataset:
                dataset.variables["microbe_methane_mfg_biomass_min"][:] = 2.0e-15
            with self.assertRaisesRegex(SystemExit, "microbe_methane_mfg_biomass_min"):
                phase3.validate_parameter_file(output_path, MANIFEST)


if __name__ == "__main__":
    unittest.main()
