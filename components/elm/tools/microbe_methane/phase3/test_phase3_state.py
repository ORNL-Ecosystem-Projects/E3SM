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
        self.assertEqual(len(manifest_names), 64)
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
        self.assertEqual(aom[3], "corrected_celsius_kelvin_source_bug")
        self.assertEqual(thaw[3], "corrected_celsius_kelvin_source_bug")

    def test_state_is_partitioned_without_duplicate_decomp_pools(self) -> None:
        pointers = set(re.findall(
            r"real\(r8\), pointer :: ([a-z0-9_]+)\(:,:\)", self.state_source
        ))
        self.assertEqual(len(pointers), 18)
        for forbidden in ("dom", "bacteria", "fungi"):
            self.assertFalse(any(forbidden in name for name in pointers))
        for state in ("acetate", "acetate_methanogen", "h2_methanogen",
                      "aerobic_methanotroph", "anaerobic_methanotroph",
                      "conc_ch4", "conc_o2", "conc_co2", "conc_h2"):
            self.assertIn(f"{state}_c_unsat_col" if not state.startswith("conc_")
                          else f"{state}_unsat_col", pointers)
            self.assertIn(f"{state}_c_sat_col" if not state.startswith("conc_")
                          else f"{state}_sat_col", pointers)
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
        self.assertIn("if (.not. use_microbe_methane) return", self.state_source)

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


if __name__ == "__main__":
    unittest.main()
