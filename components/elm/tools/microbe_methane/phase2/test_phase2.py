#!/usr/bin/env python3
"""Unit tests for the Phase 2 microbial decomposition cascade."""

from __future__ import annotations

import re
import sys
import tempfile
import unittest
from pathlib import Path


PHASE2_DIR = Path(__file__).resolve().parent
ELM_DIR = PHASE2_DIR.parents[2]
sys.path.insert(0, str(PHASE2_DIR))

import phase2  # noqa: E402


DEFINITION = PHASE2_DIR / "phase2_test_parameters.json"
FORTRAN = ELM_DIR / "src" / "biogeochem" / "MicrobeDecompMod.F90"
TRANSPORT_FORTRAN = ELM_DIR / "src" / "biogeochem" / "SoilLittVertTranspMod.F90"
CONTROL_FORTRAN = ELM_DIR / "src" / "main" / "controlMod.F90"
BUILD_NAMELIST = ELM_DIR / "bld" / "ELMBuildNamelist.pm"
CH4_FORTRAN = ELM_DIR / "src" / "biogeochem" / "CH4Mod.F90"
READ_PARAMS_FORTRAN = ELM_DIR / "src" / "main" / "readParamsMod.F90"
RESTART_FORTRAN = ELM_DIR / "src" / "main" / "restFileMod.F90"

EXPECTED_PFT_PARAMETERS = {
    "k_dom", "k_bacteria", "k_fungi",
    "m_rf_s1m", "m_rf_s2m", "m_rf_s3m", "m_rf_s4m",
    "m_batm_f", "m_fatm_f", "m_bdom_f", "m_fdom_f",
    "m_bs1_f", "m_bs2_f", "m_bs3_f",
    "m_fs1_f", "m_fs2_f", "m_fs3_f",
    "m_domb_f", "m_domf_f",
    "m_doms1_f", "m_doms2_f", "m_doms3_f",
    "cn_bacteria", "cn_fungi",
}

EXPECTED_SCALAR_PARAMETERS = {
    "cn_dom", "bacteria_pool_cn", "fungi_pool_cn",
    "cp_bacteria", "cp_fungi", "cp_dom",
    "CUEmax", "microbe_cue_cn_target", "microbe_allocation_cn_exponent",
    "bacteria_initial_c", "fungi_initial_c", "dom_initial_c",
    "dom_som_diffusion_multiplier",
    "microbe_som2_q10", "microbe_som3_q10", "microbe_som4_q10",
    "microbe_dom_q10",
    "l1dom_f", "l2dom_f", "l3dom_f",
    "s1dom_f", "s2dom_f", "s3dom_f", "s4dom_f",
    "l1s1_f", "l2s2_f", "l3s3_f",
    "s1s2_f", "s2s3_f", "s3s4_f",
}


class Phase2CascadeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.parameters, cls.document = phase2.load_parameters(DEFINITION)
        cls.pool_cn = [0.0, 90.0, 90.0, 90.0, 500.0,
                       12.0, 12.0, 10.0, 10.0, 5.0, 15.0, 10.0]
        cls.pool_cp = [0.0, 900.0, 900.0, 900.0, 5000.0,
                       360.0, 360.0, 500.0, 500.0, 150.0, 450.0, 300.0]
        cls.path, cls.respiration = phase2.build_pathways(
            cls.parameters, cls.pool_cn
        )

    def test_parameter_manifest_is_complete_and_named(self) -> None:
        self.assertEqual(set(self.document["pft"]), EXPECTED_PFT_PARAMETERS)
        self.assertEqual(set(self.document["scalar"]), EXPECTED_SCALAR_PARAMETERS)
        self.assertIn("test", self.document["warning"].lower())
        provenance = phase2.parameter_provenance(self.document)
        self.assertEqual(
            set(provenance), EXPECTED_PFT_PARAMETERS | EXPECTED_SCALAR_PARAMETERS
        )
        self.assertEqual(provenance["k_dom"], "clm_spruce_unverified_text_value")
        self.assertEqual(provenance["CUEmax"], "clm_spruce_source_literal")
        self.assertEqual(provenance["microbe_dom_q10"],
                         "clm_spruce_source_literal")
        self.assertEqual(provenance["cn_bacteria"],
                         "test_assumption_from_pool_cn_literal")
        self.assertEqual(provenance["cp_bacteria"], "new_cnp_test_assumption")

        source = FORTRAN.read_text(encoding="utf-8")
        pft_reads = set(re.findall(
            r"call read_pft_parameter\(ncid, '([^']+)'", source
        ))
        scalar_reads = set(re.findall(
            r"call read_scalar_parameter\(ncid, '([^']+)'", source
        ))
        self.assertEqual(pft_reads, EXPECTED_PFT_PARAMETERS)
        self.assertEqual(scalar_reads, EXPECTED_SCALAR_PARAMETERS)

    def test_fortran_and_oracle_transition_topology_match(self) -> None:
        source = FORTRAN.read_text(encoding="utf-8")
        calls = re.findall(
            r"call set_transition\(\s*(\d+),\s*'([^']+)',\s*([^,]+),\s*([^\)]+)\)",
            source,
        )
        self.assertEqual(len(calls), 47)

        index_value = {
            "i_atm": 0,
            "i_met_lit": 1,
            "i_cel_lit": 2,
            "i_lig_lit": 3,
            "i_cwd": 4,
            "i_soil1": 5,
            "i_soil2": 6,
            "i_soil3": 7,
            "i_soil4": 8,
            "i_bacteria": 9,
            "i_fungi": 10,
            "i_dom": 11,
        }
        parsed = []
        for expected_index, (index, name, donor, receiver) in enumerate(calls, 1):
            self.assertEqual(int(index), expected_index)
            parsed.append((name, index_value[donor.strip()], index_value[receiver.strip()]))
        self.assertEqual(tuple(parsed), phase2.TRANSITIONS)

    def test_each_donor_path_closes_and_fractions_are_bounded(self) -> None:
        self.assertEqual(len(self.path), 47)
        self.assertEqual(len(self.respiration), 47)
        for value in self.path + self.respiration:
            self.assertGreaterEqual(value, 0.0)
            self.assertLessEqual(value, 1.0)
        for total in phase2.donor_sums(self.path).values():
            self.assertAlmostEqual(total, 1.0, places=14)

    def test_feedback_paths_connect_all_three_new_pools(self) -> None:
        routes = {(phase2.POOLS[donor], phase2.POOLS[receiver])
                  for _, donor, receiver in phase2.TRANSITIONS}
        self.assertTrue({("litr1", "dom"), ("soil4", "dom")} <= routes)
        self.assertTrue({("dom", "bacteria"), ("dom", "fungi")} <= routes)
        self.assertTrue({("bacteria", "dom"), ("fungi", "dom")} <= routes)
        self.assertTrue({("bacteria", "soil4"), ("fungi", "soil4")} <= routes)
        self.assertTrue({("bacteria", "atmosphere"),
                         ("fungi", "atmosphere")} <= routes)

        # Solubilization is a conservative transfer. The legacy source set
        # respiration to one on these paths, which would discard all receiver C
        # under the current ELM cascade semantics.
        self.assertEqual(self.respiration[40:47], [0.0] * 7)

    def test_closed_box_conserves_c_n_p_and_stays_nonnegative(self) -> None:
        carbon = [0.0, 80.0, 120.0, 100.0, 30.0,
                  200.0, 400.0, 800.0, 1200.0, 5.0, 3.0, 20.0]
        nitrogen = [0.0] + [carbon[i] / self.pool_cn[i]
                            for i in range(1, len(carbon))]
        phosphorus = [0.0] + [carbon[i] / self.pool_cp[i]
                              for i in range(1, len(carbon))]
        mineral_n = 100.0
        mineral_p = 20.0
        respiration_total = 0.0
        rates = [0.0, 0.020, 0.015, 0.010, 0.001,
                 0.010, 0.005, 0.002, 0.001, 0.10, 0.10, 0.02]

        initial_c = sum(carbon)
        initial_n = sum(nitrogen) + mineral_n
        initial_p = sum(phosphorus) + mineral_p

        for _ in range(365):
            (carbon, nitrogen, phosphorus, mineral_n, mineral_p,
             respiration_total) = phase2.closed_box_step(
                carbon, nitrogen, phosphorus, mineral_n, mineral_p,
                respiration_total, rates, self.pool_cn, self.pool_cp,
                self.path, self.respiration,
            )
            self.assertGreaterEqual(min(carbon), -1.0e-12)
            self.assertGreaterEqual(min(nitrogen), -1.0e-12)
            self.assertGreaterEqual(min(phosphorus), -1.0e-12)
            self.assertGreaterEqual(mineral_n, -1.0e-12)
            self.assertGreaterEqual(mineral_p, -1.0e-12)

        self.assertAlmostEqual(sum(carbon) + respiration_total, initial_c, places=9)
        self.assertAlmostEqual(sum(nitrogen) + mineral_n, initial_n, places=9)
        self.assertAlmostEqual(sum(phosphorus) + mineral_p, initial_p, places=9)
        self.assertGreater(carbon[9], 0.0)
        self.assertGreater(carbon[10], 0.0)
        self.assertGreater(carbon[11], 0.0)

    def test_dom_transport_is_distinct_and_parameterized(self) -> None:
        source = TRANSPORT_FORTRAN.read_text(encoding="utf-8")
        self.assertIn("is_dissolved(s)", source)
        self.assertIn("dom_som_diffusion_multiplier", source)
        self.assertEqual(self.parameters["dom_som_diffusion_multiplier"], 10.0)

    def test_site_mode_retains_established_methane_coupling(self) -> None:
        build_source = BUILD_NAMELIST.read_text(encoding="utf-8")
        control_source = CONTROL_FORTRAN.read_text(encoding="utf-8")
        ch4_source = CH4_FORTRAN.read_text(encoding="utf-8")
        read_params_source = READ_PARAMS_FORTRAN.read_text(encoding="utf-8")
        restart_source = RESTART_FORTRAN.read_text(encoding="utf-8")

        self.assertIn(
            "Phase 2 microbial decomposition requires use_lch4=.true.",
            build_source,
        )
        self.assertIn(
            "to retain methane oxygen and nitrification/denitrification coupling.",
            control_source,
        )
        self.assertIn("if (use_lch4) then", ch4_source)
        self.assertNotIn("use_lch4 .and. .not. use_microbe_methane", ch4_source)
        self.assertIn("if (use_lch4) then", read_params_source)
        self.assertNotIn("use_lch4 .and. .not. use_microbe_methane", read_params_source)
        # Phase 3 adds its independent state lifecycle while preserving the
        # established backend as a permanent capability.
        self.assertIn("microbe_methane_vars%Restart", restart_source)
        self.assertIn("call ch4_vars%restart", restart_source)

    def test_parameter_file_injection_round_trip(self) -> None:
        try:
            import netCDF4  # type: ignore
        except ImportError:
            self.skipTest("netCDF4 is unavailable; Docker runs this test")

        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            input_path = temporary / "base.nc"
            output_path = temporary / "phase2.nc"
            with netCDF4.Dataset(input_path, "w") as dataset:
                dataset.createDimension("pft", 4)
                dataset.createDimension("allpfts", 1)
                sentinel = dataset.createVariable("sentinel", "f8", ("allpfts",))
                sentinel[:] = [42.0]

            phase2.inject_parameters(input_path, output_path, DEFINITION)
            phase2.validate_parameter_file(output_path, DEFINITION)

            with netCDF4.Dataset(output_path) as dataset:
                self.assertEqual(dataset.variables["k_dom"].dimensions, ("pft",))
                self.assertEqual(
                    dataset.variables["dom_som_diffusion_multiplier"].dimensions,
                    ("allpfts",),
                )
                self.assertEqual(float(dataset.variables["sentinel"][0]), 42.0)
                self.assertEqual(
                    dataset.variables["cp_dom"].source_provenance,
                    "new_cnp_test_assumption",
                )


if __name__ == "__main__":
    unittest.main()
