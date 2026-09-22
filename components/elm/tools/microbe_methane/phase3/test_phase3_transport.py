#!/usr/bin/env python3
"""Conservation and boundary tests for Phase 3 Step 3 transport."""

from __future__ import annotations

import unittest
from pathlib import Path

import phase3
from transport_oracle import (
    aerenchyma_transport,
    aqueous_tracer_transport,
    bulk_concentration,
    dom_cnp_relaxation,
    effective_aqueous_diffusivity,
    gas_transport,
    implicit_vertical_diffusion,
    legacy_fickian_gas_diffusivity,
    methane_ebullition,
    repartition,
    topounit_lateral_diffusion,
    transport_residual,
    vertical_diffusion,
)


PHASE3_DIR = Path(__file__).resolve().parent
ELM_DIR = PHASE3_DIR.parents[2]
MANIFEST = PHASE3_DIR / "phase3_reference_parameters.json"
TRANSPORT = ELM_DIR / "src" / "biogeochem" / "MicrobeGasTransportMod.F90"
STATE = ELM_DIR / "src" / "biogeochem" / "MicrobeMethaneMod.F90"
STATE_UPDATE = ELM_DIR / "src" / "biogeochem" / "MicrobeMethaneStateUpdateMod.F90"
DRIVER = ELM_DIR / "src" / "main" / "elm_driver.F90"
VARCTL = ELM_DIR / "src" / "main" / "elm_varctl.F90"
CONTROL = ELM_DIR / "src" / "main" / "controlMod.F90"
NAMELIST_DEFINITION = ELM_DIR / "bld" / "namelist_files" / "namelist_definition.xml"
BUILD_NAMELIST = ELM_DIR / "bld" / "ELMBuildNamelist.pm"
ALLOCATION = ELM_DIR / "src" / "biogeochem" / "AllocationMod.F90"
DESIGN = ELM_DIR / "docs" / "dev-guide" / "spruce-microbe-methane-design.md"


class Phase3TransportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        document = phase3.load_manifest(MANIFEST)
        cls.parameters = {
            name: float(details[1])
            for name, details in phase3.parameter_map(document).items()
        }
        cls.transport_source = TRANSPORT.read_text(encoding="utf-8")
        cls.state_source = STATE.read_text(encoding="utf-8")
        cls.state_update_source = STATE_UPDATE.read_text(encoding="utf-8")
        cls.varctl_source = VARCTL.read_text(encoding="utf-8")
        cls.control_source = CONTROL.read_text(encoding="utf-8")
        cls.namelist_definition = NAMELIST_DEFINITION.read_text(encoding="utf-8")
        cls.build_namelist = BUILD_NAMELIST.read_text(encoding="utf-8")
        cls.allocation_source = ALLOCATION.read_text(encoding="utf-8")

    def test_repartition_conserves_bulk_at_limits_and_intermediate_fractions(self) -> None:
        for old_fraction, new_fraction in (
            (0.0, 0.0),
            (0.0, 1.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.2, 0.7),
            (0.8, 0.3),
        ):
            with self.subTest(old=old_fraction, new=new_fraction):
                unsaturated, saturated = repartition(
                    old_fraction, new_fraction, 2.5, 9.5
                )
                before = bulk_concentration(old_fraction, 2.5, 9.5)
                after = bulk_concentration(new_fraction, unsaturated, saturated)
                self.assertAlmostEqual(after, before, places=14)
                self.assertGreaterEqual(unsaturated, 0.0)
                self.assertGreaterEqual(saturated, 0.0)

    def test_state_repartition_covers_every_partitioned_pool(self) -> None:
        self.assertIn("procedure, public  :: Repartition", self.state_source)
        for field in (
            "acetate_c",
            "acetate_methanogen_c",
            "h2_methanogen_c",
            "aerobic_methanotroph_c",
            "anaerobic_methanotroph_c",
            "conc_ch4",
            "conc_o2",
            "conc_co2",
            "conc_h2",
        ):
            self.assertIn(f"this%{field}_unsat_col(c,j)", self.state_source)
            self.assertIn(f"this%{field}_sat_col(c,j)", self.state_source)
        self.assertIn("this%sat_fraction_previous_col(c) = new_fraction", self.state_source)

    def test_zero_gradient_diffusion_is_zero(self) -> None:
        interfaces, tendency, surface_flux = vertical_diffusion(
            [3.0, 3.0, 3.0], [0.1, 0.2, 0.3], [1.0e-9] * 3, 3.0, 2.0e-5, 1800.0
        )
        self.assertEqual(interfaces, [0.0, 0.0, 0.0, 0.0])
        self.assertEqual(tendency, [0.0, 0.0, 0.0])
        self.assertEqual(surface_flux, 0.0)

    def test_two_layer_diffusion_matches_finite_volume_solution(self) -> None:
        interfaces, tendency, surface_flux = vertical_diffusion(
            [2.0, 0.0], [1.0, 1.0], [1.0e-6, 1.0e-6], 2.0, 0.0, 10.0
        )
        self.assertAlmostEqual(interfaces[1], 2.0e-6)
        self.assertAlmostEqual(tendency[0], -2.0e-6)
        self.assertAlmostEqual(tendency[1], 2.0e-6)
        self.assertEqual(surface_flux, 0.0)
        self.assertAlmostEqual(transport_residual(tendency, [1.0, 1.0], surface_flux), 0.0)

    def test_dom_cnp_relaxation_matches_equal_layer_timescale(self) -> None:
        rate = 1.8e-7
        dt = 1800.0
        updated_c, updated_n, updated_p, residuals = dom_cnp_relaxation(
            [10.0, 2.0], [1.0, 0.2], [0.5, 0.1], [0.1, 0.1], [rate, rate], dt
        )
        expected_difference = 8.0 / (1.0 + 2.0 * rate * dt)
        self.assertAlmostEqual(updated_c[0] + updated_c[1], 12.0)
        self.assertAlmostEqual(updated_c[0] - updated_c[1], expected_difference)
        for carbon, nitrogen, phosphorus in zip(updated_c, updated_n, updated_p):
            self.assertAlmostEqual(carbon / nitrogen, 10.0)
            self.assertAlmostEqual(carbon / phosphorus, 20.0)
        for residual in residuals:
            self.assertAlmostEqual(residual, 0.0, places=14)

    def test_dom_cnp_relaxation_conserves_unequal_layer_inventories(self) -> None:
        thickness = [0.05, 0.15, 0.4]
        initial = ([30.0, 3.0, 0.1], [2.0, 0.4, 0.02], [0.5, 0.1, 0.005])
        updated_c, updated_n, updated_p, residuals = dom_cnp_relaxation(
            *initial, thickness, [1.8e-7, 1.7e-7, 1.6e-7], 86400.0
        )
        self.assertTrue(all(value >= 0.0 for value in updated_c + updated_n + updated_p))
        for residual in residuals:
            self.assertAlmostEqual(residual, 0.0, places=14)

    def test_dom_relaxation_is_explicit_default_off_parity_option(self) -> None:
        self.assertIn('id="use_clm_microbe_dom_relaxation"', self.namelist_definition)
        self.assertIn('value=".false."', self.namelist_definition.split(
            'id="use_clm_microbe_dom_relaxation"', 1
        )[1].split("</entry>", 1)[0])
        self.assertIn("advanceMicrobeMethaneDOMRelaxation", self.state_update_source)
        self.assertIn("if (use_clm_microbe_dom_relaxation) then", self.state_source)
        self.assertIn("use_clm_microbe_dom_relaxation=.true. requires", self.control_source)
        self.assertIn("use_clm_microbe_dom_relaxation=.true. requires", self.build_namelist)

    def test_aqueous_diffusion_is_grid_aware_and_conservative(self) -> None:
        updated, advective, diffusive, _tendency, export, residual = (
            aqueous_tracer_transport(
                [10.0, 0.0], [0.1, 0.3], [0.5, 0.5],
                [5.0e-11, 5.0e-11], [0.0, 0.0, 0.0], 1.0, 1.0e-4, 86400.0
            )
        )
        self.assertLess(updated[0], 10.0)
        self.assertGreater(updated[1], 0.0)
        self.assertEqual(advective, [0.0, 0.0, 0.0])
        self.assertGreater(diffusive[1], 0.0)
        self.assertEqual(export, 0.0)
        self.assertAlmostEqual(residual, 0.0, places=14)

    def test_aqueous_advection_upwinds_and_exports_at_bottom(self) -> None:
        dt = 86400.0
        updated, advective, _diffusive, _tendency, export, residual = (
            aqueous_tracer_transport(
                [5.0, 0.0], [0.1, 0.1], [0.5, 0.5], [0.0, 0.0],
                [0.0, 1.0e-7, 1.0e-7], 1.0, 1.0e-4, dt
            )
        )
        self.assertGreater(advective[1], 0.0)
        self.assertGreater(advective[2], 0.0)
        self.assertGreater(export, 0.0)
        self.assertTrue(all(value >= 0.0 for value in updated))
        self.assertAlmostEqual(residual, 0.0, places=14)

    def test_aqueous_operator_preserves_uniform_dom_stoichiometry(self) -> None:
        inputs = [12.0, 3.0, 0.5]
        common = ([0.08, 0.2, 0.4], [0.35, 0.45, 0.5],
                  [4.0e-11, 7.0e-11, 9.0e-11],
                  [0.0, 2.0e-8, -1.0e-8, 3.0e-8], 1.0, 1.0e-4, 1800.0)
        carbon, *_ = aqueous_tracer_transport(inputs, *common)
        nitrogen, *_ = aqueous_tracer_transport([value / 10.0 for value in inputs], *common)
        phosphorus, *_ = aqueous_tracer_transport([value / 300.0 for value in inputs], *common)
        for c_value, n_value, p_value in zip(carbon, nitrogen, phosphorus):
            self.assertAlmostEqual(c_value / n_value, 10.0)
            self.assertAlmostEqual(c_value / p_value, 300.0)

    def test_aqueous_transport_immobilizes_dry_layers(self) -> None:
        updated, advective, diffusive, _tendency, export, residual = (
            aqueous_tracer_transport(
                [5.0, 0.0], [0.1, 0.1], [1.0e-8, 0.5], [1.0e-9, 1.0e-9],
                [0.0, 1.0e-6, 0.0], 1.0, 1.0e-4, 1800.0
            )
        )
        self.assertEqual(updated, [5.0, 0.0])
        self.assertEqual(advective, [0.0, 0.0, 0.0])
        self.assertEqual(diffusive, [0.0, 0.0, 0.0])
        self.assertEqual(export, 0.0)
        self.assertEqual(residual, 0.0)

    def test_aqueous_budget_check_has_transport_specific_roundoff_tolerance(self) -> None:
        self.assertIn("state_tolerance = 1.e-12_r8", self.state_update_source)
        self.assertIn("aqueous_budget_tolerance = 1.e-10_r8", self.state_update_source)
        self.assertIn(
            "aqueousResidualIsClosed(residual, initial_inventory, final_inventory)",
            self.state_update_source,
        )

    def test_physical_aqueous_transport_is_default_off_and_exclusive(self) -> None:
        entry = self.namelist_definition.split(
            'id="use_microbe_aqueous_transport"', 1
        )[1].split("</entry>", 1)[0]
        self.assertIn('value=".false."', entry)
        self.assertIn("advanceMicrobeAqueousTracerTransport", self.state_update_source)
        self.assertIn("use_microbe_aqueous_transport and", self.control_source)
        self.assertIn("are mutually exclusive", self.build_namelist)
        self.assertIn("DOM carbon isotopes", self.control_source)
        self.assertIn("DOM carbon isotopes", self.build_namelist)

    def test_peatland_root_access_is_pft_gated_and_normalized(self) -> None:
        self.assertNotIn("use_microbe_unsaturated_root_n_access", self.varctl_source)
        self.assertIn("use_peatland_roots", self.allocation_source)
        self.assertIn("veg_vp%nonvascular(ivt(p)) < 0.5_r8", self.allocation_source)
        self.assertIn("h2osoi_vol(c,j) < watsat(c,j)", self.allocation_source)
        self.assertIn("adaptive_profile_sum", self.allocation_source)
        self.assertIn("use TopounitType        , only : top_pp", self.allocation_source)
        self.assertIn(
            "top_pp%peat_depth(col_pp%topounit(c)) > 0._r8",
            self.allocation_source,
        )

    def test_diffusion_surface_budget_and_donor_limiting(self) -> None:
        dt = 86400.0
        _interfaces, tendency, surface_flux = vertical_diffusion(
            [1.0e-8, 0.0, 0.0],
            [0.1, 0.1, 0.1],
            [1.0, 1.0, 1.0],
            0.0,
            1.0,
            dt,
        )
        for concentration, rate in zip([1.0e-8, 0.0, 0.0], tendency):
            self.assertGreaterEqual(concentration + dt * rate, -1.0e-20)
        self.assertAlmostEqual(
            transport_residual(tendency, [0.1, 0.1, 0.1], surface_flux),
            0.0,
            places=20,
        )

    def test_capacity_aware_diffusion_preserves_phase_equilibrium(self) -> None:
        updated, interfaces, tendency, surface_flux = implicit_vertical_diffusion(
            concentration=[0.02, 0.002],
            layer_thickness=[0.1, 0.2],
            effective_diffusivity=[1.0e-6, 1.0e-9],
            transport_capacity=[0.2, 0.02],
            atmospheric_mobile_concentration=0.1,
            surface_conductance=1.0e-5,
            dt=1800.0,
        )
        self.assertEqual(updated, [0.02, 0.002])
        self.assertEqual(interfaces, [0.0, 0.0, 0.0])
        self.assertEqual(tendency, [0.0, 0.0])
        self.assertEqual(surface_flux, 0.0)

    def test_implicit_surface_uptake_is_positive_and_conservative(self) -> None:
        dt = 1800.0
        updated, _interfaces, tendency, surface_flux = implicit_vertical_diffusion(
            concentration=[0.0],
            layer_thickness=[0.1],
            effective_diffusivity=[1.0e-6],
            transport_capacity=[0.2],
            atmospheric_mobile_concentration=0.1,
            surface_conductance=1.0e-5,
            dt=dt,
        )
        self.assertGreater(updated[0], 0.0)
        self.assertGreaterEqual(updated[0], -1.0e-20)
        self.assertLess(surface_flux, 0.0)
        self.assertAlmostEqual(updated[0], dt * tendency[0])
        self.assertAlmostEqual(
            transport_residual(tendency, [0.1], surface_flux), 0.0, places=18
        )

    def test_aerenchyma_sign_units_and_inventory_limit(self) -> None:
        dt = 2.0
        layer_flux, tendency, surface_flux = aerenchyma_transport(
            [1.0, 0.5], [0.0, 1.0], [0.2, 0.4], [10.0, 0.1], dt
        )
        self.assertEqual(layer_flux[0], 0.5)
        self.assertLess(layer_flux[1], 0.0)
        self.assertEqual(tendency, [-value for value in layer_flux])
        self.assertGreaterEqual(1.0 + dt * tendency[0], 0.0)
        self.assertAlmostEqual(
            transport_residual(tendency, [0.2, 0.4], surface_flux), 0.0
        )

    def test_aerenchyma_influx_cannot_cross_atmospheric_equilibrium(self) -> None:
        dt = 1800.0
        layer_flux, tendency, _surface_flux = aerenchyma_transport(
            [0.0], [0.3], [0.01], [4.0], dt, apply_donor_limit=False
        )
        self.assertAlmostEqual(layer_flux[0], -0.3 / dt)
        self.assertAlmostEqual(0.0 + dt * tendency[0], 0.3)

    def test_one_way_plant_threshold_cannot_import_atmospheric_methane(self) -> None:
        dt = 1800.0
        layer_flux, tendency, surface_flux = aerenchyma_transport(
            [0.0, 0.06],
            [3.0e-6, 3.0e-6],
            [0.1, 0.1],
            [4.0, 4.0],
            dt,
            apply_donor_limit=False,
            minimum_emission_concentration=[0.05, 0.05],
            allow_influx=False,
        )
        self.assertEqual(layer_flux[0], 0.0)
        self.assertEqual(tendency[0], 0.0)
        self.assertGreater(layer_flux[1], 0.0)
        self.assertGreater(surface_flux, 0.0)

    def test_ebullition_threshold_and_budget(self) -> None:
        dt = 10.0
        loss, tendency, surface_flux = methane_ebullition(
            [2.0, 0.5], [1.0, 1.0], [1.0, 1.0], [0.2, 0.4], dt
        )
        self.assertAlmostEqual(loss[0], 0.1)
        self.assertEqual(loss[1], 0.0)
        self.assertAlmostEqual(2.0 + dt * tendency[0], 1.0)
        self.assertAlmostEqual(
            transport_residual(tendency, [0.2, 0.4], surface_flux), 0.0
        )

    def test_effective_diffusivity_uses_named_controls(self) -> None:
        reference = 2.0e-9
        temperature = self.parameters["microbe_methane_aqueous_diffusion_t_ref"]
        expected = (
            reference
            * self.parameters["microbe_methane_aqueous_gas_diffusion_multiplier"]
            * 0.25
        )
        self.assertAlmostEqual(
            effective_aqueous_diffusivity(reference, temperature, 0.25, self.parameters),
            expected,
        )
        self.assertEqual(
            effective_aqueous_diffusivity(reference, temperature, 0.0, self.parameters),
            0.0,
        )

    def test_default_legacy_fickian_coefficients_match_executed_source(self) -> None:
        temperature = self.parameters["microbe_methane_aqueous_diffusion_t_ref"]
        multiplier = self.parameters[
            "microbe_methane_aqueous_gas_diffusion_multiplier"
        ]
        for gas, source_value in enumerate((1.49e-5, 2.10e-5, 1.92e-5, 4.50e-5)):
            with self.subTest(gas=gas):
                self.assertAlmostEqual(
                    legacy_fickian_gas_diffusivity(gas, temperature, self.parameters),
                    source_value * 1.0e-3 * multiplier,
                )

    def test_two_topounit_lateral_diffusion_applies_shared_area_weight(self) -> None:
        diffusivity = 2.0e-4
        updated, residual = topounit_lateral_diffusion(
            concentration=[[4.0], [0.0]],
            storage_weight=[[0.75], [0.25]],
            layer_top_elevation=[[1.0], [1.0]],
            layer_bottom_elevation=[[0.0], [0.0]],
            effective_diffusivity=[[diffusivity], [diffusivity]],
            edges=[(0, 1, 1.0)],
            dt=1.0,
        )
        transfer = (4.0 - 0.0) * diffusivity * min(0.75, 0.25)
        self.assertAlmostEqual(updated[0][0], (4.0 * 0.75 - transfer) / 0.75)
        self.assertAlmostEqual(updated[1][0], transfer / 0.25)
        self.assertAlmostEqual(residual, 0.0, places=15)

    def test_lateral_diffusion_scales_flux_for_tiny_partition_footprint(self) -> None:
        updated, residual = topounit_lateral_diffusion(
            concentration=[[1.0], [0.0]],
            storage_weight=[[0.5], [1.0e-8]],
            layer_top_elevation=[[1.0], [1.0]],
            layer_bottom_elevation=[[0.0], [0.0]],
            effective_diffusivity=[[4.5e-8], [4.5e-8]],
            edges=[(0, 1, 1.0)],
            dt=3600.0,
        )
        self.assertLessEqual(updated[1][0], updated[0][0])
        self.assertGreater(updated[1][0], 0.0)
        self.assertAlmostEqual(residual, 0.0, places=15)

    def test_multitopounit_lateral_diffusion_is_conservative_and_donor_limited(self) -> None:
        concentration = [[1.0, 0.5], [0.0, 0.0], [0.0, 0.0]]
        weights = [[0.05, 0.05], [0.3, 0.3], [0.65, 0.65]]
        updated, residual = topounit_lateral_diffusion(
            concentration=concentration,
            storage_weight=weights,
            layer_top_elevation=[[1.0, 0.5], [0.9, 0.4], [0.8, 0.3]],
            layer_bottom_elevation=[[0.5, 0.0], [0.4, -0.1], [0.3, -0.2]],
            effective_diffusivity=[[1.0, 1.0]] * 3,
            edges=[(0, 1, 1.0), (0, 2, 2.0)],
            dt=86400.0,
        )
        before = sum(
            concentration[node][layer] * weights[node][layer]
            for node in range(3)
            for layer in range(2)
        )
        after = sum(
            updated[node][layer] * weights[node][layer]
            for node in range(3)
            for layer in range(2)
        )
        self.assertAlmostEqual(after, before, places=14)
        self.assertAlmostEqual(residual, 0.0, places=14)
        self.assertTrue(all(value >= 0.0 for profile in updated for value in profile))
        self.assertTrue(any(value > 0.0 for profile in updated[1:] for value in profile))

    def test_lateral_diffusion_respects_absolute_vertical_overlap(self) -> None:
        updated, residual = topounit_lateral_diffusion(
            concentration=[[2.0], [0.0]],
            storage_weight=[[0.5], [0.5]],
            layer_top_elevation=[[2.0], [1.0]],
            layer_bottom_elevation=[[1.1], [0.0]],
            effective_diffusivity=[[1.0], [1.0]],
            edges=[(0, 1, 1.0)],
            dt=1.0,
        )
        self.assertEqual(updated, [[2.0], [0.0]])
        self.assertEqual(residual, 0.0)

    def test_transport_mode_is_namelist_controlled_and_defaults_to_fickian(self) -> None:
        option = "use_elm_microbe_methane_transport"
        self.assertIn(
            f'entry id="{option}" type="logical" category="bgc"',
            self.namelist_definition,
        )
        entry = self.namelist_definition.split(f'<entry id="{option}"', 1)[1].split(
            "</entry>", 1
        )[0]
        self.assertIn('value=".false."', entry)
        self.assertIn(f"logical, public :: {option} = .false.", self.varctl_source)
        self.assertIn(f"namelist /elm_inparm/ &\n         use_nofire", self.control_source)
        self.assertIn(option, self.control_source)
        self.assertIn(
            f"{option}=.true. requires use_microbe_methane=.true.",
            self.build_namelist,
        )
        self.assertIn(f"if (.not. {option}) then", self.state_source)
        self.assertIn(f"if ({option}) then", self.state_source)

    def test_clm_humhol_saturation_parity_is_explicit_and_off_by_default(self) -> None:
        option = "use_clm_microbe_humhol_saturation"
        self.assertIn(
            f'entry id="{option}" type="logical" category="bgc"',
            self.namelist_definition,
        )
        entry = self.namelist_definition.split(f'<entry id="{option}"', 1)[1].split(
            "</entry>", 1
        )[0]
        self.assertIn('value=".false."', entry)
        self.assertIn(f"logical, public :: {option} = .false.", self.varctl_source)
        self.assertIn(f"{option}=.true. requires use_microbe_methane=.true.", self.build_namelist)
        self.assertIn(f"{option}=.true. requires use_humhol=.true.", self.build_namelist)
        self.assertIn("if (use_clm_microbe_humhol_saturation) fraction = 0.99_r8", self.state_source)

    def test_combined_four_gas_interface_jointly_limits_and_closes(self) -> None:
        dt = 86400.0
        concentration = [
            [1.0e-8, 0.1, 0.2, 0.3],
            [0.0, 0.1, 0.0, 0.0],
        ]
        _interfaces, _aerenchyma, _ebullition, tendency, surface_flux = gas_transport(
            concentration=concentration,
            layer_thickness=[0.1, 0.2],
            effective_diffusivity=[[1.0] * 4, [1.0] * 4],
            transport_capacity=[[1.0] * 4, [1.0] * 4],
            surface_equilibrium_concentration=[0.0, 0.2, 0.0, 0.0],
            surface_conductance=[1.0, 0.1, 0.0, 0.0],
            aerenchyma_equilibrium_concentration=[[0.0, 0.2, 0.0, 0.0]] * 2,
            aerenchyma_exchange_rate=[[1.0] * 4, [1.0] * 4],
            aerenchyma_minimum_emission_concentration=[[0.0] * 4] * 2,
            aerenchyma_allows_influx=[False, True, True, False],
            ch4_ebullition_threshold=[0.0, 0.0],
            ch4_ebullition_activation=[1.0, 1.0],
            dt=dt,
        )
        for layer in range(2):
            for gas in range(4):
                self.assertGreaterEqual(
                    concentration[layer][gas] + dt * tendency[layer][gas],
                    -1.0e-15,
                )
        for gas in range(4):
            gas_tendency = [layer[gas] for layer in tendency]
            self.assertAlmostEqual(
                transport_residual(gas_tendency, [0.1, 0.2], surface_flux[gas]),
                0.0,
                places=18,
            )

    def test_joint_limiter_preserves_relative_pathway_rates(self) -> None:
        interfaces, aerenchyma, ebullition, tendency, surface_flux = gas_transport(
            concentration=[[1.0, 0.0, 0.0, 0.0]],
            layer_thickness=[1.0],
            effective_diffusivity=[[0.0] * 4],
            transport_capacity=[[1.0] * 4],
            surface_equilibrium_concentration=[0.0] * 4,
            surface_conductance=[0.0] * 4,
            aerenchyma_equilibrium_concentration=[[0.0] * 4],
            aerenchyma_exchange_rate=[[1.0, 0.0, 0.0, 0.0]],
            aerenchyma_minimum_emission_concentration=[[0.0] * 4],
            aerenchyma_allows_influx=[False, True, True, False],
            ch4_ebullition_threshold=[0.0],
            ch4_ebullition_activation=[1.0],
            dt=1.0,
        )
        self.assertAlmostEqual(-interfaces[0][0], 0.0)
        self.assertAlmostEqual(aerenchyma[0][0], 0.5)
        self.assertAlmostEqual(ebullition[0], 0.5)
        self.assertAlmostEqual(tendency[0][0], -1.0)
        self.assertAlmostEqual(surface_flux[0], 1.0)

    def test_transport_is_pure_and_not_dispatched(self) -> None:
        for routine in (
            "repartitionMicrobeMethaneScalar",
            "computeMicrobeGasTransport",
            "computeMicrobeVerticalDiffusion",
            "computeMicrobeImplicitVerticalDiffusion",
            "computeMicrobeAerenchymaTransport",
            "computeMicrobeMethaneEbullition",
            "computeMicrobeTopounitLateralDiffusion",
        ):
            self.assertIn(f"pure subroutine {routine}", self.transport_source)
            self.assertNotIn(routine, DRIVER.read_text(encoding="utf-8"))
        self.assertNotIn("use CH4Mod", self.transport_source)
        self.assertNotIn("use MicrobeGasTransportMod", DRIVER.read_text(encoding="utf-8"))

    def test_adapter_retains_optional_elm_multiphase_transport_and_split_diagnostics(self) -> None:
        for contract in (
            "gasTransportProperties",
            "transport_capacity",
            "CH4ParamsInst%scale_factor_gasdiff",
            "CH4ParamsInst%scale_factor_liqdiff",
            "MM_CH4_SURF_DIFF_UNSAT",
            "MM_CH4_SURF_DIFF_SAT",
            "MM_CH4_SURF_AERE_UNSAT",
            "MM_CH4_SURF_AERE_SAT",
            "MM_CH4_SURF_EBUL_UNSAT",
            "MM_CH4_SURF_EBUL_SAT",
        ):
            self.assertIn(contract, self.state_source)
        for contract in (
            "legacyFickianGasDiffusivity",
            "legacyFickianAcetateDiffusivity",
            "legacyFickianSurfaceConductance",
            "source_fick_d_w = 1.49e-5_r8",
            "capacity = 1._r8",
            "top_pp%regional_target_ti",
            "top_pp%lateral_dist",
            "legacyFickianLateralGasDiffusivity",
            "lateral_layer_top",
            "lateral_storage_weight",
            "MM_LATERAL_C_FLUX",
            "lateral_carbon_flux_col",
        ):
            self.assertIn(contract, self.state_source)

    def test_design_records_transport_contract_and_source_repairs(self) -> None:
        design = DESIGN.read_text(encoding="utf-8")
        self.assertIn("area-transfer repartition", design)
        self.assertIn("equal-and-opposite interface fluxes", design)
        self.assertIn("regional_target_ti", design)
        self.assertIn("multi-topounit", design)
        self.assertIn("Step 3", design)


if __name__ == "__main__":
    unittest.main()
