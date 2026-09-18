#!/usr/bin/env python3
"""Conservation and boundary tests for Phase 3 Step 3 transport."""

from __future__ import annotations

import unittest
from pathlib import Path

import phase3
from transport_oracle import (
    aerenchyma_transport,
    bulk_concentration,
    effective_aqueous_diffusivity,
    gas_transport,
    methane_ebullition,
    repartition,
    transport_residual,
    vertical_diffusion,
)


PHASE3_DIR = Path(__file__).resolve().parent
ELM_DIR = PHASE3_DIR.parents[2]
MANIFEST = PHASE3_DIR / "phase3_reference_parameters.json"
TRANSPORT = ELM_DIR / "src" / "biogeochem" / "MicrobeGasTransportMod.F90"
STATE = ELM_DIR / "src" / "biogeochem" / "MicrobeMethaneMod.F90"
DRIVER = ELM_DIR / "src" / "main" / "elm_driver.F90"
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
            surface_equilibrium_concentration=[0.0, 0.2, 0.0, 0.0],
            surface_conductance=[1.0, 0.1, 0.0, 0.0],
            aerenchyma_equilibrium_concentration=[[0.0, 0.2, 0.0, 0.0]] * 2,
            aerenchyma_exchange_rate=[[1.0] * 4, [1.0] * 4],
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
            surface_equilibrium_concentration=[0.0] * 4,
            surface_conductance=[2.0, 0.0, 0.0, 0.0],
            aerenchyma_equilibrium_concentration=[[0.0] * 4],
            aerenchyma_exchange_rate=[[1.0, 0.0, 0.0, 0.0]],
            ch4_ebullition_threshold=[0.0],
            ch4_ebullition_activation=[1.0],
            dt=1.0,
        )
        self.assertAlmostEqual(-interfaces[0][0], 0.5)
        self.assertAlmostEqual(aerenchyma[0][0], 0.25)
        self.assertAlmostEqual(ebullition[0], 0.25)
        self.assertAlmostEqual(tendency[0][0], -1.0)
        self.assertAlmostEqual(surface_flux[0], 1.0)

    def test_transport_is_pure_and_not_dispatched(self) -> None:
        for routine in (
            "repartitionMicrobeMethaneScalar",
            "computeMicrobeGasTransport",
            "computeMicrobeVerticalDiffusion",
            "computeMicrobeAerenchymaTransport",
            "computeMicrobeMethaneEbullition",
        ):
            self.assertIn(f"pure subroutine {routine}", self.transport_source)
            self.assertNotIn(routine, DRIVER.read_text(encoding="utf-8"))
        self.assertNotIn("use CH4Mod", self.transport_source)
        self.assertNotIn("use MicrobeGasTransportMod", DRIVER.read_text(encoding="utf-8"))

    def test_design_records_transport_contract_and_source_repairs(self) -> None:
        design = DESIGN.read_text(encoding="utf-8")
        self.assertIn("area-transfer repartition", design)
        self.assertIn("equal-and-opposite interface fluxes", design)
        self.assertIn("Step 3", design)


if __name__ == "__main__":
    unittest.main()
