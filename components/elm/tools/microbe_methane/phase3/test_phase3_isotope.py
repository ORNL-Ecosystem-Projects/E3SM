#!/usr/bin/env python3
"""Tests for passive C14 routing through revised methane reactions."""

from __future__ import annotations

import re
import unittest
from dataclasses import replace
from pathlib import Path

import phase3
from isotope_oracle import advance_carbon_tracer_layer
from reaction_oracle import Environment, State
from state_update_oracle import advance_reaction_layer


PHASE3_DIR = Path(__file__).resolve().parent
ELM_DIR = PHASE3_DIR.parents[2]
MANIFEST = PHASE3_DIR / "phase3_reference_parameters.json"
ISOTOPE = ELM_DIR / "src" / "biogeochem" / "MicrobeMethaneIsotopeMod.F90"
ADAPTER = ELM_DIR / "src" / "biogeochem" / "MicrobeMethaneMod.F90"
LEGACY_METHANE = ELM_DIR / "src" / "biogeochem" / "CH4Mod.F90"
CARBON_ISO_FLUX = ELM_DIR / "src" / "biogeochem" / "CarbonIsoFluxMod.F90"
VERTICAL_TRANSPORT = ELM_DIR / "src" / "biogeochem" / "SoilLittVertTranspMod.F90"
BETR_ECOSYSTEM = ELM_DIR / "src" / "biogeochem" / "CNEcosystemDynBetrMod.F90"
C14_DECAY = ELM_DIR / "src" / "biogeochem" / "C14DecayMod.F90"
GRIDCELL_STATE = ELM_DIR / "src" / "data_types" / "GridcellDataType.F90"
COLUMN_STATE = ELM_DIR / "src" / "data_types" / "ColumnDataType.F90"
VEGETATION_STATE = ELM_DIR / "src" / "data_types" / "VegetationDataType.F90"
DRIVER = ELM_DIR / "src" / "main" / "elm_driver.F90"
CONTROL = ELM_DIR / "src" / "main" / "controlMod.F90"
BUILD_NAMELIST = ELM_DIR / "bld" / "ELMBuildNamelist.pm"
HISTORY_FILE = ELM_DIR / "src" / "main" / "histFileMod.F90"
SUBGRID_AVERAGE = ELM_DIR / "src" / "main" / "subgridAveMod.F90"
NCDIO_PIO = ELM_DIR / "src" / "main" / "ncdio_pio.F90.in"
CPL_IMPORT = ELM_DIR / "src" / "cpl" / "lnd_import_export.F90"


class Phase3IsotopeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        document = phase3.load_manifest(MANIFEST)
        cls.parameters = {
            name: float(details[1])
            for name, details in phase3.parameter_map(document).items()
        }
        cls.source = ISOTOPE.read_text(encoding="utf-8")

    def _transaction(self, state: State, fraction: float = 0.35):
        return advance_reaction_layer(
            dom_c=24.0,
            dom_n=2.4,
            dom_p=0.24,
            mineral_n=1.0,
            mineral_p=0.1,
            saturated_fraction=fraction,
            unsaturated_state=state,
            saturated_state=state,
            unsaturated_environment=Environment(286.65, 5.0, 0.7, 0.3),
            saturated_environment=Environment(286.65, 5.0, 1.0, 0.0),
            parameters=self.parameters,
            cn_dom=10.0,
            cp_dom=100.0,
            dt=1800.0,
        )

    def test_uniform_ratio_is_preserved_and_closes(self) -> None:
        bulk = State(
            acetate_c=2.0,
            acetate_methanogen_c=0.4,
            h2_methanogen_c=0.4,
            aerobic_methanotroph_c=0.4,
            anaerobic_methanotroph_c=0.4,
            conc_ch4=0.02,
            conc_o2=0.03,
            conc_co2=0.04,
            conc_h2=0.01,
        )
        transaction = self._transaction(bulk)
        ratio = 1.2e-12
        tracer = State(
            acetate_c=ratio * bulk.acetate_c,
            acetate_methanogen_c=ratio * bulk.acetate_methanogen_c,
            h2_methanogen_c=ratio * bulk.h2_methanogen_c,
            aerobic_methanotroph_c=ratio * bulk.aerobic_methanotroph_c,
            anaerobic_methanotroph_c=ratio * bulk.anaerobic_methanotroph_c,
            conc_ch4=ratio * bulk.conc_ch4,
            conc_co2=ratio * bulk.conc_co2,
        )
        dom, unsat, sat, residual = advance_carbon_tracer_layer(
            24.0, ratio * 24.0, 0.35, bulk, bulk, tracer, tracer,
            transaction, self.parameters, 1800.0,
        )
        self.assertAlmostEqual(residual, 0.0, places=24)
        self.assertAlmostEqual(dom / transaction.dom_c, ratio, places=22)
        for traced, updated in ((unsat, transaction.unsaturated_state),
                                (sat, transaction.saturated_state)):
            for name in (
                "acetate_c", "acetate_methanogen_c", "h2_methanogen_c",
                "aerobic_methanotroph_c", "anaerobic_methanotroph_c",
                "conc_ch4", "conc_co2",
            ):
                if getattr(updated, name) > 1.0e-12:
                    self.assertAlmostEqual(
                        getattr(traced, name) / getattr(updated, name), ratio, places=22
                    )

    def test_old_dom_signature_reaches_acetate_and_methane(self) -> None:
        bulk = State(
            acetate_c=0.2,
            acetate_methanogen_c=0.5,
            h2_methanogen_c=0.5,
            aerobic_methanotroph_c=0.1,
            anaerobic_methanotroph_c=0.1,
            conc_ch4=1.0e-4,
            conc_o2=0.01,
            conc_co2=0.02,
            conc_h2=0.01,
        )
        transaction = self._transaction(bulk, 1.0)
        tracer = State()
        _dom, _unsat, sat, residual = advance_carbon_tracer_layer(
            24.0, 24.0e-12, 1.0, bulk, bulk, tracer, tracer,
            transaction, self.parameters, 1800.0,
        )
        self.assertGreater(sat.conc_ch4, 0.0)
        self.assertGreater(sat.conc_co2, 0.0)
        self.assertAlmostEqual(residual, 0.0, places=24)

    def test_kernel_uses_accepted_bulk_rates_not_kinetics(self) -> None:
        self.assertIn("accepted rates", self.source)
        self.assertNotIn("computeMicrobeMethanePotentialRates", self.source)
        self.assertNotIn("computeMicrobeMethaneReactionTendencies", self.source)
        self.assertIn("initial tracer ratio", self.source)

    def test_gas_tracer_uses_linear_diffusion_and_accepted_bulk_losses(self) -> None:
        self.assertIn("computeMicrobeImplicitVerticalDiffusion", self.source)
        self.assertIn("bulk_aerenchyma_flux", self.source)
        self.assertIn("bulk_ch4_ebullition_loss", self.source)
        self.assertNotIn("computeMicrobeMethaneEbullition", self.source)
        self.assertNotIn("computeMicrobeAerenchymaTransport", self.source)

    def test_generic_microbe_pools_inherit_standard_c14_cascade_and_transport(self) -> None:
        isotope_flux = CARBON_ISO_FLUX.read_text(encoding="utf-8")
        vertical_transport = VERTICAL_TRANSPORT.read_text(encoding="utf-8")
        self.assertIn("do l = 1, ndecomp_cascade_transitions", isotope_flux)
        self.assertIn("isocol_cs%decomp_cpools_vr", isotope_flux)
        self.assertIn("c14_col_cs%decomp_cpools_vr", vertical_transport)

    def test_fire_litter_fluxes_are_mapped_for_isotopes_before_state_update(self) -> None:
        isotope_flux = CARBON_ISO_FLUX.read_text(encoding="utf-8")
        update = (ELM_DIR / "src" / "biogeochem" / "CarbonStateUpdate3Mod.F90").read_text(
            encoding="utf-8"
        )
        state_fields = set(
            re.findall(r"veg_cf%(m_[a-z0-9_]+_to_litter_fire)", update)
        )
        isotope_fields = set(
            re.findall(r"isoveg_cf%(m_[a-z0-9_]+_to_litter_fire)", isotope_flux)
        )
        self.assertTrue(state_fields)
        self.assertEqual(state_fields, isotope_fields)

    def test_betr_c14_state_update_does_not_use_c13_vegetation(self) -> None:
        source = BETR_ECOSYSTEM.read_text(encoding="utf-8")
        self.assertNotIn("c14_col_cs, c13_veg_cs, c14_col_cf", source)
        self.assertIn("c14_col_cs, c14_veg_cs, c14_col_cf", source)

    def test_isotope_gridcell_state_does_not_duplicate_bulk_history_fields(self) -> None:
        source = GRIDCELL_STATE.read_text(encoding="utf-8")
        history_block = source.split(
            "These whole-gridcell diagnostics are maintained", 1
        )[1].split("set cold-start initial values", 1)[0]
        self.assertIn("if (carbon_type == 'c12') then", history_block)
        self.assertIn("TCS_MONTH_BEGIN", history_block)
        self.assertIn("CMASS_BALANCE_ERROR", history_block)

    def test_c13_fallback_uses_matching_ratio_and_c14_restart_is_required(self) -> None:
        source = COLUMN_STATE.read_text(encoding="utf-8")
        c13 = source.split("if ( carbon_type == 'c13' ) then", 1)[1].split(
            "end if ! C13", 1
        )[0]
        c14 = source.split("if ( carbon_type == 'c14' ) then", 1)[1].split(
            "end if ! C14", 1
        )[0]
        self.assertIn("decomp_cpools_vr(i,j,k) * c3_r2", c13)
        self.assertNotIn("decomp_cpools_vr(i,j,k) * c14ratio", c13)
        self.assertIn("C14 runs require a C14 restart", c14)
        self.assertNotIn("decomp_cpools_vr(i,j,k) * c14ratio", c14)

    def test_legacy_isotope_restart_initializes_transport_boundary(self) -> None:
        source = COLUMN_STATE.read_text(encoding="utf-8")
        c13 = source.split("if ( carbon_type == 'c13' ) then", 1)[1].split(
            "end if ! C13", 1
        )[0]
        c14 = source.split("if ( carbon_type == 'c14' ) then", 1)[1].split(
            "end if ! C14", 1
        )[0]
        # SoilLittVertTransp uses nlevdecomp+1 as its lower boundary, so an
        # isotope fallback must initialize the full allocated column rather
        # than leave that boundary at the restart fill value.
        self.assertIn("do j = 1, nlevdecomp_full", c13)
        self.assertIn("this%decomp_cpools_vr(i,j,k) = 0._r8", c13)

    def test_c14_restart_requires_standard_vegetation_state(self) -> None:
        source = VEGETATION_STATE.read_text(encoding="utf-8")
        c14 = source.split("if ( carbon_type == 'c14')  then", 1)[1].split(
            "endif  ! C14 block", 1
        )[0]
        self.assertIn("C14 runs require a C14 restart; missing leafc_14", c14)

    def test_c14_cold_start_clears_complete_vegetation_and_column_state(self) -> None:
        vegetation = VEGETATION_STATE.read_text(encoding="utf-8")
        column = COLUMN_STATE.read_text(encoding="utf-8")
        veg_cold = vegetation.split("set cold-start initial values", 1)[1].split(
            "end subroutine veg_cs_init", 1
        )[0]
        col_cold = column.split("set cold-start initial values", 1)[1].split(
            "end subroutine col_cs_init", 1
        )[0]
        self.assertIn("if (carbon_type == 'c14') then", veg_cold)
        self.assertIn("this%leafc              (begp:endp) = 0._r8", veg_cold)
        self.assertIn("this%totvegc            (begp:endp) = 0._r8", veg_cold)
        self.assertIn("this%cropseedc_deficit  (begp:endp) = 0._r8", veg_cold)
        self.assertIn("carbon_type == 'c14'", vegetation)
        self.assertIn("call this%SetValues(num_patch=endp-begp+1", vegetation)
        flux_setvalues = vegetation.split("subroutine veg_cf_setvalues", 1)[1].split(
            "end subroutine veg_cf_setvalues", 1
        )[0]
        for name in (
            "cpool_leaf_gr", "cpool_leaf_storage_gr", "transfer_leaf_gr",
            "cpool_froot_gr", "cpool_froot_storage_gr", "transfer_froot_gr",
        ):
            self.assertIn(f"this%{name}(i)", flux_setvalues)
        self.assertIn("if (carbon_type == 'c14') then", col_cold)
        self.assertIn("this%decomp_cpools_vr(begc:endc,:,:) = 0._r8", col_cold)
        self.assertIn("this%ctrunc_vr(begc:endc,:) = 0._r8", col_cold)

    def test_revised_methane_c14_states_have_cold_start_history_and_restart(self) -> None:
        source = ADAPTER.read_text(encoding="utf-8")
        for name in (
            "ACETATE_C", "ACET_METH_C", "H2_METH_C",
            "AER_METHANOTROPH_C", "ANAER_METHANOTROPH_C",
            "CONC_CH4", "CONC_CO2",
        ):
            for partition in ("UNSAT", "SAT"):
                field = f"C14_MM_{name}_{partition}"
                self.assertGreaterEqual(source.count(field), 2, field)
        self.assertIn("c14ratio * biomass_seed", source)
        self.assertIn("if (use_c14) then", source)
        self.assertIn("subroutine restart_c14_state", source)
        self.assertIn("C14 runs require a C14 restart", source)
        self.assertNotIn("field = c14ratio * bulk_field", source)

    def test_revised_methane_cold_start_keeps_legacy_restart_accumulators_finite(self) -> None:
        source = LEGACY_METHANE.read_text(encoding="utf-8")
        cold_start = source.split("Initialize time varying variables", 1)[1].split(
            "end subroutine InitCold", 1
        )[0]
        self.assertIn("this%tempavg_somhr_col  (c) = 0._r8", cold_start)
        self.assertIn("this%tempavg_finrw_col  (c) = 0._r8", cold_start)

    def test_history_restart_scale_strings_match_on_file_dimension(self) -> None:
        source = HISTORY_FILE.read_text(encoding="utf-8")
        field_type = source.split("type field_info", 1)[1].split(
            "end type field_info", 1
        )[0]
        for name in (
            "p2c_scale_type", "c2l_scale_type", "l2g_scale_type", "t2g_scale_type"
        ):
            self.assertRegex(
                field_type,
                rf"character\(len=hist_dim_name_length\)\s*::\s*{name}",
            )

    def test_column_history_output_aggregates_pfts_without_overrunning_buffer(self) -> None:
        history = HISTORY_FILE.read_text(encoding="utf-8")
        subgrid = SUBGRID_AVERAGE.read_text(encoding="utf-8")
        self.assertGreaterEqual(history.count("if (type1d_out == namec"), 2)
        self.assertIn(
            "call p2c(bounds, field_pft, field_column, p2c_scale_type)",
            history,
        )
        self.assertIn(
            "call p2c(bounds, num2d, field_pft, field_column, p2c_scale_type)",
            history,
        )
        self.assertGreaterEqual(
            history.count("ERROR: incompatible history bounds for"), 2
        )
        self.assertIn("module procedure p2c_1d_scaled", subgrid)

    def test_cpl_bypass_cycles_site_forcing_and_accepts_subset_streams(self) -> None:
        source = CPL_IMPORT.read_text(encoding="utf-8")
        self.assertGreaterEqual(
            source.count("modulo(modulo(yr"), 2
        )
        for dimension in (
            "popdens_nlon", "lightng_nlon", "deposition_nlon", "aerosol_nlon"
        ):
            self.assertIn(dimension, source)
        self.assertNotIn("allocate(atm2lnd_vars%lnfm_all       (192,94,2920))", source)

    def test_one_dimensional_history_files_avoid_pnetcdf_on_reopen(self) -> None:
        history = HISTORY_FILE.read_text(encoding="utf-8")
        ncdio = NCDIO_PIO.read_text(encoding="utf-8")
        self.assertIn(
            "subroutine ncd_pio_openfile(file, fname, mode, avoid_pnetcdf)", ncdio
        )
        self.assertIn("my_io_type = pio_iotype_netcdf", ncdio)
        self.assertGreaterEqual(
            history.count("avoid_pnetcdf=.not. tape(t)%dov2xy"), 3
        )

    def test_adapter_routes_c14_through_every_enabled_transport_family(self) -> None:
        source = ADAPTER.read_text(encoding="utf-8")
        driver = DRIVER.read_text(encoding="utf-8")
        self.assertIn("advanceMicrobeMethaneCarbonTracerReactionLayer", source)
        self.assertIn("advanceMicrobeMethaneCarbonTracerGasTransport", source)
        self.assertIn("advanceMicrobeAqueousTracerTransport(c14_dom_c", source)
        self.assertIn("c14_acetate_concentration", source)
        self.assertIn("c14_col_cs%decomp_cpools_vr", source)
        self.assertIn("c14_col_cs, col_cf", driver)
        self.assertIn("C14_MM_SURFACE_CH4_FLUX", source)
        self.assertIn("C14_MM_SURFACE_CO2_FLUX", source)

    def test_revised_states_use_shared_libby_decay_and_dynamic_co2_boundary(self) -> None:
        adapter = ADAPTER.read_text(encoding="utf-8")
        decay = C14_DECAY.read_text(encoding="utf-8")
        self.assertIn("function C14DecayFactor", decay)
        self.assertIn("5568._r8 * secspday", decay)
        self.assertIn("C14DecayFactor(dt, dayspyr_mod)", adapter)
        self.assertIn("cnstate_vars%rc14_atm_patch", adapter)
        self.assertIn("atmospheric_c14_ratio_col", adapter)

    def test_configuration_allows_c14_but_still_rejects_c13(self) -> None:
        control = CONTROL.read_text(encoding="utf-8")
        build_namelist = BUILD_NAMELIST.read_text(encoding="utf-8")
        validation = control.split("subroutine validate_microbe_methane_configuration", 1)[1].split(
            "end subroutine validate_microbe_methane_configuration", 1
        )[0]
        self.assertNotIn("use_c13 .or. use_c14", validation)
        self.assertIn("if (use_c13) then", validation)
        self.assertNotIn("use_c13')) ||\n      value_is_true($nl->get_value('use_c14')", build_namelist)
        self.assertIn("does not yet support C13", build_namelist)


if __name__ == "__main__":
    unittest.main()
