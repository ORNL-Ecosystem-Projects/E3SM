# Microbe and revised-methane parameter validation

**Status:** integration-test parameter set; not approved for production science  
**Prepared:** 2026-09-18  
**Scope:** ELM-Peatlands microbial decomposition (DOM, bacteria, and fungi) and
the revised methane reaction/transport backend

## Purpose and requested review

This note isolates the scientific parameter questions that remain after the
CLM-SPRUCE/CLM-Microbe equations were ported into current ELM. It is intended
for review by biogeochemistry, methane, and soil-phosphorus colleagues. The
current parameter files are sufficient for code integration, conservation,
restart, and sensitivity tests. They are **not** a calibrated or validated
parameterization for US-MOz, SPRUCE, or global applications.

The most important review questions are:

1. In the legacy executable, were concentration controls truly evaluated as
   mol m-3 and specific rates as s-1, as the active equations imply, or were
   the parameter-file numbers intended to be mM and d-1 despite those
   equations?
2. Which exact `microbepar_in`, PFT parameter NetCDF, source revision, and
   restart were used for a scientifically accepted CLM-SPRUCE simulation?
3. Are the current microbial C:P ratios defensible, and should DOM, bacterial,
   and fungal stoichiometry be fixed or flexible?
4. Is `m_drAer = 0.002` really an O2:C stoichiometric ratio, and are the very
   large executable-equivalent growth and oxidation rates intended?
5. Which reaction and transport parameters should be site independent, and
   which require calibration against CH4 profiles and fluxes?

## Reproducible baseline

The revised-methane baseline is the checked-in runtime parameter file from
CLM-Microbe commit
`9c2e0a048bb3799669d32b91e6cc76efb36d4b75`:

`inputdata/lnd/clm2/paramdata/microbepar_in`

The local CLM-SPRUCE comparison revision is
`a85800bad2c2a57abff77af36ccec319b202166a`. Its parameter file contains
additional decomposition records and cannot be interpreted safely without the
matching positional-reader audit. The current ELM manifests are:

- microbial decomposition:
  `components/elm/tools/microbe_methane/phase2/phase2_test_parameters.json`;
- revised methane:
  `components/elm/tools/microbe_methane/phase3/phase3_reference_parameters.json`.

Every enabled parameter will ultimately live as a named variable in the
standard ELM parameter NetCDF. There will be no separate production
`microbepar_in`. The two JSON files are traceable build inputs for the present
test parameter NetCDF, not alternative runtime parameter systems.

## Executable unit contract and science-unit caveat

The port follows the units in which the active legacy equations evaluate their
operands, rather than copying inconsistent source comments or declarations.
The executable divides stored carbon in g C m-3 by 12, and gas mass
concentrations by their corresponding molecular weights, before evaluating the
kinetics. The resulting numerical concentrations are mol m-3 (equivalently
mmol L-1). It also multiplies reaction tendencies by a timestep in seconds.
The execution-parity conversion contract is therefore:

| Legacy value used by active code | ELM unit | Conversion |
| --- | --- | ---: |
| concentration or half-saturation in mol m-3 | mmol m-3 | x 1,000 |
| volumetric rate in mol m-3 s-1 | mmol m-3 d-1 | x 86,400,000 |
| specific growth, death, or oxidation rate in s-1 | d-1 | x 86,400 |
| biomass floor in mol C m-3 | g C m-3 | x 12.011 |
| Celsius threshold/reference used against Kelvin soil temperature | K | + 273.15 |

This conversion is internally consistent with every active concentration term
in the inspected CLM-Microbe equations and is covered by a one-layer
raw-versus-converted reaction parity test. In particular,
`microbe_methane_k_acetate=16,000 mmol C m-3` is the correct ELM representation
of the legacy executable's raw `m_dKAce=16`; changing it to 16 mmol C m-3 would
make the ELM Monod response 1,000 times more substrate-sensitive than the
legacy executable. This executable audit does **not** prove that the legacy
parameter-file author intended the units actually executed. The values below
make that separate scientific ambiguity consequential:

| Quantity | Raw legacy number | Current ELM value | Why it needs review |
| --- | ---: | ---: | --- |
| H2-methanogen growth | 0.01 | 864 d-1 | Implausibly fast if the raw number was meant as s-1; plausible as 0.01 d-1 |
| Acetate-methanogen growth | 0.008 | 691.2 d-1 | Same ambiguity |
| Aerobic methanotroph growth | 0.008 | 691.2 d-1 | Same ambiguity |
| Aerobic acetate oxidation | 0.05 | 4,320 d-1 | Strong indication that intent and executable units may differ |
| Aerobic CH4 half-saturation | 1.0 | 1,000 mmol CH4 m-3 | Execution parity; 1 mmol m-3 remains a candidate science interpretation of the source comment |
| Aerobic O2 half-saturation | 4.0 | 4,000 mmol O2 m-3 | Execution parity; 4 mmol m-3 remains a candidate science interpretation of the source comment |
| CH4 transport/ebullition threshold | 0.05 | 50 mmol CH4 m-3 | Numerically equivalent to 50 mM; pathway semantics still matter |
| Aerobic decomposition O2:C | 0.002 | 0.002 mol mol-1 | Runtime file conflicts with a source default of 2 |

The aerobic-oxidation constants expose a concrete documentation-versus-code
contradiction in the legacy source. `microbevarcon.F90` labels
`m_dKCH4OxidCH4` and
`m_dKCH4OxidO2` as `mMol/m3`, and `microbepar_in` supplies `1.0` and
`4.0`. However, `microbeMod.F90` first divides CH4 and O2 mass concentrations
by their molecular weights, producing numerical mol m-3 (equivalently mmol
L-1), and then compares those concentrations directly with `1.0` and `4.0`.
The executed legacy arithmetic therefore behaves as though the constants are
1 and 4 mol m-3, or 1,000 and 4,000 mmol m-3. The current ELM defaults preserve
that executed behavior. Values of 1 and 4 mmol m-3 instead preserve the units
stated in the legacy declaration; they are not a strict behavioral port.

All active runtime-file half-saturation constants were checked at their use
sites. No additional factor-of-1,000 correction is needed:

| Legacy name | Raw value | ELM execution-parity value, mmol m-3 |
| --- | ---: | ---: |
| `m_dKAce` | 16 | 16,000 |
| `m_dKAceProdO2` | 0.004 | 4 |
| `m_dKH2ProdAce` | 0.00165 | 1.65 |
| `m_dKCO2ProdAce` | 0.0825 | 82.5 |
| `m_dKH2ProdCH4` | 0.0000775 | 0.0775 |
| `m_dKCO2ProdCH4` | 0.00031 | 0.31 |
| `m_dKCH4ProdAce` | 0.05 | 50 |
| `m_dKCH4OxidCH4` | 1 | 1,000 |
| `m_dKCH4OxidO2` | 4 | 4,000 |
| `m_dKAOMCH4OxidCH4` | 1.5 | 1,500 |

`m_dKCH4ProdO2`, `m_dKAerO2`, and `m_dKe` are read from the legacy positional
file but do not participate in an active reaction in the inspected source
revision; the first is referenced only by commented equations and the latter
two have no active use site. They are therefore not part of the present ELM
reaction parameter set. Adding them would be a new science-path decision, not
a unit correction.

The declared-unit interpretation is scientifically credible enough to test.
For comparison, legacy `CH4Mod` uses CH4 and O2 half-saturation constants of
about 5 and 20 mmol m-3. At year 50 in the shared-O2 US-MOz run, the current
1,000/4,000 mmol m-3 constants leave the aerobic methanotroph at its minimum
biomass and strongly suppress both Monod factors. Repeat the coupled 2-by-2
sensitivity after the shared-O2 fix, treating 1/4 mmol m-3 as a candidate
science setting and 1,000/4,000 mmol m-3 as the execution-parity baseline.
Do not silently replace the parity baseline until this sensitivity and a
one-layer comparison using an archived legacy state have been reviewed.

Reviewers should not resolve this by choosing the more plausible number alone.
The preferred evidence is an archived successful executable plus a one-layer
state/rate vector. Published parameter definitions and observational ranges
are the next-best evidence.

## Microbial decomposition parameters

### Values copied or promoted from CLM-SPRUCE

| Parameter(s) | Current value(s) | Units | Provenance and validation need |
| --- | --- | --- | --- |
| `k_dom`, `k_bacteria`, `k_fungi` | 0.007, 0.22, 0.22 | d-1 probability | Values occur in text input but their population of the legacy PFT arrays is unverified |
| `m_rf_s1m`, `m_rf_s2m`, `m_rf_s3m`, `m_rf_s4m` | 0.28, 0.46, 0.55, 0.75 | 1 | SOM microbial C-use efficiencies; text-input provenance is unverified |
| `m_batm_f`, `m_fatm_f` | 0.20, 0.20 | 1 | Biomass-turnover respiration fractions; unverified text values |
| `m_bdom_f`, `m_fdom_f` | 0.10, 0.10 | 1 | Biomass lysis to DOM; unverified text values |
| `m_bs1_f`, `m_bs2_f`, `m_bs3_f` | 0.10, 0.12, 0.18 | 1 | Bacterial residue to SOM1-SOM3; unverified text values |
| `m_fs1_f`, `m_fs2_f`, `m_fs3_f` | 0.10, 0.12, 0.18 | 1 | Fungal residue to SOM1-SOM3; unverified text values |
| `m_domb_f`, `m_domf_f` | 0.30, 0.30 | 1 | DOM uptake to bacteria/fungi; unverified text values |
| `m_doms1_f`, `m_doms2_f`, `m_doms3_f` | 0.20, 0.15, 0.05 | 1 | DOM stabilization to SOM1-SOM3; unverified text values |
| `cn_bacteria`, `cn_fungi` | 5, 15 | g C g N-1 | Test values repeat source pool literals; missing verified PFT input |
| `cn_dom` | 10 | g C g N-1 | Active CLM-SPRUCE literal |
| `bacteria_pool_cn`, `fungi_pool_cn` | 5, 15 | g C g N-1 | Initial/authoritative pool ratios; determine whether these should be distinct from allocation C:N |
| `CUEmax` | 0.8 | 1 | Active source literal |
| `microbe_cue_cn_target` | 8 | g C g N-1 | Active source literal in litter-CUE equation |
| `microbe_allocation_cn_exponent` | 0.6 | 1 | Active source literal |
| `bacteria_initial_c`, `fungi_initial_c`, `dom_initial_c` | 1e-5, 1e-5, 0 | g C m-3 soil | Cold-start source literals; restarts take precedence |
| `dom_som_diffusion_multiplier` | 10 | 1 | Active source literal; DOM currently uses the SOM vertical operator |
| `microbe_som2_q10`, `microbe_som3_q10`, `microbe_som4_q10`, `microbe_dom_q10` | 1.5, 2.0, 2.5, 1.25 | 1 | Active source literals |
| `l1dom_f`, `l2dom_f`, `l3dom_f` | 0.10, 0.08, 0.06 | 1 | Litter solubilization literals |
| `s1dom_f`, `s2dom_f`, `s3dom_f`, `s4dom_f` | 0.18, 0.14, 0.10, 0.06 | 1 | SOM solubilization literals |
| `l1s1_f`, `l2s2_f`, `l3s3_f` | 0.19, 0.21, 0.23 | 1 | Direct litter stabilization literals |
| `s1s2_f`, `s2s3_f`, `s3s4_f` | 0.14, 0.23, 0.27 | 1 | Direct SOM stabilization literals |

The SOM4 fractions for bacterial residue, fungal residue, and DOM stabilization
are closure residuals, not independent parameters. The ELM reader rejects
negative residuals.

### New phosphorus parameters

CLM-SPRUCE did not define phosphorus for DOM, bacteria, or fungi. The current
values exist only so that ELM's CNP transfers and balance checks can run:

| Parameter | Current test value | Units | Status |
| --- | ---: | --- | --- |
| `cp_bacteria` | 150 | g C g P-1 | New assumption; science approval required |
| `cp_fungi` | 450 | g C g P-1 | New assumption; science approval required |
| `cp_dom` | 300 | g C g P-1 | New assumption; science approval required |

DOM-P is dissolved **organic** phosphorus carried with the DOM pool.
`solutionp_vr` is dissolved **inorganic** phosphate. They are distinct stocks.
The generic donor/receiver stoichiometry can mineralize excess organic P to
solution P or immobilize solution P. Review is needed on:

- fixed versus flexible microbial and DOM C:P;
- whether the three ratios should vary by PFT, soil order, or nutrient status;
- whether DOM C:N:P transport should remain tied to the SOM vertical operator;
- appropriate initial microbial and DOM P stocks; and
- whether methane guild biomass should remain carbon-only or acquire explicit
  N and P limitation in a later science change.

## Revised-methane parameters

The full 65-variable list below is the current ELM test baseline. “Runtime”
means the number came from the CLM-Microbe runtime file; “literal/default”
means it was absent there and came from active source or a declaration.

### Acetate production and acetogenesis

| ELM parameter | Legacy name | Current value and unit | Source/status |
| --- | --- | --- | --- |
| `microbe_methane_mfg_biomass_min` | `MFGbiomin` | 1.2011e-14 g C m-3 | Declaration default converted from mol C m-3 |
| `microbe_methane_k_acetate` | `m_dKAce` | 16,000 mmol C m-3 | Runtime; execution-parity concentration conversion verified |
| `microbe_methane_acetate_prod_max` | `m_dAceProdACmax` | 207.36 mmol C m-3 d-1 | Runtime; volumetric-rate conversion |
| `microbe_methane_k_acetate_prod_o2` | `m_dKAceProdO2` | 4 mmol O2 m-3 | Runtime; concentration conversion |
| `microbe_methane_dom_to_acetate_q10` | `m_dACMinQ10` | 3 | Runtime |
| `microbe_methane_acetogenesis_max` | `m_dH2ProdAcemax` | 4.32 mmol C m-3 d-1 | Runtime; volumetric-rate conversion |
| `microbe_methane_k_acetogenesis_h2` | `m_dKH2ProdAce` | 1.65 mmol H2 m-3 | Runtime; concentration conversion |
| `microbe_methane_k_acetogenesis_co2` | `m_dKCO2ProdAce` | 82.5 mmol CO2 m-3 | Runtime; concentration conversion |
| `microbe_methane_acetogenesis_q10` | `m_dH2AceProdQ10` | 2 | Runtime |

### Methanogens

| ELM parameter | Legacy name | Current value and unit | Source/status |
| --- | --- | --- | --- |
| `microbe_methane_h2_methanogen_growth_rate` | `m_dGrowRH2Methanogens` | 864 d-1 | Runtime; high-priority unit review |
| `microbe_methane_h2_methanogen_death_rate` | `m_dDeadRH2Methanogens` | 86.4 d-1 | Runtime; high-priority unit review |
| `microbe_methane_h2_methanogen_yield` | `m_dYH2Methanogens` | 0.015 | Runtime |
| `microbe_methane_k_h2_methanogenesis_h2` | `m_dKH2ProdCH4` | 0.0775 mmol H2 m-3 | Runtime; concentration conversion |
| `microbe_methane_k_h2_methanogenesis_co2` | `m_dKCO2ProdCH4` | 0.31 mmol CO2 m-3 | Runtime; concentration conversion |
| `microbe_methane_h2_methanogenesis_q10` | `m_dH2CH4ProdQ10` | 2 | Runtime |
| `microbe_methane_h2_methanogenesis_co2_inhibition_scale` | active literal | 9,200 mmol CO2 m-3 | Literal; execution-parity concentration conversion verified |
| `microbe_methane_acetate_methanogen_growth_rate` | `m_dGrowRAceMethanogens` | 691.2 d-1 | Runtime; high-priority unit review |
| `microbe_methane_acetate_methanogen_death_rate` | `m_dDeadRAceMethanogens` | 172.8 d-1 | Runtime; high-priority unit review |
| `microbe_methane_acetate_methanogen_yield` | `m_dYAceMethanogens` | 0.2 | Runtime |
| `microbe_methane_k_acetoclastic_methanogenesis_acetate` | `m_dKCH4ProdAce` | 50 mmol C m-3 | Runtime; concentration conversion |
| `microbe_methane_acetoclastic_methanogenesis_q10` | `m_dCH4ProdQ10` | 2 | Runtime |
| `microbe_methane_acetoclastic_methanogenesis_ch4_yield` | `m_drCH4Prod` | 0.5 mol CH4 mol-1 acetate-C | Runtime |

### Methanotrophs and oxidation

| ELM parameter | Legacy name | Current value and unit | Source/status |
| --- | --- | --- | --- |
| `microbe_methane_aerobic_methanotroph_growth_rate` | `m_dGrowRMethanotrophs` | 691.2 d-1 | Runtime; high-priority unit review |
| `microbe_methane_aerobic_methanotroph_death_rate` | `m_dDeadRMethanotrophs` | 172.8 d-1 | Runtime; high-priority unit review |
| `microbe_methane_aerobic_methanotroph_yield` | `m_dYMethanotrophs` | 0.40 | Runtime |
| `microbe_methane_k_aerobic_oxidation_ch4` | `m_dKCH4OxidCH4` | 1,000 mmol CH4 m-3 | Runtime; execution parity verified, scientific calibration unresolved |
| `microbe_methane_k_aerobic_oxidation_o2` | `m_dKCH4OxidO2` | 4,000 mmol O2 m-3 | Runtime; execution parity verified, scientific calibration unresolved |
| `microbe_methane_aerobic_oxidation_q10` | `m_dCH4OxidQ10` | 1.2 | Runtime |
| `microbe_methane_aerobic_oxidation_o2_ch4_ratio` | `m_drCH4Oxid` | 2 mol O2 mol-1 CH4 | Runtime; chemically expected stoichiometry |
| `microbe_methane_aerobic_decomp_o2_c_ratio` | `m_drAer` | 0.002 mol O2 mol-1 C | Runtime; conflicts with source default 2 and needs resolution |
| `microbe_methane_aerobic_acetate_oxidation_rate` | active literal | 4,320 d-1 | Literal converted from 0.05 s-1; high-priority review |
| `microbe_methane_anaerobic_methanotroph_growth_rate` | `m_dGrowRAOMMethanotrophs` | 345.6 d-1 | Runtime; high-priority unit review |
| `microbe_methane_anaerobic_methanotroph_death_rate` | `m_dDeadRAOMMethanotrophs` | 172.8 d-1 | Runtime; high-priority unit review |
| `microbe_methane_anaerobic_methanotroph_yield` | `m_dYAOMMethanotrophs` | 0.15 | Runtime |
| `microbe_methane_k_anaerobic_oxidation_ch4` | `m_dKAOMCH4OxidCH4` | 1,500 mmol CH4 m-3 | Runtime; execution-parity concentration conversion verified |
| `microbe_methane_anaerobic_oxidation_q10` | `m_dAOMCH4OxidQ10` | 1.2 | Runtime |
| `microbe_methane_aom_o2_inhibition_scale` | active literal | 4,600 mmol O2 m-3 | Literal; execution-parity concentration conversion verified |

### Environmental response

| ELM parameter | Legacy name | Current value and unit | Source/status |
| --- | --- | --- | --- |
| `microbe_methane_ph_min`, `_ph_opt`, `_ph_max` | `pHmin`, `pHopt`, `pHmax` | 4, 7, 10 | Declaration defaults |
| `microbe_methane_soil_ph_fallback` | none | 7 globally; 4.5 for the initial SPRUCE test | Temporary site-level fallback while native ELM soil pH is unavailable |
| `microbe_methane_saturation_reaction_threshold` | `micfinundated` literal | 0.99 | Active literal mapped to current ELM hydrology |
| `microbe_methane_reaction_t_ref` | temperature literal | 286.65 K | Active 13.5 degrees C reference expressed in K |
| `microbe_methane_aom_t_ref` | AOM temperature literal | 286.65 K | Repairs a Celsius literal compared directly with Kelvin state |
| `microbe_methane_acetate_ph_trigger` | active literal | 5.5 | Active literal |
| `microbe_methane_acidification_coefficient` | active literal | 4.2e-12 m3 mmol-1 C | Adjusted because ELM acetate is mmol m-3 |
| `microbe_methane_acetate_feedback_half_saturation` | active literal | 100 mmol C m-3 | Concentration conversion |
| `microbe_methane_soil_water_potential_min` | active literal | -10 legacy source units | Units unresolved; current adapter uses the ELM moisture scalar |

Native ELM currently allocates a soil-pH field but does not populate it for
this case. The adapter therefore uses the separate
`microbe_methane_soil_ph_fallback` parameter. Its reference value is 7, while
the initial SPRUCE test prescribes 4.5 based on the existing site assumption;
`ph_opt` remains 7 and continues to define the response curve. This fallback
must not be treated as a spatial soil-pH product or calibration. A spatial,
depth-resolved ELM pH input is still required before calibrating pH-sensitive
rates.

### Diffusion, plant transport, ebullition, and atmosphere

| ELM parameter | Legacy name | Current value and unit | Source/status |
| --- | --- | --- | --- |
| `microbe_methane_dom_diffusivity` | `dom_diffus` | 1.8e-7 m2 s-1 | Runtime; legacy operator is not dimensionally identical to ELM finite-volume transport |
| `microbe_methane_dom_relaxation_rate` | executed `dom_diffus` neighbor relaxation | 1.8e-7 s-1 | Parity-test mapping of the value's actual role in the CLM equation; read only when the test switch is enabled |
| `microbe_methane_aqueous_dom_molecular_diffusivity` | none | 1.0e-10 m2 s-1 | Initial hypothesis for physical aqueous transport; research/calibration required |
| `microbe_methane_aqueous_acetate_molecular_diffusivity` | none | 1.0e-9 m2 s-1 | Initial hypothesis for physical aqueous transport; research/calibration required |
| `microbe_methane_aqueous_dom_mobile_fraction` | none | 1.0 | Initial equilibrium-mobile-fraction hypothesis; research/calibration required |
| `microbe_methane_aqueous_acetate_mobile_fraction` | none | 1.0 | Initial equilibrium-mobile-fraction hypothesis; research/calibration required |
| `microbe_methane_aqueous_solute_dispersivity` | none | 0.01 m | Initial longitudinal-dispersivity hypothesis; soil/site validation required |
| `microbe_methane_aqueous_solute_tortuosity_exponent` | none | 2.0 | Initial saturation/tortuosity hypothesis; soil/site validation required |
| `microbe_methane_aqueous_solute_min_liquid_fraction` | none | 1.0e-4 m3 m-3 | Numerical immobility threshold; convergence and dry-soil validation required |
| `microbe_methane_aqueous_gas_diffusion_multiplier` | `m_Fick_ad` | 2 | Runtime |
| `microbe_methane_plant_transport_coefficient` | `m_dPlantTrans` | 0.007 m s-1 | Runtime; dimensional mapping and root-area scaling need validation |
| `microbe_methane_ch4_transport_threshold` | `m_dCH4min` | 50 mmol CH4 m-3 | Runtime; retained for one-way plant emission and ebullition only |
| `microbe_methane_h2_plant_transport_threshold` | `g_dMaxH2inWater` | 0.473 mmol H2 m-3 | Declaration default; one-way plant emission threshold |
| `microbe_methane_ch4_h2_root_efold_depth` | active literal | 0.25 m | Active literal |
| `microbe_methane_o2_root_efold_depth` | active literal | 0.15 m | Active literal |
| `microbe_methane_ebullition_efold_depth` | active literal | 0.35 m | Active literal |
| `microbe_methane_transport_thaw_threshold` | transport literal | 273.05 K | Corrects -0.1 degrees C used against Kelvin state |
| `microbe_methane_plant_o2_consumption_fraction` | active literal | 0.001 | Active literal |
| `microbe_methane_plant_co2_flux_fraction` | active literal | 0.001 | Active literal |
| `microbe_methane_aqueous_diffusion_t_ref` | active literal | 298 K | Active literal |
| `microbe_methane_aqueous_diffusion_temperature_exponent` | active literal | 1.87 | Active literal |
| `microbe_methane_atmospheric_ch4_mixing_ratio` | `atmch4` | 1.7e-6 mol mol-1 | Runtime fallback; atmospheric forcing takes precedence |
| `microbe_methane_atmospheric_o2_mixing_ratio` | `atmo2` | 0.20946 mol mol-1 | Declaration fallback |
| `microbe_methane_atmospheric_co2_mixing_ratio` | `atmco2` | 3.97e-4 mol mol-1 | Declaration fallback |
| `microbe_methane_atmospheric_h2_mixing_ratio` | `atmh2` | 5.5e-7 mol mol-1 | Declaration fallback |

The revised adapter also consumes existing controls from the standard ELM
parameter file. These are not duplicated in the Microbe namespace:

| Existing ELM parameter | Current test value | Role in revised methane | Validation need |
| --- | ---: | --- | --- |
| `f_sat` | 0.95 | Select aqueous transport when water-filled porosity crosses this fraction | Test phase switching and water-table sensitivity |
| `satpow` | 2 | Porosity exponent in aqueous diffusivity | Compare saturated profiles and flux timing |
| `scale_factor_gasdiff` | 1 | Air-phase diffusivity multiplier | Constrain with upland atmospheric uptake and soil-gas profiles |
| `scale_factor_liqdiff` | 1 | Aqueous diffusivity multiplier | Constrain with wetland concentration/emission profiles |
| `organic_max` | 130 kg m-3 | Scale for the peat/mineral gas-diffusion blend | Verify transferability across organic and mineral soils |

The soil hydraulic exponent `bsw` is also used in mineral-soil gas diffusion,
but remains an ordinary spatial soil parameter. Porosity, organic matter,
liquid/ice content, and layer thickness are state or surface-data inputs.
`d_con_w`, `d_con_g`, `c_h_inv`, `kh_theta`, `kh_tbase`, the gas constant, and
water/ice densities are shared physical constants rather than calibration
parameters. Both the migrated CLM-Microbe values and these reused ELM controls
belong to the standard ELM parameter/data system.

The atmospheric boundary condition is not a calibration parameter. Surface
diffusion is bidirectional on a common gas-equivalent basis. Henry-law capacity
maps that mobile concentration to the stored bulk-soil inventory and to the
aqueous phase below the water table. CH4 and H2 plant transport are one-way emission paths
that activate above their respective thresholds. CH4 ebullition also retains
the CH4 threshold. An earlier integration-test version incorrectly substituted
the 50 mmol m-3 CH4 emission threshold for the atmospheric equilibrium in the
plant exchange equation. Starting from zero soil CH4, that formulation
imported atmospheric carbon at an enormous rate and oxidation consumed it.
That was a boundary-condition defect, not evidence for excessive production at
depth or a parameter that should be tuned.

### Selectable transport mappings

The original CLM-Microbe implementation contains explicit Fickian vertical
transport for both carbon substrates and gases. DOC and acetate use
`dom_diffus*(T/298)^1.87` in separate saturated and unsaturated arrays. CH4,
O2, CO2, and H2 use the aqueous `Fick_D_w*m_Fick_ad` coefficient in saturated
and unsaturated soil. With `HUM_HOL`, additional lateral concentration-gradient
terms exchange material between hummock and hollow columns. The source also
resets a selected near-surface layer to atmospheric Henry equilibrium each
timestep instead of representing a finite atmosphere/snow/pond/topsoil
resistance.

The revised backend now defaults to the CLM-Microbe aqueous Fickian mapping.
Acetate uses `dom_diffus*(T/298)^1.87`; CH4, O2, CO2, and H2 use the active
source's `Fick_D_w*m_Fick_ad*(T/298)^1.87` coefficients and direct stored-
concentration gradients in both area partitions. The source's `1e-3` factor on
the vertical `Fick_D_w` equation is retained as executed. Conservative
backward-Euler transport replaces its sequential layer mutation, and a finite
top-half-layer Fickian boundary replaces its timestep reset to atmospheric
Henry equilibrium.

When `use_humhol=.true.`, the default mapping also generalizes the source's
two-column lateral gas diffusion to ELM's full topounit graph. It uses
`regional_target_ti` and `lateral_dist`, actual topounit/column area weights,
and overlap between absolute-elevation layer intervals. All edges are evaluated
from the same beginning state and one limiter covers each donor's simultaneous
outgoing transfers. Each pairwise flux is scaled by the smaller participating
horizontal footprint before it is divided by topounit- and partition-weighted
storage. This is required when a saturated or unsaturated partition approaches
zero area; omitting it conserves gridcell inventory but creates an unphysical
local concentration spike in the vanishing partition. This avoids hard-coded
hummock/hollow indices and the
source's fixed 75/25 area split while preserving its factor-of-ten smaller
lateral Fickian coefficient. `MM_LATERAL_C_FLUX` records the net CH4-plus-CO2
carbon transfer into each column; its area-weighted gridcell sum must be zero.

Setting `use_elm_microbe_methane_transport=.true.` selects the alternative ELM
multiphase mapping:

- gas-phase diffusion through air-filled pores using ELM `CH4Mod` soil-
  structure physics;
- aqueous diffusion below the water table, retaining `m_Fick_ad` and its
  temperature dependence;
- Henry capacity and diffusivity transforms on one gas-equivalent mobile basis;
- a finite surface resistance; and
- a conservative backward-Euler tridiagonal diffusion solve.

DOM C/N/P normally remains on ELM's standard SOM vertical operator with the
10x DOM multiplier under both gas-transport choices. In the peat accumulation
branch the base SOM diffusion is zero, so this multiplier does not mix DOM.
The separate, default-off `use_clm_microbe_dom_relaxation` parity test applies
the executed CLM adjacent-layer timescale to DOM C/N/P after reactions. It uses
a conservative backward-Euler matrix and closed boundaries rather than the
source's sequential, non-conservative update. The gas-transport choice itself
still changes acetate and gas transport, not DOM transport or reactions.
The generalized lateral operator currently moves gases only. The source's
water-flux-driven lateral DOM/nutrient/acetate advection and its special
SPRUCE vertical remapping require a separate conservative C/N/P design. The
optional ELM multiphase gas mapping likewise retains its prior no-lateral
behavior until an appropriate multiphase conductance is validated.
It does not add a high-affinity
methanotroph pathway, suppress existing unsaturated production, or retune the
reaction parameters. Previous coupled results in this document used the
aqueous-only adapter and remain reaction/O2-coupling baselines, not valid
surface-transport validation. Repeat them with `MM_CH4_SURF_DIFF`,
`MM_CH4_SURF_AERE`, and `MM_CH4_SURF_EBUL` and their `_UNSAT`/`_SAT`
contributions before drawing conclusions about upland atmospheric uptake.

The first post-change native-ARM smoke case is
`20260918MM07_US-MOz_ICB1850CNRDCTCBC_ad_spinup_multiphase_transport_kch4_1_ko2_4_1y`.
Its Docker output is under
`/output/e3sm_run/20260918MM07_US-MOz_ICB1850CNRDCTCBC_ad_spinup_multiphase_transport_kch4_1_ko2_4_1y/run`.
The one-year cold-start run completed with a mean carbon-balance error of
`-1.26e-18`. Annualized mean CH4 production was `0.03605 g C m-2 yr-1`,
oxidation was `0.01842`, and diffusive surface exchange was a small source of
`0.01656`. The unsaturated and saturated diffusion contributions were
`0.01646` and `0.0000969 g C m-2 yr-1`, respectively, and summed to the bulk
field within single-precision history output. On the final day, diffusion was
a net sink at an annualized instantaneous rate of `-0.00790 g C m-2 yr-1`,
dominated by the unsaturated contribution. This verifies signed atmospheric
uptake and the diagnostic identities, but a cold-start year is not a flux
validation.

The mature low-`K` repeat is
`20260918MM08_US-MOz_ICB1850CNRDCTCBC_ad_spinup_multiphase_transport_kch4_1_ko2_4_50y`.
It differs from the earlier `MM06` low-`K` case by the multiphase transport
implementation and its diagnostics; forcing, parameters, and hydrology are the
same. Year-50 annual means are:

| Diagnostic | MM06 aqueous-only | MM08 multiphase |
| --- | ---: | ---: |
| gross CH4 production, g C m-2 yr-1 | 5.70512 | 5.79828 |
| gross CH4 oxidation, g C m-2 yr-1 | 5.72977 | 5.30549 |
| aerobic oxidation, g C m-2 yr-1 | 5.67891 | 5.30549 |
| anaerobic oxidation, g C m-2 yr-1 | 0.050862 | 2.23e-12 |
| net `FCH4`, g C m-2 yr-1 | +0.014063 | +0.494940 |
| diffusive CH4 exchange, g C m-2 yr-1 | not separated | +0.494940 |
| unsaturated diffusion contribution, g C m-2 yr-1 | not separated | +0.492280 |
| saturated diffusion contribution, g C m-2 yr-1 | not separated | +0.002660 |
| aerenchyma / ebullition, g C m-2 yr-1 | not separated | 0 / 0 |
| maximum unsaturated CH4, mmol m-3 soil | 2.6819 | 0.03812 |
| maximum saturated CH4, mmol m-3 soil | 4.2313 | 1.7955 |
| maximum unsaturated O2, mmol m-3 soil | 337.64 | 2,110.58 |
| revised-methane additional C, g C m-2 | 13.6305 | 2.17091 |
| DOM C / litter C / SOM C, g C m-2 | 62.478 / 265.972 / 540.354 | 61.371 / 242.566 / 537.520 |
| NEE, g C m-2 yr-1 | -114.359 | -115.592 |
| saturated fraction / ZWT, m | 0.174149 / 3.12476 | 0.174152 / 3.12471 |

Maximum absolute annual carbon-balance error over the 50-year MM08 run was
`6.47e-15 g C m-2`. The three MM08 surface pathways sum to net `FCH4` after
the documented g-to-kg conversion, and the saturated/unsaturated diffusion
terms sum to bulk diffusion within history precision.

Gross production changes by only 1.6%, confirming that this remains primarily
a transport perturbation. Gas-phase transport raises the maximum unsaturated
O2 inventory but rapidly exports unsaturated CH4: maximum unsaturated CH4 falls
by about 70-fold. That substrate loss reduces unsaturated oxidation, so the
annual column remains a small source even though the boundary can produce an
atmospheric sink when soil CH4 is low. The source is almost entirely the
unsaturated diffusive pathway; plant transport and ebullition remain inactive.
This result is physically more interpretable than the aqueous-only flux, but it
does not by itself establish that the remaining unsaturated production,
oxidation kinetics, or `0.174` saturated-area fraction are realistic for an
upland forest.

## Superseded pre-shared-stress US-MOz integration signal

The `P3S` 50-year accelerated-spinup run is retained as a pre-shared-stress
baseline. It uses the corrected plant boundary, corrected gross-production
diagnostic, and explicit saturated/unsaturated reaction diagnostics, but it
predates the common-O2 and root-respiration bridge corrections documented
below. Earlier `P3J` and `P3K` gross-production values must not be used: those
histories counted the complete acetoclastic process extent as CH4 rather than
multiplying by both the non-biomass fraction and CH4 yield. The reaction state
and carbon budget were correct; only that gross history field was wrong.

Annual year-50 means from `P3S` are:

| Diagnostic | Year-50 value | Interpretation |
| --- | ---: | --- |
| gross CH4 production | 21.8421 g C m-2 yr-1 | Correct CH4-C yield accounting |
| production, saturated contribution | 6.3012 g C m-2 yr-1 | Area-weighted; 28.85% of total production |
| production, unsaturated contribution | 15.5409 g C m-2 yr-1 | Area-weighted; 71.15% of total production |
| gross CH4 oxidation | 1.26990 g C m-2 yr-1 | Aerobic plus anaerobic oxidation |
| oxidation, saturated contribution | 0.66595 g C m-2 yr-1 | Area-weighted; 52.44% of total oxidation |
| oxidation, unsaturated contribution | 0.60396 g C m-2 yr-1 | Area-weighted; 47.56% of total oxidation |
| net `FCH4`, positive to atmosphere | +20.5565 g C m-2 yr-1 | US-MOz remains a strong modeled source |
| production - oxidation - emission | +0.0157 g C m-2 yr-1 | Small net CH4-storage tendency and annual averaging residual |
| mean revised saturated-area fraction | 0.17398 | Subgrid topographic/surface-water fraction, not water-table depth |
| mean water-table depth (`ZWT`) | 3.12765 m below surface | Confirms conventionally upland column hydrology |
| revised-methane additional C stock | 10.5876 g C m-2 | Acetate, guild biomass, and represented gases |
| carbon mass-balance error | 9.24e-16 g C m-2 | Numerical closure is near machine precision |

The saturated subarea produces 36.22 g C per square metre of saturated area
per year, versus 18.81 g C per square metre of unsaturated area. Thus the
saturated subarea is about 1.9 times as productive per unit area, but the much
larger unsaturated area supplies most grid-cell production. This is a direct
warning that the moisture-scaled unsaturated fermentation pathway, as well as
the saturated-area estimate, must be validated before tuning production-rate
parameters.

### Aerobic half-saturation sensitivity and oxygen-ordering result

A 2-by-2 50-year sensitivity changed only the aerobic-methanotroph Monod
constants:

| `K_CH4` | `K_O2` | production | oxidation | `FCH4` |
| ---: | ---: | ---: | ---: | ---: |
| 1,000 | 4,000 mmol m-3 | 21.8421 | 1.26990 | 20.5565 |
| 1 | 4,000 mmol m-3 | 21.8421 | 1.26990 | 20.5565 |
| 1,000 | 4 mmol m-3 | 21.8421 | 1.26990 | 20.5565 |
| 1 | 4 mmol m-3 | 21.8421 | 1.26990 | 20.5565 |

All reported fluxes are g C m-2 yr-1. Every numeric year-50 history variable
was identical among the four cases. The run directories contain the intended
four parameter combinations, and a separate runtime-instrumented check logged
`K_CH4=1` and `K_O2=4`; this is not a parameter staging or NetCDF-read failure.

This table records the superseded sequential-demand implementation. Its
sensitivity was inactive because of integration ordering. The adapter
first removes standard heterotrophic respiration, root respiration, and
nitrification O2 demand, then evaluates methane reactions, and only afterward
replenishes O2 through gas transport. In the low-`K` runtime check, when upper-
layer CH4 reached 0.5749 mmol m-3, pre-reaction O2 was exactly zero, so aerobic
CH4 oxidation and its biomass-growth tendency were exactly zero. Annual
history can nevertheless show about 0.3 mol O2 m-3 because it records the
post-transport state. The existing 1.27 g C m-2 yr-1 oxidation is consequently
consistent with the K-independent anaerobic pathway. The adapter now records
`MM_CH4_OXID_AER` and `MM_CH4_OXID_AOM`, together with `_UNSAT` and `_SAT`
contributions; the post-fix coupled runs below use those fields directly.

Therefore the earlier explanation that high `K_CH4` and `K_O2` alone keep the
aerobic guild at its biomass floor is incomplete. Those unit choices remain
scientifically unresolved, but they cannot be calibrated in the coupled case
until methanotrophs receive a nonzero, physically defensible share of O2. The
adapter now includes potential standard respiration, root respiration,
nitrification, aerobic acetate oxidation, and aerobic methane oxidation in one
CH4Mod-style shared stress. The stress scales microbial reactions immediately
and is lagged to ELM decomposition through `o_scalar`; the independent
nitrification inventory cap was removed. The post-fix paired sensitivity below
separates that coupling artifact from the real Monod response. Compare against
bounded subcycling if the lagged scheme is timestep-sensitive.

### Root-respiration bridge correction and post-fix K sensitivity

The first shared-stress implementation still suppressed aerobic oxidation for
an implementation reason. It read `col_cf%rr_vr`, a vertically resolved root-
respiration profile that is initialized to a fill value and populated only
inside legacy `CH4Mod`. Because revised mode deliberately bypasses that
solver, the fill value entered the potential O2 demand and reduced the raw
shared stress to approximately `1e-40`. A mature-state diagnostic showed
nonzero aerobic oxidation before the O2 limiter but exactly zero realized
aerobic oxidation, with aerobic biomass pinned at `1.2011e-14 g C m-3`.

The adapter now reconstructs root respiration from the authoritative patch
fields with the same legacy weighting,
`veg_cf%rr * rootfr_patch * veg_pp%wtcol`, over the soil-patch filter. It does
not mutate or depend on the legacy `rr_vr` work array. After this correction,
the same mature-state continuation with `K_CH4=1` and `K_O2=4 mmol m-3`
produced `9.993 g C m-2 yr-1` of aerobic oxidation, and active-layer aerobic
biomass reached `1e-4` to `1e-2 g C m-3`. Raw annual/daily O2 stress is now
finite and physical (`~1e-5` to 1 in active layers), not a fill-value artifact.

A clean paired 50-year AD-spinup retest then changed only the two Monod arrays:

| Year-50 diagnostic | 1,000/4,000 mmol m-3 | 1/4 mmol m-3 |
| --- | ---: | ---: |
| gross CH4 production, g C m-2 yr-1 | 5.55181 | 5.70512 |
| total CH4 oxidation, g C m-2 yr-1 | 0.292840 | 5.72977 |
| aerobic oxidation, g C m-2 yr-1 | 2.83e-11 | 5.67891 |
| anaerobic oxidation, g C m-2 yr-1 | 0.292840 | 0.0508623 |
| net `FCH4`, g C m-2 yr-1 | +5.45275 | +0.0140628 |
| maximum saturated aerobic biomass, g C m-3 | 1.38e-14 | 0.01168 |
| maximum saturated CH4, mmol m-3 | 64.199 | 4.2313 |
| NEE, g C m-2 yr-1 | -115.300 | -114.359 |
| DOM C, g C m-2 | 60.8745 | 62.4781 |
| SOM C, g C m-2 | 536.620 | 540.354 |
| maximum absolute C-balance error, g C m-2 | 3.99e-16 | 7.73e-16 |

The run-local parameter files differ only in
`microbe_methane_k_aerobic_oxidation_ch4` and
`microbe_methane_k_aerobic_oxidation_o2`. Hydrology is identical: both have a
year-50 saturated fraction of `0.17415` and mean `ZWT=3.1248 m`. Maximum O2 is
also effectively identical at about `338 mmol m-3`. Thus the large response is
now a real kinetic sensitivity. The execution-parity 1,000/4,000 values make
the Monod product too small for biomass growth to exceed mortality at the
functional floor, whereas the stated-unit 1/4 values activate and sustain the
aerobic guild. The low-K result consumes almost all locally produced CH4 but
still gives a small positive annual surface flux, not demonstrated atmospheric
CH4 uptake. Neither setting should be called calibrated from this test alone.

## SPRUCE three-topounit integration result

The first revised-methane SPRUCE integration is
`20260918MM11_US-SPR_ICB1850CNRDCTCBC_ad_spinup_microbe_methane_3topounit_kch4_1_ko2_4_50y`.
It is a 50-year AD-spinup run generated from the current `SPRUCE.cfg`, with
`K_CH4=1 mmol m-3`, `K_O2=4 mmol m-3`, and the run-local soil-pH fallback set
to 4.5. The surface has three topounits: 50% boardwalk/fen bare ground, 17%
bog hollow, and 33% bog hummock. Hollow and hummock both use peatland PFTs
3/5/14/18 at 36/14/25/25% (evergreen needleleaf tree, deciduous needleleaf
tree, deciduous broadleaf shrub, and moss). The timestamped `0051` files hold
the annual mean for simulation year 50.

| Year-50 diagnostic | Boardwalk/fen | Bog hollow | Bog hummock | Grid mean |
| --- | ---: | ---: | ---: | ---: |
| gross CH4 production, g C m-2 yr-1 | 0.00310 | 1.87636 | 1.00786 | 0.65313 |
| gross CH4 oxidation, g C m-2 yr-1 | 0.01363 | 1.21405 | 0.51216 | 0.38222 |
| surface CH4 flux, g C m-2 yr-1 | -0.01053 | +0.69479 | +0.44662 | +0.26023 |
| revised saturated-area fraction | 0.73348 | 0.29818 | 0.05218 | 0.43465 |
| water-table depth (`ZWT`), m | 0.03513 | 0.07415 | 0.22067 | 0.10299 |
| DOM C, g C m-2 | 0.0152 | 99.954 | 131.614 | 60.432 |
| bacteria C, g C m-2 | 0.00202 | 10.476 | 11.816 | 5.681 |
| fungi C, g C m-2 | 0.00487 | 19.438 | 22.045 | 10.582 |
| SOM C, g C m-2 | 0.384 | 1357.492 | 1251.740 | 644.040 |

Year-50 `MM_CH4_SURF_AERE` is exactly zero. This is not caused by a zero
transport coefficient: `m_dPlantTrans` is 0.007. The one-way CH4 plant pathway
also requires stored CH4 above `m_dCH4min = 0.05 mol m-3` and a nonzero root
fraction in the same thawed layer. In this case the effective root profile is
confined to the top four decomposition layers, whose annual CH4 means are only
about `2.5e-5` to `1.3e-4 mol m-3`. The larger annual CH4 concentrations,
about `0.003-0.0105 mol m-3`, occur in deeper effectively unrooted layers.
Ebullition can still operate there because it is not multiplied by root
fraction. Seasonal layer-resolved maxima should be archived before deciding
whether the legacy threshold or the ELM root mapping needs revision.

The grid-mean production-minus-oxidation residual is
`0.27091 g C m-2 yr-1`; subtracting the `0.26023` surface emission leaves
about `0.01067 g C m-2 yr-1` as the annual CH4-storage tendency and averaging
residual. Aerobic oxidation accounts for effectively all oxidation; annual
AOM is `2.1e-11 g C m-2 yr-1`. The maximum year-50 grid carbon-balance error
is `2.14e-17 g C m-2`.

The long run uses `do_budgets=.false.` only to disable the older global
monthly `CNPBudgetMod` print/abort diagnostic. With that diagnostic enabled,
the three-topounit reduction accumulated a June residual of approximately
`1.1e-5 g C m-2`, which exceeded its hard-coded relative tolerance even though
the independent column/grid ecosystem balance checks remained near machine
precision. `ColCBalanceCheck`, `GridCBalanceCheck`, and the
`CMASS_BALANCE_ERROR` history field remain active in this run. The legacy
global reducer should be audited separately before it is used as a required
multi-topounit pass criterion.

### Default-Fickian rerun and lateral-footprint correction

The first default-Fickian integration attempt,
`20260918MM15_US-SPR_ICB1850CNRDCTCBC_ad_spinup_microbe_methane_3topounit_fickian_kch4_1_ko2_4_50y`,
stopped on 14 June of year 1. A saturated layer in one topounit reached
unbounded H2, CO2, and guild-biomass concentrations. The vertical implicit
Fickian transaction remained conservative. The failure was in the newly
generalized lateral adapter: a per-interface-area flux was divided by storage
that included topounit area and saturated fraction, but the flux itself had not
been multiplied by a corresponding horizontal footprint. A receiving
partition approaching zero area could therefore acquire finite inventory at
an arbitrarily large local concentration.

The lateral transaction now scales every source-layer/target-layer flux by
the smaller participating horizontal area fraction before applying the common
donor limiter. A focused regression uses a receiving storage fraction of
`1e-8` and requires finite, non-overshooting concentrations and exact weighted
inventory closure.

The corrected clean OLMT case is
`20260918MM16_US-SPR_ICB1850CNRDCTCBC_ad_spinup_microbe_methane_3topounit_fickian_kch4_1_ko2_4_50y`.
It completed 50 years with the same three-topounit surface, pH 4.5, and
`K_CH4=1`, `K_O2=4 mmol m-3` configuration. `use_elm_microbe_methane_transport`
was explicitly false. The `0051` files contain the year-50 annual mean.

| Year-50 diagnostic | Boardwalk/fen | Bog hollow | Bog hummock | Grid mean |
| --- | ---: | ---: | ---: | ---: |
| gross CH4 production, g C m-2 yr-1 | 0.000109 | 1.28987 | 0.644232 | 0.431930 |
| gross CH4 oxidation, g C m-2 yr-1 | 0.000155 | 1.23173 | 0.600920 | 0.407775 |
| surface CH4 flux, g C m-2 yr-1 | 0.000131 | 0.052025 | 0.036383 | 0.020916 |
| revised saturated-area fraction | 0.73769 | 0.30948 | 0.04534 | 0.43642 |
| water-table depth (`ZWT`), m | 0.03507 | 0.07392 | 0.22165 | 0.10325 |
| DOM C, g C m-2 | 0.0041 | 93.134 | 122.450 | 56.243 |
| bacteria C, g C m-2 | 0.0018 | 9.372 | 11.686 | 5.450 |
| fungi C, g C m-2 | 0.0045 | 19.559 | 25.561 | 11.762 |
| SOM C, g C m-2 | 0.362 | 1441.891 | 1824.082 | 847.249 |

Year-50 lateral CH4-plus-CO2 carbon tendencies are `+0.03811`, `-0.13869`,
and `+0.01370 g C m-2 yr-1` for boardwalk/fen, hollow, and hummock. Their
0.50/0.17/0.33 area-weighted sum is zero to output precision; the grid history
value is `-8.46e-22 g C m-2 s-1`. Across all 51 history records, the largest
absolute grid `MM_LATERAL_C_FLUX` is `1.81e-21 g C m-2 s-1`, the largest
absolute `CMASS_BALANCE_ERROR` is `6.66e-16 g C m-2`, and every active-column
revised-methane history value is finite. The maximum saturated H2 concentration
is `8.98e-5 mol m-3`. Year-50 aerenchyma and ebullition are zero, so all
surface CH4 exchange is Fickian diffusion in this run.

### CLM-SPRUCE SOM initialization sensitivity

The CLM-Microbe cold initializer assigns all four methane guilds in both
partitions `1.e-15 mol C m-3`, consistent with `MFGbiomin` and ELM's original
`1.2011e-14 g C m-3` seed. Its reaction loop separately raises saturated
acetoclastic methanogens below `1.e-5 mol C m-3` before every reaction call.
The 100-year reset experiment showed that this nonconservative operation
increased year-50 grid-mean production by only 4.17%, so it was removed from
ELM rather than retained as a namelist capability.

No high-biomass initialization or floor is retained. ELM cold-starts every
methane guild at `MFGbiomin`, exactly as the CLM-Microbe cold initializer does.

The decomposition-substrate sensitivity used
`SPRUCE-finalspinup-peatland-carbon-initial.nc`. CLM column 2 (hollow) was
mapped to ELM topounits 1 (boardwalk/fen) and 2 (hollow), and CLM column 1
(hummock) to ELM topounit 3 (hummock). Only SOM1--SOM4 C and N were copied.
The corresponding ELM P pools were rescaled at fixed ELM C:P because the CLM
source case is CN-only. Litter, DOM, bacteria, fungi, vegetation, hydrology,
and methane state were held identical to the control. Both treatments branched
from the same normal-mode ELM restart, avoiding a second application of ELM's
AD exit-spinup factors. A field-by-field restart comparison confirmed that the
12 SOM C/N/P fields were the only changed variables.

The paired 100-year OLMT cases are `20260920MM23` (ELM SOM control) and
`20260920MM24` (CLM SOM). Initial grid-mean SOM increased from 122.41 to
317.99 kg C m-2, principally because the CLM hollow profile filled ELM's
nearly empty 50%-area boardwalk/fen topounit. Year-50 grid-mean CH4 production
increased from 0.3156 to 0.5862 g C m-2 yr-1 and surface emission from 0.0210
to 0.0285 g C m-2 yr-1. This does not explain CLM/ELM parity: in the shared
hollow, production changed from 0.9917 to 0.8875, and in the shared hummock
from 0.4454 to 0.4456 g C m-2 yr-1. At year 100 the corresponding control to
CLM-SOM values were 0.6085 to 0.5850 for the hollow and 0.2805 to 0.3025 for
the hummock. The experiment therefore rules out SOM initialization as the
main cause of CLM-SPRUCE's much larger methane production.

A more important structural difference is the saturated-area treatment. In
CLM-SPRUCE, the `HUM_HOL` preprocessor path assigns `finundated(c) = 0.99` and
overrides `micfinundated = 0.99` during state/rate aggregation. Both reference
columns therefore run methane chemistry as approximately 99% saturated even
though their water tables differ. ELM intentionally uses
`max(FSAT, frac_h2osfc)` instead. In `20260920MM24` at year 50,
`MM_SAT_FRACTION` is 0.313 for the hollow and 0.0526 for the hummock. The CLM
reference produces 14.18 and 13.18 g C m-2 yr-1 in the hollow and hummock,
respectively, versus 0.8875 and 0.4456 in the matched-SOM ELM case. A
parity-only 0.99 saturated-fraction experiment is required before attributing
the remaining difference to reaction parameters or unit conversions.

### Conservative DOM C/N/P relaxation diagnostic

The next source-parity experiment isolates the CLM-Microbe DOM vertical
update. The source executes `dom_diffus=1.8e-7` as an adjacent-layer rate in
s-1, not as the m2 s-1 coefficient stated in its metadata. ELM now exposes a
separate test-only parameter, `microbe_methane_dom_relaxation_rate`, and the
default-off namelist switch `use_clm_microbe_dom_relaxation`. After each
reaction update, one backward-Euler finite-volume matrix is applied to DOM C,
N, and P with closed vertical boundaries. Each element has an independent
inventory residual check. This preserves the source's equal-layer relaxation
timescale without copying its sequential, non-conservative update.

Paired one-year OLMT SPRUCE cases started from the same CLM-SOM restart and
used the 0.99 saturated-fraction parity option. `20260920MM27` enabled DOM
relaxation; `20260920MM28` disabled it while sharing MM27's newly built native
ARM executable. MM28 exactly reproduced the earlier MM26 annual output for all
24 common numeric history fields, confirming that the default-off path is
bit-for-bit unchanged for the tested outputs.

The enabled update moved DOM sharply downward. In the hollow, layer-10 DOM C
increased from 0.0417 to 3.726 g C m-3; in the hummock it increased from
0.00210 to 4.650 g C m-3. C, N, and P moved together: their layer-10 changes
and total-profile responses retained the prescribed DOM stoichiometry. The
one-year active-layer DOM inventories were 5.44% larger in the hollow and
4.92% larger in the hummock than in the switch-off run. This inventory
difference is an indirect reaction response to the redistributed substrate;
the relaxation operator itself is conservative every timestep.

Grid-mean annual DOM C increased from 43.567 to 44.834 g C m-2. Gross CH4
production increased from 2.879 to 4.479 g C m-2 yr-1, oxidation from 2.733 to
3.226 g C m-2 yr-1, and surface emission from 0.0426 to 0.7170 g C m-2 yr-1.
The enabled annual `CMASS_BALANCE_ERROR` was `1.59e-14 g C m-2`. These are
diagnostic sensitivity results, not validation of the legacy mixing rate. The
large response supports continuing with a longer parity comparison before any
decision to make this pathway scientific default behavior.

The 50-year continuation is OLMT case `20260920MM29`. It shares the MM27
native ARM executable and otherwise matches the MM26 0.99-saturation,
CLM-SOM configuration. The frozen comparison is the corrected pH-4.5
CLM-SPRUCE case
`US-SPR_I1850CLM45CN_microbe_humhol_ph45_ref50`; CLM column 2 maps to the ELM
hollow and CLM column 1 maps to the ELM hummock. CLM depth rates were integrated
over the first ten soil-layer thicknesses and converted with
`catomw=12.011 g C mol-1`.

| Year-50 diagnostic | ELM off hollow | ELM relax hollow | CLM hollow | ELM off hummock | ELM relax hummock | CLM hummock |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| DOM C, g C m-2 | 59.301 | 52.733 | 99.852 | 64.860 | 56.496 | 108.548 |
| bacteria C, g C m-2 | 15.840 | 15.668 | 15.417 | 14.165 | 14.146 | 14.240 |
| fungi C, g C m-2 | 33.831 | 32.561 | 29.028 | 30.983 | 29.701 | 26.684 |
| gross CH4 production, g C m-2 yr-1 | 6.128 | 9.694 | 14.182 | 5.520 | 9.045 | 13.185 |
| gross CH4 oxidation, g C m-2 yr-1 | 6.058 | 5.310 | 11.708 | 5.418 | 4.525 | 10.772 |
| surface CH4 flux, g C m-2 yr-1 | 0.0691 | 4.417 | 1.217 | 0.0716 | 4.525 | 1.236 |

Relaxation therefore raises ELM production by 58.2% in the hollow and 63.9%
in the hummock, closing the production ratio from approximately 42--43% to
68.4--68.6% of CLM. It does not close the full difference. It also reduces
oxidation, leaving ELM at only 45.4% and 42.0% of the CLM hollow and hummock
oxidation, respectively. Surface emission becomes about 3.6 times CLM because
the additional production is not balanced by oxidation.

The profile comparison explains part of this response. Without relaxation,
only 4.5--4.7% of the first-ten-layer DOM inventory is in layers 6--10. The
relaxation experiment raises this to 54.5% in the hollow and 56.9% in the
hummock, but CLM retains 89.9% and 87.2%, respectively. ELM layer-10 DOM C is
5.83 and 6.67 g C m-3, versus 24.00 and 24.38 in CLM. Conversely, ELM retains
much more surface DOM: 188.6 and 184.4 g C m-3 versus CLM's 38.7 and 60.5.
The source CLM restart already begins with substantially more deep DOM than
the ELM restart: approximately 4 g C m-3 in layer 10 versus near zero in ELM.
The 50-year comparison therefore combines transport behavior with unequal DOM
initial conditions; only SOM C/N/P had previously been mapped from CLM.

The gas profiles identify an additional transport/oxidation difference. ELM
saturated CH4 is about 40--42% of CLM in layer 1 but 3.8 times CLM in layer 10.
ELM acetate is only about 4--9% of CLM across the profile. Thus the remaining
production gap is not explained by bacteria or fungi stocks, which are already
close. The combination of weak deep-to-surface gas redistribution, lower
oxidation, different acetate state, and unmatched initial DOM must be separated
before altering the relaxation rate. The ELM case completed all 50 years; its
year-50 `CMASS_BALANCE_ERROR` was `-4.14e-14 g C m-2` and the largest absolute
annual value, including initialization output, was `4.07e-11 g C m-2`.

### Physical aqueous advection and diffusion

The science-oriented alternative is now implemented behind the default-off
`use_microbe_aqueous_transport` switch. It is mutually exclusive with the CLM
relaxation diagnostic. When enabled, the standard decomposition-pool operator
skips dissolved pools and the revised-methane adapter transports authoritative
DOM C/N/P and acetate with one conservative finite-volume operator. Bacteria
and fungi are not transported. The operator uses ELM layer geometry, liquid
water, and interlayer liquid flux; its implicit upwind-advection and centered-
diffusion matrix is nonnegative under the tested cases and reports an elemental
inventory residual each timestep.

The physical implementation is intentionally not execution parity with
CLM-Microbe. The seven parameters in the table above are hypotheses that need
measurements or priors. DOM C/N/P share the same coefficients so a uniform
stoichiometric profile remains uniform, while nonuniform ratios can redistribute
without creating an element. Clean infiltration has zero solute. Negative top
soil-water flow is not allowed to export solute because ELM's top flux includes
ground evaporation; true exfiltration and runoff require a surface aqueous pool.
Downward flow at the bottom of the decomposition domain is an external export.
Drainage diagnosed later in ELM's hydrologic sequence is not yet connected and
must be added explicitly before interpreting long-term leaching.

OLMT case `20260920MM30` was the first one-year coupled SPRUCE smoke test. The
uncorrected top boundary completed with `CMASS_BALANCE_ERROR=-7.08e-14 g C
m-2`; mean DOM transport residuals were `-2.47e-17 g C m-2`, `-8.95e-18 g N
m-2`, and `-2.44e-19 g P m-2`. It also diagnosed an apparent aqueous C export
of `9.67e-8 g C m-2 s-1` while bottom DOM export was zero. The corresponding
C:N:P ratio identified this as DOM loss through negative top flow, revealing
the evaporation error. The adapter now clips negative top solute flow; these
MM30 fluxes are retained only as a failed-boundary diagnostic and must not be
used scientifically. A corrected rerun and explicit drainage/runoff tests are
required next. Corrected OLMT case `20260920MM31` completed one year with zero
aqueous boundary export for this hydrologic realization, a whole-column carbon
balance error of `-1.36e-14 g C m-2`, and mean DOM C/N/P operator residuals of
`3.07e-17`, `-7.13e-18`, and `3.19e-19 g element m-2`, respectively. Nonzero
opposing advective and diffusive interface fluxes confirm that the internal
operator was active; their near cancellation also makes dispersivity and
porewater-concentration scaling priority sensitivity tests.

The first 50-year extension, OLMT case `20260920MM32`, stopped at
`0003-05-23_22:00:00`. The implicit aqueous solve remained nonnegative, but its
C and N inventory residuals (`1.72e-10 g C m-2` and `1.72e-11 g N m-2`) exceeded
the generic `1e-12` reaction-state tolerance. The identical C:N:P scaling of
the residual identified conditioned-solve roundoff rather than elemental mass
creation. A separate `1e-10` relative aqueous-inventory tolerance now handles
this numerical check; the `1e-12` chemistry and nonnegative-state tolerance is
unchanged, and no transport coefficient or flux was altered.

OLMT case `20260920MM33` rebuilt with that check and completed 50 years. The
largest absolute annual DOM transport residuals were `2.25e-14 g C m-2`,
`2.24e-15 g N m-2`, and `7.47e-17 g P m-2`. Excluding the initialization
record, the largest absolute annual whole-column carbon-balance error was
`1.98e-13 g C m-2`. Diagnosed top and bottom aqueous exports were zero in every
annual record, as expected until explicit drainage and runoff coupling is
implemented. At the final timestamp, DOM, bacteria, and fungi contained
`35.53`, `12.95`, and `27.74 g C m-2`; the final mean methane production,
oxidation, and surface flux were `1.79e-7`, `1.00e-7`, and `7.93e-11` in their
history-variable flux units. Advective and diffusive DOM interface fluxes
remained large and opposing, reinforcing the need to test the assumed mobile
fraction, dispersivity, and porewater concentration mapping scientifically.

A paired 50-year sensitivity, OLMT case `20260920MM34`, removed both special
initial-condition assumptions: it used the native ELM SOM restart
`elm_normal.elm.r.0002-01-01-00000.nc` and set
`use_clm_microbe_humhol_saturation=.false.`. All other MM33 configuration,
forcing, parameters, and aqueous transport choices were retained. The native
microbial saturated fraction is `max(FSAT, frac_h2osfc)`, not `FSAT` alone; its
50-year mean was `0.398` (final `0.442`) versus the fixed `0.99` in MM33, while
the two cases had essentially identical hydrologic `FSAT` and water-table
depth.

Relative to MM33, MM34 reduced mean methane production by 25.4% (`4.72` versus
`6.33 g C m-2 yr-1`) but reduced oxidation by 63.1% (`1.23` versus `3.32 g C
m-2 yr-1`). Consequently, its mean net microbial methane production and mean
surface `FCH4` increased to `3.50` and `3.49 g C m-2 yr-1`, respectively,
versus `3.00` and `2.98 g C m-2 yr-1` in MM33. Final-year `FCH4` was `2.98`
versus `2.50 g C m-2 yr-1`. MM34 final DOM increased 29.6%, bacteria and fungi
decreased about 49--51%, and `TOTSOMC` was 61.5% lower because the paired test
also removed the CLM-mapped SOM initialization. These responses therefore
cannot separate saturation from initial-SOM effects; single-factor cases are
required for attribution. MM34 completed with maximum post-initialization
annual carbon-balance error `1.26e-13 g C m-2`, maximum DOM C/N/P transport
residuals `1.81e-14`, `1.80e-15`, and `6.07e-17 g element m-2`, and zero
diagnosed aqueous boundary export.

The clean spinup-history experiment used OLMT cases `20260920MM36`: 100 years
of AD spinup from a cold start followed automatically by 10 years of normal
final spinup. It retained native ELM microbial saturation, enabled physical
aqueous transport for both phases, used the three-topounit SPRUCE surface and
PFT fractions, and restored the SPRUCE PFT parameter overrides used to create
the earlier `elm_normal` restart. A generated-parameter-file comparison
confirmed that every numeric parameter common to MM36 and the earlier MM20 AD
case was identical. MM36 adds only the aqueous-transport parameters and the
inactive DOM-relaxation parameter. The final case read the AD
`0101-01-01` restart in normal mode and reported one AD-exit conversion of the
SOM, N, P, and deadwood pools.

At the end of 10 normal years, MM36 contained `17.21`, `3.06`, and `5.97 g C
m-2` in DOM, bacteria, and fungi; `TOTSOMC` was `54.40 kg C m-2`. Its 10-year
mean methane production, oxidation, and surface flux were `1.815`, `1.031`,
and `0.785 g C m-2 yr-1`. Compared at the same 10-normal-year timestamp, MM34
contained `68.66`, `6.39`, and `13.40 g C m-2` in DOM, bacteria, and fungi;
`TOTSOMC` was `122.39 kg C m-2`; and its corresponding methane rates were
`4.839`, `1.238`, and `3.562 g C m-2 yr-1`. Hydrology did not explain the
difference: the 10-year mean microbial saturated fractions were `0.399` and
`0.393`, respectively. MM36's maximum annual carbon-balance error after the
initialization record was `5.45e-14 g C m-2`, and its maximum annual DOM C/N/P
transport residuals were `6.11e-15`, `6.11e-16`, and `2.04e-17 g element
m-2`.

This result shows that using the physical aqueous pathway throughout spinup
substantially changes the state supplied to final spinup. It is not direct
leaching: MM36's mean AD aqueous C export was only about `0.00173 g C m-2
yr-1` and was zero by AD year 100. Redistribution changes the substrate
profile and therefore the cumulative reaction and vegetation-soil feedbacks.
Because MM20 used an older executable as well as the switch-off pathway, a
paired 100-year AD run with `use_microbe_aqueous_transport=.false.` under the
current executable is still required to attribute the difference specifically
to aqueous transport. An initial MM35 attempt omitted the five SPRUCE PFT
overrides and is retained only as a parameterization sensitivity, not as the
primary clean-chain result.

### Standalone SPRUCE non-bog gas-edge sensitivity

The three-topounit SPRUCE surface assigns 50% of the gridcell to a bare
fen/boardwalk topounit, 17% to hollow, and 33% to hummock. Because the generated
surface has no explicit `TopounitRegionalTarget`, ELM's fallback graph connects
hollow to fen and hummock to hollow. Lateral groundwater and revised-methane
gas exchange therefore include the bare fen, whereas DOM C/N/P and acetate have
vertical transport only and litter, SOM, bacteria, and fungi are immobile.

A one-year continuation from the clean MM36 final-spinup restart quantified the
existing edge behavior. In the edge-on control (`20260920MM39`), the fen gained
`0.1703 g C m-2 yr-1` through lateral revised-methane gas transport while its
local methane production was only `0.00036 g C m-2 yr-1`. Its final DOM,
bacteria, and fungi inventories remained very small. The area-weighted lateral
carbon transfer closed at roundoff, confirming that this was redistribution
rather than a gridcell source.

The paired case `20260920MM38` used the same executable and restart but set
`use_microbe_nonbog_lateral_gas_transport=.false.`. This removed the fen edge
while retaining hollow-hummock exchange. Relative to the edge-on control, the
one-year gridcell methane surface flux rose from `0.7756` to `0.7895 g C m-2
yr-1` (`+1.80%`), production rose from `1.6424` to `1.6544 g C m-2 yr-1`
(`+0.73%`), and oxidation fell from `0.8676` to `0.8608 g C m-2 yr-1`
(`-0.78%`). Fen methane oxidation fell by `78.7%`, and hollow and hummock
surface methane fluxes rose by `2.07%` and `1.72%`, respectively. Carbon
balances and the gridcell-area sum of lateral carbon transfer remained at
roundoff. A same-executable edge-on run exactly reproduced the earlier MM37
diagnostic output for the compared numeric fields, validating the default-on
compatibility path.

The non-bog gas edge is therefore a real but modest one-year gridcell effect,
not the dominant explanation for the ELM-versus-CLM methane discrepancy. The
larger structural uncertainty is that lateral water moves among these
topounits, including into the bare fen, but dissolved DOM C/N/P and acetate do
not yet accompany it. A conservative, donor-limited lateral dissolved-solute
operator should be the next transport experiment before interpreting longer
SPRUCE methane simulations.

### Paired 100-year AD plus 10-year final-spinup backend comparison

OLMT cases `20260920MM40` and `20260920MM41` provide a controlled ELM
comparison of the Microbe/revised-methane option and legacy `CH4Mod`. Both use
the same current native-ARM executable, three-topounit SPRUCE surface, plot-07
forcing, CNP-RD-CTC configuration, PFT fractions and overrides, and 100 years
of accelerated decomposition spinup followed by 10 years of normal final
spinup. The Microbe case uses native ELM saturation, the default Fickian gas
mapping, full topounit gas graph, physical aqueous DOM/acetate transport, pH
4.5, and the candidate `K_CH4/K_O2=1/4 mmol m-3` setting. The control sets
`use_microbe_methane=.false.` and executes legacy `CH4Mod`. The initialization
history record was excluded from all 10-year means.

| Diagnostic | Microbe | CH4Mod | Microbe difference |
| --- | ---: | ---: | ---: |
| final-year litter C, g C m-2 | 9,414 | 9,502 | -0.9% |
| final-year soil C, g C m-2 | 54,400 | 64,579 | -15.8% |
| final-year total column C, g C m-2 | 64,731 | 75,745 | -14.5% |
| final-year vegetation C, g C m-2 | 761 | 1,355 | -43.8% |
| 10-year mean HR, g C m-2 yr-1 | 108.36 | 260.54 | -58.4% |
| 10-year mean GPP, g C m-2 yr-1 | 505.86 | 691.00 | -26.8% |
| 10-year mean NPP, g C m-2 yr-1 | 116.80 | 257.16 | -54.6% |
| 10-year mean CH4 production, g C m-2 yr-1 | 1.815 | 5.543 | -67.3% |
| 10-year mean CH4 oxidation, g C m-2 yr-1 | 1.031 | 2.858 | -63.9% |
| 10-year mean FCH4, g C m-2 yr-1 | 0.785 | 2.678 | -70.7% |
| final-year FCH4, g C m-2 yr-1 | 0.859 | 0.894 | -3.9% |

The final-year methane fluxes are much closer than their 10-year means because
legacy `CH4Mod` has large transient emission pulses in final-spinup years 2, 3,
and 8; its annual FCH4 spans `0.89--8.66 g C m-2 yr-1`. Revised methane spans
only `0.70--0.86 g C m-2 yr-1`. Thus a single final year would conceal the
large difference in transient behavior. Over the decade, Microbe soil C falls
by about `31 g C m-2`, whereas CH4Mod soil C rises by about `61 g C m-2`, so
neither endpoint should yet be interpreted as equilibrium.

The paired result reflects the complete microbial decomposition plus revised
methane option, not methane kinetics alone. In particular, the large
vegetation, HR, NPP, and soil-C differences arise from the altered DOM and
microbial cascade and its nutrient feedbacks. Methane state and carbon-flux
accounting also differ structurally between the backends, so their native
storage diagnostics are not directly interchangeable. Both cases completed
with maximum normal-spinup annual carbon-balance errors below `5.5e-14 g C
m-2`.

### Phase 2 bridge productivity attribution

OLMT case `20260920MM42` repeats the 100-year AD plus 10-year final-spinup
experiment with
`use_microbe_methane=.true., use_legacy_ch4_with_microbe=.true.`. This retains
the DOM-bacteria-fungi cascade but substitutes legacy `CH4Mod` for the revised
methane backend. The case uses the same site, forcing, surface data, PFTs, and
plant parameters as MM40/MM41. Revised aqueous-solute, saturation, and gas-
transport suboptions are off. A second annual history stream retains native
PFT dimensions for NPP, GPP, vegetation carbon, LAI, BTRAN, and plant N/P
demand and allocation.

| 10-year/final diagnostic | Microbe + legacy CH4 | Standard + legacy CH4 | Microbe + revised CH4 |
| --- | ---: | ---: | ---: |
| mean NPP, g C m-2 yr-1 | 251.26 | 257.16 | 116.80 |
| mean GPP, g C m-2 yr-1 | 688.03 | 691.00 | 505.86 |
| mean HR, g C m-2 yr-1 | 254.80 | 260.54 | 108.36 |
| final vegetation C, g C m-2 | 1,344 | 1,355 | 761 |
| final soil C, g C m-2 | 103,785 | 64,579 | 54,400 |

The bridge therefore reproduces control productivity closely: mean NPP is
2.3% lower, mean GPP is 0.4% lower, and final vegetation C is 0.8% lower than
the standard-cascade control. At the final January restart, area-weighted
evergreen and moss LAI differ from the control by -0.1% and +0.9%,
respectively. Deciduous leaf area is zero in all three January restarts because
those PFTs are dormant, not dead. All eight vegetated hollow/hummock patches
have positive annual NPP and vegetation carbon.

This result confirms that the severe productivity reduction in MM40 is not an
inevitable consequence of adding DOM, bacteria, and fungi. At this stage the
complete revised-methane configuration remained implicated because MM40 also
enabled physical aqueous DOM/acetate transport. The matched attribution test
below separates that pathway from the rest of revised methane.

The similarity in vegetation does not validate the bridge's soil spinup. Its
final soil C is about 61% above the standard-cascade control, and its final DOM,
bacteria, and fungi stocks are 91.0, 6.9, and 14.9 g C m-2. Accelerated-spinup
mapping for the microbial cascade therefore requires separate validation even
though the vegetation result is reassuring. Both phases completed normally;
the final-spinup maximum absolute annual carbon-balance error was
`3.3e-14 g C m-2`.

### Matched aqueous-transport attribution

OLMT cases `20260920MM43` and `20260920MM44` repeat the 100-year AD plus
10-year final-spinup experiment with revised methane active in both cases. They
share the same MM42 native-ARM executable and differ only in
`use_microbe_aqueous_transport`: MM43 is off and MM44 is on. Separate native-
PFT and native-column streams retain plant N/P demand and allocation, DOM C/N/P,
mineral N/P, nutrient uptake, nitrification/denitrification, and O2 profiles.
The transport-on case matches MM40 exactly across all 404 common numeric
history file-variable pairs, confirming that the expanded diagnostics did not
change the solution.

| 10-year/final diagnostic | Transport off | Transport on | On-minus-off response |
| --- | ---: | ---: | ---: |
| mean NPP, g C m-2 yr-1 | 234.70 | 116.80 | -50.2% |
| mean GPP, g C m-2 yr-1 | 678.84 | 505.86 | -25.5% |
| mean HR, g C m-2 yr-1 | 228.60 | 108.36 | -52.6% |
| mean plant N uptake, g N m-2 yr-1 | 4.294 | 1.845 | -57.0% |
| mean plant P uptake, g P m-2 yr-1 | 0.2163 | 0.0988 | -54.3% |
| final vegetation C, g C m-2 | 1,269 | 761 | -40.0% |
| final soil C, g C m-2 | 122,190 | 54,400 | -55.5% |
| mean DOM C, g C m-2 | 46.10 | 16.90 | -63.3% |
| mean mineral N, g N m-2 | 29.97 | 110.99 | +270.3% |
| mean CH4 production, g C m-2 yr-1 | 0.362 | 1.815 | +401.6% |
| mean CH4 surface flux, g C m-2 yr-1 | 0.0203 | 0.785 | +3,760% |

Turning transport off recovers 87.7% of the NPP gap between the transport-on
run and the Microbe-plus-legacy-CH4 bridge. It recovers 95% of the corresponding
GPP gap and 87% of the final vegetation-C gap. Physical aqueous transport is
therefore the dominant cause of MM40's vegetation response; revised methane/O2
coupling explains only the smaller residual in this experiment.

The mechanism is nutrient redistribution rather than loss of total mineral
nutrient. In the final-year hollow profile, transport reduces mean mineral N in
layers 1--5 from 3.01 to 0.55 g N m-3 while raising the layers 6--10 mean from
0.22 to 88.94 g N m-3. In the hummock the corresponding change is 1.17 to 0.27
near the surface and 0.15 to 131.46 at depth. Most of the deep accumulation is
nitrate. Upper-profile solution P also declines: hollow layers 1--5 fall from
2.13 to 1.04 g P m-3, and hummock layers 1--5 from 1.14 to 0.26. Thus the high
whole-column mineral inventories in the transport-on run are poorly aligned
with active roots.

PFT allocation diagnostics confirm nutrient limitation. Transport-on N
allocation/demand ratios are 0.36, 0.48, 0.32, and 0.35 for evergreen tree,
deciduous tree, shrub, and moss, compared with 0.65, 0.93, 0.60, and 0.61 when
transport is off. P allocation ratios show nearly the same reductions. BTRAN
is nearly unchanged, so water stress does not explain the vegetation response.
The same redistribution moves DOM into deeper anoxic layers and raises methane
production, explaining why vegetation suppression and methane enhancement
occur together.

This is a numerical attribution, not validation of the transport-on science.
The current mobile fractions are one, and the molecular diffusivity,
dispersivity, tortuosity exponent, and minimum liquid fraction are uncalibrated
hypotheses. The extreme deep nitrate accumulation and strong removal of upper-
profile DOM indicate that the present transport formulation or its parameter
values are too aggressive, or that a missing process such as sorption,
immobilization, root-accessible dissolved uptake, or coupled lateral/surface
export must be represented. Both runs completed with maximum final-spinup
absolute annual carbon-balance errors below `5.6e-14 g C m-2`.

### Experimental unsaturated deep-root N access

The MM44 final profile showed that plants already had nonzero nominal access
to layers 6--7: those layers supplied 7.3% of hollow and 10.0% of hummock
column-integrated N uptake. The problem was not a missing flux path but a fixed
HUMHOL RD demand profile. It followed prescribed roots only and did not adapt
to the 209--356 g N m-3 mineral-N concentrations produced in layers 6--7.

The reconstruction implements that coupling through the general
`use_peatland_roots` capability rather than a methane-specific switch. It
retains the prescribed PFT root profile, multiplies it by local NH4 plus NO3
only where liquid water is below porosity, normalizes over layer thickness,
and passes the redistributed demand through the existing RD competition and
plant-uptake fluxes. It creates neither N nor an independent uptake term. The
adaptive profile is restricted to columns with positive peat depth and to
vascular PFTs; uplands retain the ordinary uptake profile, while moss and
other nonvascular PFTs retain their prescribed shallow profiles. Fully
saturated layers receive no adaptive demand, and an empty eligible set falls
back to the unchanged root profile. HUMHOL enables peatland roots by default.

The one-year MM52 continuation from the MM44 restart built and completed with
the new native ARM executable. At the annual endpoint, hollow layers 6 and 7
had WFPS 0.954 and 0.956 and supplied 81.9% of integrated plant N uptake;
hummock values were WFPS 0.980 and 0.989 and 88.6% of uptake. Layers at WFPS
1.000 received essentially no uptake. This verifies the intended gate and
reveals a deliberately strong response: a binary unsaturated test allows
nearly saturated layers, and concentration weighting responds sharply to the
large existing deep-N reservoir.

Matched 50-year MM53 (switch off) and MM54 (switch on) continuations used the
same MM44 restart and MM52 executable. Full adaptive access raised mean NPP
from 111.74 to 274.33 g C m-2 yr-1, plant N uptake from 1.859 to
5.453 g N m-2 yr-1, plant P uptake from 0.0999 to 0.2720 g P m-2 yr-1, and
final vegetation C from 755 to 1,211 g C m-2. The final layer-6/7 share of
integrated plant N uptake increased from 6.9% to 25.0% in the hollow and from
9.7% to 39.6% in the hummock. In those layers the control retained hundreds
of g N m-3, whereas the adaptive case reduced mineral N to approximately
0.6--2.6 g N m-3. The mechanism therefore restores nutrient acquisition, but
the unblended concentration weighting mines the inherited deep reservoir
strongly and should not yet be treated as a calibrated root-uptake model.

The vegetation response feeds back on methane substrate. Mean methane
production increased from 1.794 to 4.227 g C m-2 yr-1 and net surface flux
from 0.766 to 2.986 g C m-2 yr-1; final DOM C increased from 16.3 to
63.4 g C m-2. Maximum absolute annual carbon-balance errors remained below
`5.9e-14 g C m-2`, and hydrology changed little. The next scientific test
should blend fixed-root and concentration-adaptive profiles, or replace the
binary saturation gate with a calibrated root-activity/moisture response,
rather than accepting the full adaptive result as the default. Exact annual
values and final profiles are archived in
`phase3-artifacts/spruce-20260918/root_n_access_20260920` outside the source
tree.

### Methanogenic half-saturation screening

MM45--MM51 branch from the MM43 transport-off `0011-01-01` restart and run
matched 10-year continuations. The reference and six one-at-a-time parameter
perturbations used the same native ARM64 executable and ran concurrently in
isolated OLMT cases. Generated parameter files were compared variable by
variable; each case differs from the reference only in its intended K value.

| Perturbation | CH4 production (g C m-2 yr-1) | Change | CH4 oxidation | FCH4 | FCH4 change | NPP |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Reference | 0.3561 | -- | 0.3359 | 0.02025 | -- | 237.75 |
| `microbe_methane_k_acetate=5000` | 0.4977 | +39.8% | 0.4753 | 0.02186 | +8.0% | 240.49 |
| `microbe_methane_k_acetate=1600` | 0.5673 | +59.3% | 0.5422 | 0.02368 | +16.9% | 241.89 |
| `microbe_methane_k_acetoclastic_methanogenesis_acetate=10` | 0.6133 | +72.2% | 0.5820 | 0.02266 | +11.9% | 237.85 |
| `microbe_methane_k_acetoclastic_methanogenesis_acetate=5` | 0.7052 | +98.0% | 0.6701 | 0.02343 | +15.7% | 237.23 |
| `microbe_methane_k_h2_methanogenesis_h2=0.04` | 0.3599 | +1.1% | 0.3386 | 0.02122 | +4.8% | 237.55 |
| `microbe_methane_k_h2_methanogenesis_h2=0.02` | 0.3620 | +1.7% | 0.3400 | 0.02184 | +7.9% | 237.89 |

Acetoclastic substrate affinity is the strongest tested control. Lowering its
K from 50 to 10 or 5 mmol C m-3 reduces the fraction of saturated
acetate-methanogen cells near the biomass floor from 90% to 37%. Lowering the
DOM-to-acetate K also increases production, but depletes the final DOM stock
by 9.6--14.2% over the continuation. The H2 K sensitivity is weak. Oxidation
rises nearly in step with gross production, so none of these changes produces
a proportional increase in net atmospheric flux. NPP changes by less than
1.8%, and maximum absolute carbon-balance error remains below
`1.1e-13 g C m-2`.

The moderate `K_acetoclastic=10 mmol C m-3` case is the leading candidate for
a full spinup sensitivity, not a validated replacement for the legacy value.
Exact configurations, annual values, restart-state diagnostics, and the
reduction script are archived in
`phase3-artifacts/spruce-20260918/sensitivity_k_20260920` outside this source
tree.

## Validation experiments and diagnostics

Generate coupled site experiments with `elm_olmt` by default. Avoid using a
CIME clone as the normal test path because OLMT's generated domain, surface,
and modified parameter files can remain run-local and clone configuration can
retain absolute references to the parent. Executable sharing is still allowed
when the fresh OLMT cases have compatible build configuration.

For Docker runs, invoke CIME case setup with `--disable-git` (set
`disable_git = True` under OLMT `[case_options]`). Case-local Git bookkeeping
is unnecessary for bind-mounted development runs and can fail because host
and container ownership differ.

OLMT writes generated domain, surface, parameter, and dynamic-landuse inputs
directly into each case's known run directory. The prefix/date and complete
case name are already part of that path. Do not reintroduce a shared
`elm-olmt/temp` staging file: concurrent setup of two cases can otherwise
corrupt `surfdata.nc` or fail with an HDF5 writer collision.

Parameter decisions should be tested in a staged hierarchy:

1. **One-layer equation parity.** Use an archived legacy state vector to
   compare every gross reaction rate before transport or ELM coupling. Vary
   one concentration and one biomass stock at a time to expose mmol/mol and
   daily/second errors.
2. **Closed-column transport.** Turn reactions off and verify diffusion,
   plant transport, and ebullition independently. Check sign, threshold,
   inventory closure, and response to atmospheric CH4.
3. **Process-isolation sensitivities.** Run acetoclastic production,
   hydrogenotrophic production, aerobic oxidation, and AOM separately. Report
   guild biomass and limiting factors, not only net `FCH4`.
4. **US-MOz upland test.** Require plausible atmospheric CH4 uptake when soil
   supply is small; compare gross production, gross oxidation, soil CH4
   profiles, and each surface pathway. Do not calibrate only to net flux.
5. **SPRUCE wetland test.** Compare water-table response, vertical CH4 and O2
   profiles, seasonal emissions, and pathway partitioning against observations
   and the archived CLM-SPRUCE run.
6. **CNP sensitivity.** Vary `cp_bacteria`, `cp_fungi`, and `cp_dom`; evaluate
   microbial stocks, DOM accumulation, solution-P demand, litter/SOM spinup,
   HR, NEE, and mass balance.
7. **Long spinup and transient tests.** Evaluate both unaccelerated pool
   magnitudes and transient fluxes. Accelerated decomposition remains useful
   numerically, but equivalence to final spinup must be demonstrated for the
   additional nonlinear microbial state.

At minimum, archive annual and layer-resolved values for `FCH4`, gross CH4
production and oxidation, pathway-resolved surface fluxes, CH4/O2/CO2/H2 and
acetate concentrations, all four guild biomasses, DOM/bacteria/fungi C-N-P,
solution P, mineral N, litter/SOM pools, HR, NEE, GPP, NPP, and balance errors.

## Evidence requested from reviewers

Please attach or cite, for each proposed change:

- the exact source/runtime parameter name and unit;
- an archived input file or publication/table supporting the value;
- whether the value is universal, PFT-specific, soil-specific, or site-calibrated;
- the valid range and expected temperature/moisture/pH context;
- whether it applies to a concentration per water volume, air volume, or bulk
  soil volume; and
- a benchmark state and expected gross rate when possible.

A production parameter set should not be frozen until the legacy rate-unit
ambiguity, `m_drAer`, the three new C:P ratios, soil pH input, and transport
geometry have explicit science-owner decisions.
