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

1. Which intended-unit parameter track should replace software execution
   parity for science runs? Xu et al. (2015) and the local source support raw
   concentration controls in mmol m-3 and specific rates on hourly/daily—not
   second-based—timescales, but several later processes changed operator or
   lack an unambiguous published mapping.
2. Which exact `microbepar_in`, PFT parameter NetCDF, source revision, and
   restart were used for a scientifically accepted CLM-SPRUCE simulation?
3. Are the current microbial C:P ratios defensible, and should DOM, bacterial,
   and fungal stoichiometry be fixed or flexible?
4. Is `m_drAer = 0.002` really an O2:C stoichiometric ratio, and are the very
   large executable-equivalent growth and oxidation rates intended?
5. Which reaction and transport parameters should be site independent, and
   which require calibration against CH4 profiles and fluxes?

For short reaction/transport calibration experiments only, an observed-DOM
restoring harness is available behind `use_microbe_observed_dom_calibration`.
Its named parameters are
`microbe_methane_observed_dom_relaxation_timescale`,
`microbe_methane_observed_dom_deep_concentration`,
`microbe_methane_observed_dom_surface_amplitude`, and
`microbe_methane_observed_dom_efold_depth`. These describe a smooth fit to the
2013 SPRUCE DOC profile and are not proposed as production-science parameters.
The imposed signed carbon source or sink is diagnosed separately as
`MM_DOM_PROFILE_RESTORE`.

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

### Published Xu et al. (2015) unit contract

Xu et al. (2015), the paper describing the original microbial-functional-group
CH4 kernel, provides stronger evidence about scientific intent than the later
CLM-Microbe source declarations alone:

- Appendix A states that the reaction time step is hourly and that all state
  variables are expressed as mmol per m3 of soil/water.
- Table 1 gives concentration controls in mmol m-3 (with two aerobic-oxidation
  constants reported in mmol L-1), microbial specific rates in d-1, and the
  maximum acetate-production rate in mmol m-3 h-1.
- The original application was a one-layer, sealed incubation calculation.
  Transport and ecosystem coupling added later therefore cannot be validated
  from the incubation paper alone.
- Xu et al. explicitly identify limited growth/death data, unproven anaerobic
  oxidation in the incubation, theoretical pH feedback, and absence of
  vegetation as validation limitations. Those parameters should not be
  treated as independently field-validated merely because they appear in
  Table 1.

Reference: Xu et al. (2015),
<https://doi.org/10.1002/2015JG002935>.

This evidence makes it unlikely that a single legacy-file conversion rule was
scientifically intended. The following examples distinguish three different
objects: the Xu publication value, the later CLM-SPRUCE/CLM-Microbe raw value,
and the value needed to reproduce the later executable's arithmetic in ELM.

| Process parameter | Xu et al. (2015), converted to ELM units | Later raw value | Current execution-parity ELM value | Interpretation |
| --- | ---: | ---: | ---: | --- |
| Available-C half saturation | 12 mmol m-3 | 16 | 16,000 mmol m-3 | Later raw value is close to Xu only if it is already mmol m-3 |
| Maximum acetate production | 0.84 mmol m-3 d-1 | 2.4e-6 to 2.48e-6 | 207.36 to 214.272 mmol m-3 d-1 | Treating the raw rate as mmol m-3 s-1 gives 0.207 to 0.214 mmol m-3 d-1 and removes the concentration x1000 |
| H2-methanogen growth | 0.25 d-1 | 0.01 | 864 d-1 | Nearby source comments imply h-1; raw x24 gives 0.24 d-1, almost exactly Xu |
| Acetate-methanogen growth | 0.035 d-1 | 0.007034 to 0.008 | 607.738 to 691.2 d-1 | Raw x24 gives 0.169 to 0.192 d-1, a plausible later calibration |
| Aerobic-methanotroph growth | 0.15 d-1 | 0.008 | 691.2 d-1 | Raw x24 gives 0.192 d-1, close to Xu |
| Acetoclastic CH4 half saturation | 5 mmol m-3 | 0.05 | 50 mmol m-3 | The later raw and executable-parity values bracket, but do not match, Xu |
| Aerobic CH4 half saturation | 2.5 mmol m-3 (0.0025 mmol L-1) | 1 | 1,000 mmol m-3 | The candidate 1 mmol m-3 setting is a later declared-unit hypothesis, not the Xu value |
| Aerobic O2 half saturation | 500 mmol m-3 (0.5 mmol L-1) | 4 | 4,000 mmol m-3 | The candidate 4 mmol m-3 setting is far below Xu; 4,000 is eightfold above it |
| CH4 bubbling threshold | 0.0005 mmol m-3 | 0.05 | 50 mmol m-3 | Threshold semantics and units changed substantially after Xu |
| Plant transport coefficient | 0.68 d-1 | 0.007 | 0.007 m s-1 in ELM | The later transport operator is dimensionally different and needs independent validation |

The matching CLM-SPRUCE source strengthens this conclusion. Its declarations
label the microbial growth and death parameters as d-1 and the Monod constants
as mmol m-3. Comments beside several growth/death constants nevertheless
multiply them by 24 and cite hourly ranges, indicating that the active raw
values were probably intended as h-1. In `microbeMod.F90`, `dt` is assigned directly from
`get_step_size()` in seconds and is multiplied by the unconverted reaction and
mortality rates. The gas and substrate states are also formed by dividing
g m-3 by molecular weight, which numerically produces mol m-3 (or mmol L-1),
and are compared directly with constants labeled mmol m-3. Thus the later
source contains two independent dimensional inconsistencies:

1. hour- or day-based specific rates are stepped as if they were s-1; and
2. mmol m-3 half-saturation constants are compared with mol m-3 state values.

Supply limiters keep these mistakes from necessarily causing a numerical
blow-up, but they can turn reactions into near-instantaneous, substrate-limited
events and can pin biomass at its lower bound. That is consistent with the
behavior seen in the current execution-parity tests.

The repository history identifies the transition. CLM-SPRUCE commit
`8daeedb3ee20969a4bbf30350d5e59a45a6c5a26` is titled “Update to consider
model time steps (1800 seconds).” It added `dt=get_step_size()` and multiplied
reaction, oxidation, mortality, and transport tendencies by that second-based
`dt`, while leaving the nominal guild growth/death constants unchanged. It
also retuned selected parameters by unrelated factors (for example maximum
acetate production changed from `0.0015` to `0.00005`, a factor of 30, rather
than an hourly-to-second factor of 3,600). This is the clearest point at which
the original hourly incubation formulation became a mixed-unit ecosystem
implementation. Later optimized values were therefore calibrated through
that implementation and cannot all be repaired by one algebraic conversion.

The Xu values do **not** establish a ready-to-use SPRUCE parameter set: the
later ecosystem model was restructured, recalibrated, and optimized against
field observations. They do establish that the current executable-parity set
must not be described as the likely scientific unit interpretation. We will
therefore retain three explicitly named validation tracks:

1. **execution parity** -- reproduces the later legacy arithmetic, including
   its apparent unit mistakes;
2. **later-source intended units** -- treats later raw concentration values as
   mmol m-3, volumetric rates as mmol m-3 s-1, and specific guild rates as
   h-1 where the source's x24 comments support that interpretation; and
3. **Xu-published science values** -- maps the 2015 Table 1 values through the
   current ELM equations, changing only parameters for which the process
   mapping is unambiguous.

The later-source intended-unit track is a reconstruction hypothesis, not a
claim to reproduce the calibrated Ricciuto simulation. Its purpose is to
restore physically interpretable units before recalibration.

Inhibition scales, pH-feedback scales, transport thresholds, and ebullition
thresholds must be audited individually. They should not automatically follow
either the concentration half-saturation conversion or the specific-rate
conversion.

All active runtime-file half-saturation constants were checked at their use
sites. No additional factor-of-1,000 correction is needed **for execution
parity**:

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
| `microbe_methane_k_acetate_prod_o2` | `m_dKAceProdO2` | 4 mmol O2 m-3 | Runtime; half-saturation for aerobic acetate oxidation; concentration conversion. The commented CLM-SPRUCE fermentation-inhibition expression is not enabled. |
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
| `microbe_methane_denitrification_rate_multiplier` | none | 1 reference; 0.2 provisional SPRUCE sensitivity | Revised-backend-only potential-rate multiplier; unity preserves the inherited equation, while 0.2 produced an approximately 50% reduction in realized denitrification in the short continuation test |
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

The native nitrogen module previously retained its independent pH placeholder
of 6.5 during revised-methane runs. At SPRUCE this made the Parton
nitrification pH scalar 0.920, whereas the methane kernel used the configured
fallback pH of 4.5; the corresponding nitrification scalar is 0.364. Revised
methane now passes the same fallback pH to native nitrification while all
other backends retain the inherited 6.5 behavior. The Del Grosso/Parton
denitrification potential-rate equation has no explicit pH response, so it
does not reuse the methanogen pH curve. Its pH sensitivity remains indirect,
through the nitrate supplied by nitrification.

A one-year continuation from the completed 100-year AD plus 10-year final
spinup reduced hummock/hollow gross nitrification from 1.51/1.16 to 0.95/0.89
g N m-2 yr-1 after aligning pH, but denitrification remained 1.24/1.03 versus
1.23/1.04 g N m-2 yr-1. The existing nitrate inventory and the carbon/anoxia
limits therefore buffer denitrification on this timescale. The independent
`microbe_methane_denitrification_rate_multiplier` scales potential
denitrification only for the revised backend. Its reference value of 1
preserves the inherited equation and is not a pH response. Because nitrate is
allocated competitively among denitrification, plant uptake, and microbial
immobilization, halving the potential rate did not halve the realized flux: a
multiplier of 0.5 reduced the hummock/hollow weighted flux by only 20.4%, from
1.098 to 0.874 g N m-2 yr-1. A provisional multiplier of 0.2 produced 0.552
g N m-2 yr-1, a 49.7% realized reduction in the same one-year continuation.

An 18-year continuation with the aligned pH and a unity denitrification
multiplier did not produce the expected progressive denitrification decline.
Hummock/hollow weighted nitrification was 0.911 and denitrification was 1.098
g N m-2 yr-1 in year 1, but their year 15--18 means were 1.595 and 1.477 g N
m-2 yr-1, respectively. The mineral-soil nitrate stock increased from 3.82
to a year 15--18 mean of 4.95 g N m-2. In this simulation, ammonium supply,
gross mineralization, competition, and interannual environmental variability
allow nitrification to rebound; the pre-existing nitrate reservoir is not
drawn down. Thus the pH correction changes short-term kinetics but does not,
by itself, deliver a sustained 50% reduction in denitrification.

A cold-start 100-year accelerated plus 10-year final-spinup experiment with
the aligned pH and a 0.2 multiplier did produce the expected ecosystem
feedback. Relative to the existing 100+10 reference, the final-decade
hummock/hollow denitrification flux decreased from 1.419 to 1.103 g N m-2
yr-1 (22.3%), plant mineral-N uptake increased from 3.057 to 4.131 g N m-2
yr-1 (35.1%), and NPP increased from 193.1 to 254.9 g C m-2 yr-1 (32.0%).
Mineral N more than doubled, from 5.61 to 12.10 g N m-2, and `FPG` increased
from 0.796 to 0.811. The realized denitrification reduction is much smaller
than the one-year response because greater mineralization and substrate
stocks compensate for the lower potential-rate coefficient. This is not a
strict single-parameter attribution: the earlier reference retained native
nitrification's pH 6.5 placeholder, while the treatment uses the corrected
pH 4.5. A matched pH-4.5, unity-multiplier 100+10 control is required before
assigning the productivity response specifically to denitrification.

Oxygen handling is shared but not equation-identical. The revised adapter
publishes its saturated and unsaturated O2 concentrations and accepted O2
demands through the legacy `ch4_vars` carrier. Native denitrification converts
those quantities to an anaerobic fraction with the Arah/Vinten expression;
it does not consume O2 directly. With `anoxia_wtsat=.false.`, as in the
100-year AD plus 10-year final-spinup N-deposition experiment, that calculation
uses only the unsaturated O2 state. Enabling `anoxia_wtsat` would include the
saturated state and must be tested separately because it is expected to
increase, rather than reduce, the inferred anaerobic fraction at SPRUCE.

### Diffusion, plant transport, ebullition, and atmosphere

| ELM parameter | Legacy name | Current value and unit | Source/status |
| --- | --- | --- | --- |
| `microbe_methane_dom_diffusivity` | `dom_diffus` | 1.8e-7 m2 s-1 | Runtime; legacy operator is not dimensionally identical to ELM finite-volume transport |
| `microbe_methane_dom_relaxation_rate` | executed `dom_diffus` neighbor relaxation | 1.8e-7 s-1 | Parity-test mapping of the value's actual role in the CLM equation; read only when the test switch is enabled |
| `microbe_methane_aqueous_dom_molecular_diffusivity` | none | 1.0e-10 m2 s-1 | Initial hypothesis for physical aqueous transport; research/calibration required |
| `microbe_methane_aqueous_dom_saturated_macrodispersion` | none | 0 m2 s-1 | Default zero preserves the current aqueous operator; positive values test unresolved thawed saturated-zone mixing and require profile validation |
| `microbe_methane_aqueous_acetate_molecular_diffusivity` | none | 1.0e-9 m2 s-1 | Initial hypothesis for physical aqueous transport; research/calibration required |
| `microbe_methane_aqueous_nh4_molecular_diffusivity` | none | 1.98e-9 m2 s-1 | Initial 25 C aqueous value; temperature, peat tortuosity, and profile validation required |
| `microbe_methane_aqueous_no3_molecular_diffusivity` | none | 1.90e-9 m2 s-1 | Initial 25 C aqueous value; temperature, peat tortuosity, and profile validation required |
| `microbe_methane_aqueous_nh4_partition_coefficient` | none | 0.005 m3 water kg-1 dry soil (5 L kg-1) | Initial linear-equilibrium sorption hypothesis; strongly site-, pH-, and substrate-dependent and must be calibrated |
| `microbe_methane_aqueous_mineral_n_advection_multiplier` | none | 1.0 | Experimental sensitivity only. Unity keeps mineral-N advection coupled to the modeled water flux; non-unity represents unresolved preferential mass flow or hydraulic redistribution and must not be treated as calibrated. |
| `microbe_methane_aqueous_dom_mobile_fraction` | none | 1.0 | Initial equilibrium-mobile-fraction hypothesis; research/calibration required |
| `microbe_methane_aqueous_dom_mobile_saturation_exponent` | none | 0.0 | Zero preserves constant mobility; test 1.0 applies linear liquid-saturation scaling; research/calibration required |
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

`dom_diffus` does not originate in the Xu et al. (2015) incubation model. That
study ran each topsoil, mineral-soil, and permafrost sample as an independent
one-layer microcosm, so it contains no vertical DOM transport parameter. The
parameter appears in the later 10-layer ELM-SPRUCE integration reported by
Ricciuto et al. (2021), which describes vertical DOC and acetate transport as
Fickian diffusion and gives an optimized value of `1.8e-7` over a
`1.44e-7`--`2.16e-7` sensitivity range.

The next source-parity experiment isolates the CLM-Microbe DOM vertical
update. The source executes `dom_diffus=1.8e-7` as an adjacent-layer rate in
s-1, not as the m2 s-1 coefficient stated in its metadata. For similarly sized
adjacent layers at 298 K, the value corresponds to an approximately 64-day
e-folding time; the actual coefficient also contains the layer-thickness ratio
and `(T/298)^1.87`. It therefore cannot be compared numerically with a
finite-volume coefficient in m2 s-1. ELM now exposes a
separate test-only parameter, `microbe_methane_dom_relaxation_rate`, and the
default-off namelist switch `use_clm_microbe_dom_relaxation`. After each
reaction update, one backward-Euler finite-volume matrix is applied to DOM C,
N, and P and independently to the saturated and unsaturated acetate stores,
with closed vertical boundaries. Each transported inventory has an independent
residual check. This matches the source's use of `dom_diffus` for both DOC and
acetate while preserving the equal-layer relaxation timescale without copying
its sequential, non-conservative update.

A clean 80-year AD-spinup comparison tested DOM-only relaxation against the
same relaxation applied independently to saturated and unsaturated acetate.
The years 74--80 enclosure means were nearly unchanged: methane production was
7.8470 versus 7.8479 g C m-2 yr-1, net surface flux was 7.8563 versus
7.8569 g C m-2 yr-1, and DOM inventory was 42.5066 versus 42.5844 g C m-2.
Endpoint median model/observation ratios changed from 0.174 to 0.169 for DOC,
0.203 to 0.204 for acetate, and 0.0779 to 0.0838 for CH4. Thus omission of
acetate relaxation was a real source-structure difference, but it does not
explain the large profile mismatch. Parameter-unit interpretation and the
still-evolving spinup state remain much larger suspects.

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

### Native-ELM aqueous-flux coupling regression

The reconstructed native-ELM adapter initially read `col_wf%qflx_adv` without
restoring the hydrology code that populates that field. Native ELM allocated
the array at its `1e36 mm s-1` fill value; only the BeTR path otherwise wrote
it. The conservative aqueous solver therefore interpreted the fill value as
an enormous downward water flux, transported DOM through the entire column,
and accounted for its removal as bottom-boundary export. The transport balance
check correctly closed that mathematically valid but physically meaningless
export, so it did not expose the uninitialized input.

Native soil hydrology now saves the infiltration and interlayer `qin`/`qout`
fluxes in `qflx_adv`, with interfaces below bedrock explicitly zeroed. The
revised-methane adapter also aborts if an aqueous run receives a non-finite or
fill-valued interface flux. This is required ELM coupling infrastructure, not
a scientific departure from the intended CLM-Microbe aqueous transport.

The clean 50-year three-topounit SPRUCE regression
`20260922SPRMIC50E_..._aqueous_fluxfix` used the reference parameter file
(`K_acetoclastic=50 mmol C m-3`, not the exploratory value of 10). Over years
44--50, gridcell-weighted CH4 production was `3.580 g C m-2 yr-1`, oxidation
was `1.743 g C m-2 yr-1`, net surface flux was `1.838 g C m-2 yr-1`, and DOM
stock was `65.92 g C m-2`. The largest absolute DOM transport residual was
`2.20e-13 g C m-2` per timestep. The preceding fill-value run had zero DOM by
year 11 and only `1.09e-8 g C m-2 yr-1` CH4 production, confirming that its
collapsed methane cycle was a coupling error rather than a kinetic response.
The corrected run had a closed lower boundary in the native hydrology solve,
so diagnosed bottom DOM export was zero; lateral/runoff solute export remains
a separate future development.

### Saturation-dependent DOM mobility sensitivity

The aqueous operator now accepts a layer-dependent DOM mobile fraction. The
standard parameter file contains
`microbe_methane_aqueous_dom_mobile_saturation_exponent`; the layer value is

`f_mobile = f_mobile,max * S_liq ** exponent`.

An exponent of zero is an explicit exact-preservation branch for the existing
constant-mobile-fraction formulation. Positive exponents are experimental.
The mobile fraction multiplies the porewater concentration used by both
upwind advection and diffusion/dispersion. This is distinct from the existing
saturation/tortuosity factor in molecular conductivity: the former partitions
the DOM inventory between mobile and immobile material, while the latter
represents connectivity of the aqueous pathway. With exponent one and
`theta = porosity*S_liq`, the mobile concentration becomes `DOM/porosity`
rather than `DOM/theta`, avoiding an inverse-water-content concentration
increase as soil dries. This remains an instantaneous equilibrium assumption,
not a calibrated sorption or mobile/immobile exchange model.

Three native-ARM OLMT members started from the same year-51 restart and used
the same one-year SPRUCE forcing, executable, surface data, and methane
parameters. The enclosure mean excludes fen and renormalizes to 34% hollow
and 66% hummock.

| DOM transport | Final DOM stock (g C m-2) | CH4 production (g C m-2 yr-1) | CH4 oxidation | Net FCH4 |
| --- | ---: | ---: | ---: | ---: |
| Current aqueous, exponent 0 | 103.15 | 3.172 | 2.772 | 0.392 |
| Saturation-dependent mobile fraction, exponent 1 | 99.19 | 2.659 | 2.591 | 0.147 |
| CLM adjacent-layer relaxation | 105.84 | 7.238 | 3.330 | 2.183 |

Relative to the current aqueous operator, linear saturation dependence lowered
the final DOM stock by 3.8%, gross methane production by 16.2%, and net methane
flux by 62.5%. It did more than uniformly slow transport: because it changes
the transported concentration profile, the annual net face flux reversed from
downward to upward around 0.12--0.21 m. Standard bacteria/fungi uptake rose
from 29.81 to 33.44 g C m-2 yr-1, while revised-methane fermentation fell from
19.47 to 17.71 g C m-2 yr-1.

The CLM relaxation reference mixed substantially more DOM below 1 m and more
than doubled methane production. It is therefore useful for source parity but
is not interchangeable with water-flux-driven transport. All three diagnosed
DOM budgets closed: maximum absolute layer residuals were below
`8.2e-12 g C m-3 yr-1`. Exact configurations and plot/CSV artifacts are under
`olmt_runs/20260922_spruce_microbe_50yr/dom_mobile_sensitivity` outside the
source tree. A longer continuation is needed before selecting a mobility
exponent; this one-year test establishes direction and mechanism, not an
equilibrated profile or calibration.

### Saturated-zone DOM macrodispersion sensitivity

The physical operator retains molecular diffusion and water-flux-driven
mechanical dispersion but now permits an additional default-zero DOM
macrodispersion coefficient,
`microbe_methane_aqueous_dom_saturated_macrodispersion`. This term represents
unresolved mixing by connected macropores, preferential flow, and water-table
motion rather than faster molecular diffusion. Its layer conductivity is

`theta_liq * D_sat * f_thaw * f_sat`,

where `f_sat` ramps linearly from zero at
`microbe_methane_saturation_reaction_threshold` (0.99 in the reference file)
to one at complete liquid saturation. Harmonic face conductance requires both
adjoining layers to have a positive contribution, so the added pathway cannot
cross an unsaturated or frozen intervening layer. DOM C, N, and organic P use
the same matrix. Acetate is intentionally unchanged pending evidence that its
effective mobility requires the same unresolved process.

The coefficient is added in an explicit positive-value branch. A value of
zero therefore executes the pre-existing conductivity expression without a
floating-point reassociation, preserving the default path bit for bit.

This is distinct from the CLM adjacent-layer relaxation. It remains metric,
grid-aware, and conservative, while relaxation has an equivalent diffusivity
that grows with squared layer thickness. The initial sensitivity values are
`1e-9`, `3e-9`, and `1e-8 m2 s-1`; none is a calibrated default. These tests
retain a constant DOM mobile fraction (saturation exponent zero) so the two
hypotheses are not combined.

The process hypothesis is supported by peat tracer studies rather than by a
direct calibration of `D_sat`. Ours et al. (1997) found preferential transport
through active macropores, exchange with dead pore space, and longitudinal
dispersivities ranging from centimeters to tens of centimeters in one-meter
bog-peat cores
([doi:10.1016/S0022-1694(96)03247-7](https://doi.org/10.1016/S0022-1694(96)03247-7)).
Hoag and Price (1997) examined matrix diffusion and solute retardation in peat
([manuscript](https://uwaterloo.ca/wetlands-hydrology/sites/default/files/uploads/files/1997_-_rob_s_hoag_-_theeffectsofmatrixdiffusiononsolutetransportandretretrieved-2016-07-19.pdf)),
and Rezanezhad et al. (2012) resolved mobile--immobile exchange and diffusion
into dead-end pores in breakthrough experiments
([doi:10.4141/CJSS2011-050](https://doi.org/10.4141/CJSS2011-050)). Liu et al.
(2020) related preferential transport to peat structure, anisotropy,
macropores, and saturated conductivity
([doi:10.1002/hyp.13717](https://doi.org/10.1002/hyp.13717)). Conversely,
Simhayov et al. (2018) found a conventional convection--dispersion equation
adequate for conservative chloride in one constructed-fen peat, demonstrating
that enhanced or dual-domain transport is peat dependent
([doi:10.5194/soil-4-63-2018](https://doi.org/10.5194/soil-4-63-2018)).
Ricciuto et al. (2021) remains the SPRUCE profile and flux evaluation target
([doi:10.1029/2019JG005468](https://doi.org/10.1029/2019JG005468)).

Four native-ARM OLMT members used the same year-51 SPRUCE restart, forcing,
executable, pH-4.5 parameter file, and constant mobile fraction. A fifth CLM
relaxation member provides a reference. Enclosure means exclude fen and use
34% hollow plus 66% hummock weights.

| Saturated mixing | Final DOM stock (g C m-2) | CH4 production (g C m-2 yr-1) | CH4 oxidation | Net FCH4 |
| --- | ---: | ---: | ---: | ---: |
| `D_sat=0` | 103.153 | 3.172 | 2.772 | 0.392 |
| `D_sat=1e-9 m2 s-1` | 103.149 | 3.199 | 2.786 | 0.395 |
| `D_sat=3e-9 m2 s-1` | 103.155 | 3.233 | 2.801 | 0.405 |
| `D_sat=1e-8 m2 s-1` | 103.210 | 3.302 | 2.848 | 0.413 |
| CLM relaxation | 105.839 | 7.238 | 3.330 | 2.183 |

The largest tested coefficient increased final DOM from 1.01 to
2.00 g C m-3 at 1.04 m and from 0.228 to 0.497 g C m-3 at 1.73 m. It raised
CH4 production by 4.1% and net flux by 5.4%, establishing the expected
direction without approaching relaxation's deep DOM or methane response in
one year. The zero member matched all 103 common numeric history variables
from the pre-feature control bit for bit. Every positive-coefficient DOM
budget also closed, with maximum absolute layer residual below
`8.4e-12 g C m-3 yr-1`. Longer runs are required to determine whether the
profile response accumulates and whether `D_sat` should remain fixed or be
linked to diagnosed water-table motion or hydraulic flux.

#### Water-table-overlap gate audit

The initial macrodispersion gate used layer-mean relative liquid saturation.
In the five-year mature-restart control, neither the layer intersected by the
water table nor the layer immediately below it ever reached the 0.99 reaction
threshold in the hummock or hollow. The extra pathway therefore switched on
only near 0.6 m and deeper, despite an enclosure-mean water-table depth of
approximately 0.126 m. This is a geometric inconsistency: a layer can contain
a connected saturated interval below the water table while its layer-mean
liquid saturation remains below 0.99.

The default-off `use_microbe_zwt_macrodispersion` sensitivity replaces only
that gate with the fraction of layer thickness below the connected water
table,

```text
f_zwt = clamp((z_bottom - max(z_top, max(0, ZWT))) /
              (z_bottom - z_top), 0, 1).
```

The existing liquid-volume and thawed-fraction multipliers remain in the
conductivity. Thus the option does not create water, does not bypass ice
impedance, and does not alter the conservative transport matrix. With the
option disabled, the previous 0.99 liquid-saturation threshold remains
unchanged. `MM_DOM_MACRODISP_SCALAR` and `MM_DOM_DIFF_CONDUCTIVITY` expose the
applied layer scalar and total diffusion/dispersion conductivity.

Three five-year native-ARM OLMT continuations started from the same mature
SPRUCE restart. Preferential flow was disabled. The enclosure mean excludes
fen and combines 34% hollow with 66% hummock.

| Gate and coefficient | Final DOM stock | DOC/obs median | Acetate/obs median | CH4/obs median | CH4 production | CH4 oxidation | Net FCH4 | NPP |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.99 threshold, `D_sat=1e-7` | 45.52 | 0.018 | 0.867 | 0.605 | 3.82 | 3.68 | 1.54 | 297.3 |
| ZWT overlap, `D_sat=1e-8` | 38.77 | 0.193 | 1.250 | 0.826 | 6.51 | 5.94 | 1.36 | 287.5 |
| ZWT overlap, `D_sat=1e-7` | 13.54 | 0.131 | 1.540 | 1.041 | 11.99 | 6.94 | 2.50 | 264.6 |

Stocks and NPP are `g C m-2` and `g C m-2 yr-1`, respectively; methane fluxes
are `g C m-2 yr-1`. The profile ratios are medians of modeled/observed values
matched by sampling depth and day of year.

The geometric gate corrected the intended activation. At the 0.119 m layer,
the mean scalar changed from 0.014 in the threshold control to 0.505; at
0.212 m it changed from 0.025 to 0.823. With `D_sat=1e-7`, DOC at 1.04 m rose
from effectively zero to 5.48 g C m-3 water and at 1.73 m from zero to
1.01 g C m-3 water. Those values remain far below the approximately
50--60 g C m-3 observed at those depths. The fifth-year downward diffusive
DOM-C flux was 9.53 g C m-2 yr-1 at 1.04 m and 1.17 at 1.73 m, while local
fermentation losses were 9.73 and 7.94 g C m-2 yr-1. The enhanced transport is
therefore delivering DOM to depth, but biology consumes most of it instead of
allowing the observed deep concentration to accumulate.

This experiment supports a ZWT-aware gate as the more internally consistent
geometry, but does not select `D_sat`. The `1e-7` member strongly depleted the
column DOM stock over five years, depressed NPP, overshot the acetate profile,
and still did not reproduce deep DOC. The next attribution should separate
deep fermentation and standard microbial uptake from DOM supply using the
full layer budget, then test their kinetics before increasing transport
further. A longer run is also required because none of these five-year members
is demonstrably equilibrated.

### Methane transport and storage attribution

A five-member, one-year OLMT experiment isolated the high-flux/low-porewater-
CH4 result from the observed-DOM calibration. All members started from the
same productive year-6 restart, retained DOM and acetate relaxation, disabled
aqueous transport, used identical microbial kinetics, excluded fen from the
enclosure mean, and changed only the CH4 threshold, plant coefficient, and gas
diffusion multiplier. The surface history fields separated diffusion,
aerenchyma, and ebullition. Initial and final restart inventories provided an
independent storage term.

| Member | Threshold input | Plant coefficient | Diffusion multiplier | Production | Oxidation | Diffusion | Plant | Ebullition | Net flux | Delta storage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| current | 0.05 | 0.007 | 2 | 5.732 | 0.00060 | 0.077 | 0.098 | 5.561 | 5.736 | -0.005 |
| corrected threshold | 50 | 0.007 | 2 | 5.732 | 0.00067 | 0.814 | 0 | 2.789 | 3.604 | 2.127 |
| corrected, no plant | 50 | 0 | 2 | 5.732 | 0.00067 | 0.815 | 0 | 2.789 | 3.604 | 2.127 |
| diffusion only | 1e9 | 0 | 2 | 5.732 | 0.00069 | 0.877 | 0 | 0 | 0.877 | 4.854 |
| nearly sealed | 1e9 | 0 | 1e-8 | 5.732 | 0.00069 | 3.0e-7 | 0 | 0 | 3.0e-7 | 5.731 |

Carbon terms are g C m-2 over the experiment year. Every member closed
`production - oxidation - surface pathways - storage change` to better than
`1e-7 g C m-2`, and NPP was identical (`614.7 g C m-2`). The state-to-stock
mapping and pathway diagnostics are therefore internally consistent.

The `0.05` input sends 97% of production directly to ebullition, explaining
the high surface flux and nearly unchanged, extremely low methane inventory.
The `50 mmol m-3` parity value cuts ebullition approximately in half and allows
storage to rise by `2.13 g C m-2`. Corrected-threshold members with and without
the plant coefficient are indistinguishable because the same threshold also
acts as the minimum concentration for one-way plant transport. The code should
eventually use separate parameters for the ebullition threshold and the plant-
transport minimum; their physical meanings and calibration targets differ.

The original profile analysis incorrectly treated the conserved
`MM_CONC_CH4_*` bulk-soil inventory density as a dissolved porewater
concentration. A non-prognostic diagnostic now reconstructs dissolved CH4 from
the bulk inventory, actual liquid and air volume fractions, and temperature-
dependent air--water partition capacity. It reports `MM_CH4_POREWATER`, plus
separate saturated- and unsaturated-subarea values, in `mol m-3 water`
(numerically `mmol L-1`). This conversion does not alter methane state,
transport, reactions, restart data, or the carbon budget.

The five cases were rerun with that diagnostic. Their median modeled/observed
SPRUCE porewater-CH4 ratios were `0.00083`, `0.0920`, `0.0920`, `0.1038`, and
`0.1269`, respectively. Thus the phase conversion raises the inferred
concentrations modestly but does not remove the discrepancy: even the nearly
sealed member remains about eightfold low after one year. All fluxes and stocks
reproduced the original experiment to numerical precision, and direct versus
partition-reconstructed porewater diagnostics agreed within
`5.8e-6 mol m-3`.

Because the corrected-threshold member gained `2.13 g C m-2` of methane during
its first year, it is not near a repeating seasonal storage state. The site
forcing spans 2015--2023, so continue it in complete nine-year forcing cycles
rather than choosing an arbitrary endpoint or comparing meteorologically
different adjacent years. Treat it as stabilized when two successive cycles
have end-of-cycle methane stocks within `1%` (and `0.05 g C m-2`), cycle-mean
annual production and surface fluxes within `2%`, and corresponding seasonal
porewater profiles within `5%` at the observed depths and dates. Only a
stabilized profile should be used to tune production or oxidation. The
immediate finding remains that the current high-flux/low-inventory result is
dominated by the too-low ebullition threshold, not oxidation kinetics or a
failure of inventory accounting.

An 18-year continuation from the corrected-threshold member's year-2 restart
subsequently supplied two complete forcing cycles (simulation years 2--10 and
11--19). The first and second cycle ended with `0.07716` and
`0.07575 g C m-2` of CH4. Their absolute difference (`0.00140 g C m-2`) passed
the absolute stock criterion, but the `1.82%` relative difference narrowly
failed the relative criterion. Mean production was already forcing-repeatable
(`5.962` versus `5.986 g C m-2 yr-1`, a `0.39%` difference), whereas mean
surface flux fell from `2.199` to `0.1288 g C m-2 yr-1` and the cycle-mean
porewater profile changed by `89.8%` in vector norm. Median modeled/observed
porewater CH4 declined from `0.0670` to `0.00659`.

This is a microbial-guild establishment transient, not merely slow filling of
the CH4 reservoir. From restart years 2 to 11, aerobic methanotroph carbon rose
from `1.69e-11` to `0.0339 g C m-2`, and anaerobic methanotroph carbon rose from
`4.36e-5` to `0.170 g C m-2`. At year 20 they were `0.0254` and
`0.173 g C m-2`, respectively. In the second cycle, mean production was
`5.986`, mean oxidation was `5.857`, mean surface flux was `0.1288`, and the
total nine-year storage change was only `-0.00140 g C m-2`. The second cycle
looks close to a repeating state internally, but the first-versus-second-cycle
test formally fails because the first cycle contains guild growth. A third
nine-year cycle is required to compare two post-establishment cycles before
declaring convergence or calibrating the low porewater concentrations.

### Mature methane-pathway audit and biological sensitivity

A pathway-resolved year-20 continuation established that methane production
was entirely acetoclastic (`5.764 g C m-2 yr-1`; hydrogenotrophic production
was numerically zero). Total oxidation was `5.624`, partitioned into only
`0.359` aerobic oxidation and `5.266 g C m-2 yr-1` anaerobic oxidation (AOM).
The realized/pre-O2 aerobic-oxidation ratio was `0.9993`, ruling out the shared
O2 limiter as the cause of low aerobic oxidation. Net surface flux was
`0.140 g C m-2 yr-1`. Both the production-pathway and oxidation-pathway sums
closed, and the full methane budget residual was `0.0002 g C m-2 yr-1` at
history-output precision.

An eight-member, nine-year sensitivity ensemble then started every case from
the same mature year-20 restart and used the corrected `50 mmol m-3` transport
threshold. The enclosure mean again excludes fen.

| Change | Production | Aerobic oxidation | AOM | Net flux | Final CH4 stock | Median model/observed porewater CH4 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | 6.031 | 0.413 | 5.488 | 0.130 | 0.074 | 0.0071 |
| acetate production maximum x4 | 21.775 | 1.329 | 20.160 | 0.285 | 0.079 | 0.0075 |
| `K_acetoclastic=5 mmol C m-3` | 0.257 | 0.020 | 0.201 | 0.033 | 0.105 | 0.0038 |
| AOM growth rate x1/2 | 6.031 | 1.244 | 0.523 | 3.950 | 2.903 | 0.0913 |
| AOM growth rate x1/4 | 6.031 | 1.276 | 0.079 | 4.363 | 2.901 | 0.0913 |
| AOM `K_CH4=15 mmol m-3` | 6.032 | 0.882 | 4.911 | 0.133 | 1.021 | 0.0357 |
| AOM O2-inhibition scale `0.46` | 6.031 | 0.454 | 5.449 | 0.129 | 0.076 | 0.0069 |
| AOM off | 6.032 | 1.285 | 0 | 4.434 | 2.902 | 0.0913 |

Flux terms are nine-year annual means in `g C m-2 yr-1`; stock is the final
year-29 restart value in `g C m-2`. Hydrogenotrophic production remained zero
in every member. The baseline third cycle is consistent with the prior mature
cycle. Quadrupling the acetate production ceiling alone mostly grew acetate
methanogens and anaerobic methanotrophs: AOM rose nearly in step with
production, so neither storage nor porewater concentration improved
materially. Reducing the AOM growth parameter below its persistence threshold
instead collapsed the guild, raised final-year surface flux to about
`4.65 g C m-2 yr-1`, and increased methane storage, but the porewater median
still reached only about 9% of observations. Increasing the AOM methane
half-saturation produced an intermediate storage response without increasing
the cycle-mean flux.

These tests reject an aerobic-O2-limitation explanation and show that neither
an unchanged eight-parameter production screen nor single-parameter AOM tuning
is sufficient. The next sensitivity must cross an acetoclastic-production
control with a graded AOM control and evaluate porewater profiles, seasonal
flux, storage drift, and both guild stocks together. AOM-off is a diagnostic
bound, not a scientifically acceptable calibration. The strong guild
extinction threshold also requires timestep and longer-cycle checks before any
AOM-growth value is selected.

The recommended crossed sensitivity was then run for a complete nine-year
forcing cycle. It combined acetate-production maxima of `13.44`, `26.88`, and
`53.76 mmol C m-3 d-1` with AOM CH4 half-saturation values of `1.5`, `5`, and
`15 mmol CH4 m-3`.

| Acetate-production maximum | AOM K_CH4 | Production | AOM | Net flux | Final CH4 stock | Median model/observed porewater CH4 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 13.44 | 1.5 | 6.031 | 5.488 | 0.130 | 0.074 | 0.0071 |
| 13.44 | 5 | 6.031 | 5.285 | 0.128 | 0.277 | 0.0163 |
| 13.44 | 15 | 6.032 | 4.911 | 0.133 | 1.021 | 0.0357 |
| 26.88 | 1.5 | 11.726 | 10.801 | 0.187 | 0.077 | 0.0067 |
| 26.88 | 5 | 11.726 | 10.566 | 0.184 | 0.250 | 0.0198 |
| 26.88 | 15 | 11.727 | 10.030 | 0.215 | 0.940 | 0.0426 |
| 53.76 | 1.5 | 21.775 | 20.160 | 0.285 | 0.079 | 0.0075 |
| 53.76 | 5 | 21.778 | 19.894 | 0.286 | 0.242 | 0.0217 |
| 53.76 | 15 | 21.780 | 18.988 | 0.620 | 0.738 | 0.0517 |

Flux terms are nine-year annual means in `g C m-2 yr-1`; stocks are year-29
restart values in `g C m-2`. Increasing production alone again increased AOM
biomass and oxidation almost proportionally. Increasing AOM `K_CH4` provided a
smooth increase in methane storage and porewater concentration, unlike the
growth-rate extinction experiment, but the strongest combination still
reached only `5.2%` of the observed porewater median. It also retained
`18.99 g C m-2 yr-1` of AOM, so increased production cannot overcome the
current unconstrained AOM response.

The generated baseline confirms AOM maximum growth and mortality rates of
`0.096` and `0.048 d-1`. The earlier half-growth member therefore set maximum
growth exactly equal to mortality before CH4, temperature, pH, and O2
limitations; extinction was a mathematical consequence, not a useful smooth
calibration response. More importantly, the current AOM equation has no
represented electron acceptor: it is limited by CH4, temperature, pH, and O2
inhibition but not sulfate, nitrate, or ferric iron supply. Further tuning of
production or `K_CH4` should pause until the intended AOM rate units are
reconfirmed and a scientifically defensible electron-acceptor constraint,
empirical AOM capacity, or documented AOM-off mode is selected.

### Narrow AOM-growth sweep and literature provenance

A four-member sweep tested whether growth rates moderately above the active
`0.048 d-1` mortality rate could reduce AOM without eliminating its biomass.
All cases retained the baseline acetate-production maximum and AOM
`K_CH4=1.5 mmol m-3` and started from the mature year-20 restart.

| AOM maximum growth (d-1) | AOM | Aerobic oxidation | Net flux | Final CH4 stock | Final AOM biomass | Median model/observed porewater CH4 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.060 | 4.980 | 0.811 | 0.132 | 1.050 | 0.150 | 0.0321 |
| 0.072 | 5.293 | 0.583 | 0.129 | 0.309 | 0.182 | 0.0132 |
| 0.084 | 5.419 | 0.480 | 0.130 | 0.105 | 0.184 | 0.0100 |
| 0.096 | 5.488 | 0.413 | 0.130 | 0.074 | 0.181 | 0.0071 |

Flux terms are first-cycle annual means in `g C m-2 yr-1`; stocks are year-29
restart values in `g C m-2`. The response is smooth and none of these rates
caused guild extinction. Because the `0.060` and `0.072 d-1` members gained
`0.975` and `0.233 g C m-2` CH4 in the first cycle, both were continued through
a repeated forcing cycle. Their second-cycle storage changes fell to `0.068`
and `0.015 g C m-2`, while final AOM biomass was `0.143` and
`0.182 g C m-2`. The corresponding second-cycle AOM rates were `5.149` and
`5.422`, net fluxes were `0.143` and `0.134 g C m-2 yr-1`, and porewater ratios
were `0.0312` and `0.0137`. Thus lowering growth moderately above mortality
does create persistent intermediate concentrations, but it does not remove
the profile discrepancy or substantially increase emissions.

This behavior follows from the biomass equation. At a persistent state,
environmentally limited specific growth balances mortality. Lowering maximum
growth therefore raises the CH4 concentration needed for that balance and can
reduce the supporting biomass, while total AOM continues to adjust toward the
available methane supply. A value extremely close to the persistence boundary
could create a sharp concentration response, but would be vulnerable to
seasonal forcing, timestep dependence, and extinction.

The parameter provenance is weaker than the phrase "adopted from Xu et al."
can imply. Xu et al. (2015) Table 1 reports generic/aerobic methanotroph growth,
death, and yield of `0.15 d-1`, `0.005 d-1`, and `0.40`; it does not give a
separate AOM growth/death pair. It also concludes that AOM was probably a minor
incubation contribution because no iron or other electron acceptor was
observed. Ricciuto et al. (2021) Table 2 gives generic methanotroph growth and
death ranges of `0.0064--0.0096` and `0.0016--0.0024`, with optima of `0.008`
and `0.002`, and gives only an AOM-yield range (`0.12--0.18`, optimum `0.15`).
It does not report separate optimized AOM growth, death, `K_CH4`, or electron-
acceptor parameters, and AOM controls were not among the five parameters found
important to surface flux. The later CLM-SPRUCE AOM pair and its ELM unit
conversion must therefore remain explicitly labeled as weakly constrained.

References: [Xu et al. (2015)](https://doi.org/10.1002/2015JG002935) and
[Ricciuto et al. (2021)](https://doi.org/10.1029/2019JG005468).

### Literal AOM-rate test and legacy effective limiters

The later CLM-SPRUCE parameter file gives AOM growth and death as `0.004` and
`0.002`, and the source declarations label both quantities `d-1`. Git history
shows that this pair was introduced together in the January 2018 parameter
change and retained thereafter. The earlier `x24` intended-unit hypothesis is
therefore not supported for this particular pair: `0.004/0.002 d-1` is the
better literal scientific interpretation, although it cannot reproduce the
legacy executable because that executable multiplied the values by a timestep
in seconds.

A four-cycle continuation from the same productive year-20 restart tested the
literal pair without changing any other parameter. It initially increased
methane storage and emission because AOM biomass established slowly, but it did
not solve the equilibrium profile discrepancy:

| Nine-year window | Production | Aerobic oxidation | AOM | Net flux | End CH4 stock | End AOM biomass | Median model/observed porewater CH4 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| years 20--28 | 6.031 | 1.223 | 0.968 | 3.557 | 2.626 | 0.716 | 0.0852 |
| years 29--37 | 6.184 | 1.022 | 3.870 | 1.534 | 0.441 | 2.783 | 0.0135 |
| years 38--46 | 6.258 | 0.542 | 5.614 | 0.138 | 0.122 | 3.606 | 0.00740 |
| years 47--55 | 6.307 | 0.488 | 5.681 | 0.139 | 0.110 | 3.764 | 0.00685 |

Flux terms are annual means in `g C m-2 yr-1`; stocks are `g C m-2`. The final
window changed total CH4 storage by only `-0.012 g C m-2`, so the delayed AOM
response had effectively caught up. Lowering both rates by 24 therefore
changes the establishment timescale and the biomass required to support the
sink, but not the eventual tendency for AOM to consume nearly all production.

The source audit found no omitted electron-acceptor state or explicit AOM
capacity limit. It did find three legacy effects that limited or deprioritized
AOM and are not reproduced by the conservative ELM kernel:

1. The CLM-SPRUCE implementation multiplied nominal per-day growth and death
   constants directly by the timestep in seconds. With a 1,800-second step,
   the `0.002` mortality term can remove more than the current biomass in one
   call; the subsequent `MFGbiomin` clamp repeatedly resets the guild to its
   floor. This is an accidental timestep-dependent limiter, not a defensible
   rate conversion.
2. The legacy substrate-consumption rate already includes `pHeffect`, and the
   biomass-growth calculation multiplies the resulting uptake by
   `pHeffect` again. Mortality receives it once. At the SPRUCE fallback pH of
   4.5, the response is approximately `0.306`, so this extra factor strongly
   suppresses net guild establishment. ELM applies the environmental response
   once to uptake and derives growth from yield, avoiding that duplicate.
3. In the legacy saturated path, aerobic oxidation and plant/ebullition losses
   are applied before AOM. In the unsaturated path, AOM potential is evaluated
   before current-step methane production, while transport losses are still
   removed before the AOM tendency. ELM instead lets aerobic oxidation and AOM
   compete for initial plus same-step production before surface transport.
   This gives AOM earlier access to methane and can matter greatly when the
   stored concentration is small.

The legacy AOM Q10 expression also subtracts a Celsius value (`13.5`) directly
from Kelvin soil temperature, which accelerates rather than limits AOM and is
already corrected in ELM. Because the old result reflects interacting unit,
pH, temperature, clipping, and operator-order defects, these behaviors should
not be restored as a bundled “parity” switch. The next clean attribution test
should preserve `0.004/0.002 d-1` and separately test (a) the duplicate legacy
pH factor on biomass growth and (b) reaction-versus-transport ordering. A
scientific production configuration still requires an electron-acceptor or
empirical AOM-capacity formulation.

### Interim AOM-off porewater calibration

Pending an electron-acceptor-limited AOM formulation, the SPRUCE profile
experiments set
`microbe_methane_anaerobic_methanotroph_growth_rate=0`. This is an explicit
interim experiment setting, not a new global parameter default. All cases
started from the same productive year-20 restart, ran for two complete
nine-year forcing cycles, and were evaluated over years 29--37. Aqueous
transport was off; CLM-style DOM and acetate relaxation and the observed-DOM
calibration were on. Porewater comparisons use the hummock--hollow enclosure
mean and exclude the fen.

The first crossed ensemble varied acetoclastic production and the shared CH4
transport threshold. Raising the threshold increased storage much more than it
changed production. At the original production maximum of
`13.44 mmol C m-3 d-1`, thresholds of `50`, `250`, and `1000 mmol m-3`
produced median modeled/observed concentration ratios of `0.091`, `0.325`,
and `0.733`. At `1000 mmol m-3`, the annual surface flux was only
`0.356 g C m-2 yr-1`; diffusion supplied `0.356`, plant transport was zero,
and ebullition supplied only `0.0003`. Doubling production increased the
ratio to `1.017` and flux to `3.949`, but also made ebullition
`3.462 g C m-2 yr-1`. The threshold therefore cannot be treated as an
independent concentration-fitting parameter: it controls both ebullition and
the minimum aerenchyma release concentration and changes pathway partitioning.

An aerobic-kinetics ensemble then compared `K_CH4/K_O2` pairs of `1/4`,
`2.5/500`, and `1000/4000 mmol m-3`. The execution-parity pair effectively
eliminated aerobic oxidation and was rejected even where it reduced profile
RMSE. The `2.5/500` pair retained aerobic oxidation and was used for a refined
production-by-threshold ensemble. Its best member used an acetate-production
maximum of `20.16 mmol C m-3 d-1` and a transport threshold of
`750 mmol m-3`:

| Diagnostic | Years 29--37 result |
| --- | ---: |
| Gross CH4 production | 9.170 g C m-2 yr-1 |
| Aerobic oxidation | 5.328 g C m-2 yr-1 |
| Surface flux | 3.831 g C m-2 yr-1 |
| Diffusive surface flux | 1.093 g C m-2 yr-1 |
| Ebullition | 2.738 g C m-2 yr-1 |
| Date-matched median model/observed porewater CH4 | 0.879 |
| Date-matched RMSE | 0.219 mmol L-1 |

Layer diagnostics showed why the residual error was concentrated near the
surface. The area-weighted O2 availability scalar was about `0.5--1` in the
upper `0.25 m` and nearly zero below about `0.4 m`. Deep concentrations were
therefore already close to observations, while upper-profile concentrations
were depleted by aerobic oxidation.

A final narrow sweep increased only aerobic `K_CH4`, holding production,
threshold, `K_O2=500`, and AOM-off fixed:

| Aerobic K_CH4 (mmol m-3) | Aerobic oxidation | Surface flux | Final aerobic biomass | Median model/observed | Date-matched RMSE |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2.5 | 5.328 | 3.831 | 0.4303 | 0.879 | 0.2186 |
| 5 | 4.904 | 4.259 | 0.3926 | 0.891 | 0.2166 |
| 10 | 3.947 | 5.259 | 0.3051 | 0.915 | 0.2132 |
| 12.5 | 3.028 | 6.229 | 0.2728 | 0.929 | 0.2116 |
| 15 | 1.514 | 7.733 | 0.2399 | 0.969 | 0.2093 |
| 17.5 | 0.288 | 8.866 | 0.0442 | 1.063 | 0.2016 |
| 20 | 0.0497 | 9.092 | 0.0060 | 1.070 | 0.2009 |
| 25 | 0.0017 | 9.138 | 0.000085 | 1.070 | 0.2008 |

Fluxes and oxidation are in `g C m-2 yr-1`; biomass is `g C m-2`. The
apparent optimum above `K_CH4=15` is caused by a sharp aerobic-methanotroph
population collapse, not a credible smooth kinetic improvement. The
`K_CH4=15` member is the provisional upper bound that improves the profile
while retaining a substantial methanotroph population. It gives date-specific
median model/observed ratios of `1.05`, `0.969`, and `0.747` for DOY 121, 182,
and 244. It still underestimates the observed upper `0.1--0.25 m` late-season
concentrations and its `7.73 g C m-2 yr-1` surface flux must be checked against
flux observations before adoption.

These results do not justify replacing the standard parameter values. The
next scientific checks are a longer continuation of the provisional member,
independent surface-flux validation, and a model formulation that prevents
guild persistence from becoming an abrupt calibration switch. AOM should
remain available through its parameter and should be re-enabled only after a
defensible electron-acceptor or empirical-capacity constraint is implemented.

### Observed acetic-acid profile calibration

The SPRUCE acetic-acid observations in
`CLM_SPRUCE/scripts/UQ/constraints_ch4/CACES.txt` cover depths of 0.15--2.0 m
on three sampling dates in 2013--2014. They range from approximately
`0.00245` to `0.252 mmol acetate L-1`. The model history variables
`MM_ACETATE_C_SAT` and `MM_ACETATE_C_UNSAT` are bulk-soil acetate carbon in
`g C m-3`, not porewater acetate molecules. A valid comparison first divides
the area-weighted bulk inventory by the matching liquid-water fraction derived
from `SOILLIQ`, then divides by `2 * 12.011 g C mol-1` because one acetate
molecule contains two carbon atoms. Directly relabeling the history field as
`mmol L-1` is invalid.

After that correction, the old test settings
`K_acetoclastic=0.05 mmol C m-3` and
`acetate_feedback_half_saturation=0.1 mmol C m-3` produced a median
modeled/observed concentration ratio of only about `0.0014`. Replacing the
forced `0.99` saturated area with ELM's diagnosed saturated fraction changed
that ratio by less than 1%, showing that saturated-area mapping was not the
cause of the acetate deficit.

The source audit identified a coupled concentration-unit error in the run
parameters. CLM-SPRUCE first converts acetate carbon to `mol C m-3` and then
compares it directly with the raw `m_dKCH4ProdAce=0.05`; the strict ELM
execution-parity value is therefore `50 mmol C m-3`. The active acetate
feedback literal has the same concentration mismatch: raw `0.1 mol C m-3`
maps to `100 mmol C m-3`. Correcting only the methanogenesis half-saturation
while retaining feedback at `0.1` starved the acetoclastic guild and eventually
pinned its biomass at the lower bound. Both concentration scales must be
changed together for a meaningful test.

An AOM-off calibration ensemble used ELM's diagnosed saturated area,
CLM-style DOM and acetate relaxation, observed-DOM restoring, an acetate
production maximum of `20.16 mmol C m-3 d-1`, and aerobic
`K_CH4/K_O2=15/500 mmol m-3`. Aqueous transport was off. The observations were
matched by day of year to a repeated 2015--2023 forcing cycle; this evaluates
seasonal profile shape but is not an exact meteorological match to the
2013--2014 samples. The following values are means over the final complete
nine-year forcing cycle:

| Acetoclastic K (mmol C m-3) | Median modeled/observed acetate | Log10 profile RMSE | CH4 production | CH4 oxidation | Surface CH4 flux | Acetoclastic biomass |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 50 | 0.828 | 0.419 | 8.641 | 2.630 | 6.026 | 0.1385 |
| 75 | 1.113 | 0.385 | 7.444 | 1.862 | 5.589 | 0.1087 |
| 100 | 1.305 | 0.398 | 6.515 | 1.195 | 5.324 | 0.0894 |

Fluxes are in `g C m-2 yr-1`; biomass is the enclosure- and profile-mean
concentration in `g C m-3`. All three members used
`acetate_feedback_half_saturation=100 mmol C m-3`. The `K=75` member had the
lowest cycle-mean profile error and provides a provisional SPRUCE calibration.
At its final year the median modeled/observed ratio was `0.803`, production was
`8.269`, oxidation was `2.496`, and surface flux was
`5.329 g C m-2 yr-1`, demonstrating material interannual variability even
after the long continuation.

These results separate two recommendations. `K=50` and feedback `=100` are the
strict later-code execution-unit mapping. `K=75` and feedback `=100` are a
site-level observational calibration under the experimental AOM-off,
observed-DOM-relaxation configuration; they are not yet justified as global
defaults. Before adoption, repeat the comparison with independent surface CH4
flux data, verify the acetate-carbon and acetic-acid measurement semantics with
the data providers, test additional sites, and re-evaluate the pair after a
scientifically constrained AOM formulation and prognostic DOM transport are
selected.

## Validation experiments and diagnostics

### Preferential-flow rainfall sensitivity

`use_microbe_dom_preferential_flow` is a default-off event-transport
hypothesis. When enabled with physical aqueous transport,
`microbe_methane_dom_preferential_flow_fraction` sets the fraction of positive
infiltration used to sample mobile DOM from all layers above the diagnosed
water table. The sampled C, N, and P are removed conservatively and deposited in the
first layer intersecting the diagnosed water table without exchange with
intermediate layers. Fractions 0.1 and 0.3 are uncalibrated sensitivity values.
An earlier diagnostic version imposed a 0.3 m donor-depth cap; that artificial
cap was removed so the donor domain follows the hydrologic unsaturated depth.
Earlier development runs called this option
`use_microbe_dom_complete_bypass`; that name was replaced because only a
specified event fraction follows the fast pathway.

The current implementation represents a simplified fast-domain case. It does not add water to the
hydrologic state and does not resolve fast-domain storage, travel time, matrix
exchange, or changes in water-table position. It must therefore be used only
to ask whether event-scale bypass could plausibly correct the modeled
surface-to-deep DOM gradient. A production implementation would need a
dual-domain water and solute budget and observational constraints on connected
mobile porosity, exchange, and event activation.

A matched five-year continuation from the year-11 cold-start/compaction
restart tested fractions 0, 0.10, and 0.30. The 0.30 case transferred 17.75 g C
m-2 yr-1 in year 5, reduced DOC by 9--15% in the upper 0.06 m, and increased
DOC by 13%, 22%, and 50% at 0.21, 0.37, and 0.62 m, respectively. The absolute
0.62 m concentration remained only 2.80 g C m-3 water and the approximately
1.04 m concentration remained near 0.003 g C m-3 water. CH4 production,
oxidation, and surface flux were 3.86, 3.96, and 1.38 g C m-2 yr-1 versus
3.82, 3.68, and 1.54 in the control; NPP changed from 297.3 to 295.6 g C m-2
yr-1. Preferential-flow fraction alone therefore improved the profile in the
expected direction but did not repair the deep-DOC deficit. During year 5 the
enclosure-mean water table remained shallower than 0.227 m, so removal of the
former 0.30 m donor cap did not affect this particular comparison.

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

For concurrent one-rank Docker runs, set
`OMPI_MCA_hwloc_base_binding_policy=none` before applying `taskset` to
`case.submit`. OpenMPI can otherwise rebind every actual `e3sm.exe` rank to
core 0 even when the wrapper processes were assigned distinct cores. Confirm
placement on the live model PIDs with `taskset -pc`; wrapper affinity alone is
not evidence that the simulations are running on separate CPUs.

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
