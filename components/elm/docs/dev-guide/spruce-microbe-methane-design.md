# CLM-SPRUCE microbial decomposition and methane integration design

Status: Phase 0 harness and Phases 1-2 implementation complete; Phase 3 Steps
1-5 are implemented on their feature branch. Step 5 makes the executing
methane backends mutually exclusive and connects revised storage and surface
fluxes to ELM balances, budgets, and atmosphere exchange. A default-off
compatibility option retains the Phase 2 microbial-cascade plus legacy-CH4Mod
mode for controlled attribution. The archived
CLM-SPRUCE scientific reference, exact-restart comparison, and SPRUCE science
validation remain pending. Docker/CIME build, five-day enabled/disabled smoke
tests, and a paired 50-year US-MOz run pass numerically with carbon balance at
roundoff. The long run exposes large DOM and soil-C shifts and unconstrained
anaerobic-methanotroph growth, which remain scientific calibration blockers.
Date: 2026-09-18

## 1. Executive summary

This change should be implemented as a selective port and refactor, not as a Git
merge or rebase. The CLM-SPRUCE source predates the current ELM architecture and
uses global `clmtype` storage, CPP selection, positional text parameters, and
direct state mutation. Current ELM already contains a substantially evolved
`CH4Mod`, modern typed state, restart and history infrastructure, C-N-P
decomposition, topounit hydrology, and atmosphere coupling. Replacing those
facilities with the older CLM-SPRUCE implementation would regress current ELM
capabilities.

The proposed design adds one opt-in science switch and independently selectable
transport suboptions:

```fortran
use_microbe_methane = .false.
use_legacy_ch4_with_microbe = .false.
use_elm_microbe_methane_transport = .false.
use_microbe_nonbog_lateral_gas_transport = .true.
```

`use_lch4` remains the existing umbrella switch for methane. The default-off
`use_legacy_ch4_with_microbe` compatibility option retains the microbial
decomposition cascade while selecting legacy `CH4Mod`; it permits controlled
separation of decomposition and revised-methane effects. When
`use_microbe_methane` is omitted or false, ELM follows its present call graph,
pool topology, parameter reads, restart schema, history registration, and
floating-point operations. When both switches are true, the new switch selects:

1. the CLM-SPRUCE-style microbial CTC cascade, adding DOM, bacteria, and fungi
   to the standard decomposition pools; and
2. a revised methane backend with acetate, explicit methanogenic and
   methanotrophic guilds, CH4/O2/CO2/H2 inventories, and gas transport.

Within revised methane, the default false transport suboption uses the
CLM-Microbe aqueous Fickian mapping. Setting it true selects the alternative ELM
multiphase mapping; it has no valid standalone meaning when revised methane is
disabled.

The existing `CH4Mod` remains intact and remains the default methane backend.
The two methane backends must never run on the same soil column and timestep.

The primary production target is the E3SM-Peatlands C-N-P, relative-demand,
vertical CTC configuration. C:P behavior for the three new decomposition pools
is therefore part of the initial integration, not a deferred enhancement.

### 1.1 Science description

The scientific purpose of this option is to connect decomposition, dissolved
organic matter, microbial biomass, redox state, and methane cycling in one
conservative soil-carbon system. Standard ELM represents litter and soil
organic matter as a cascade of pools whose turnover is controlled by
temperature, moisture, depth, and nutrient availability. Carbon leaving a
donor pool is either respired or transferred to another litter or SOM pool.
That structure is efficient for long-term carbon storage, but it does not
explicitly represent the soluble organic substrates that feed anaerobic
metabolism, nor does it distinguish living bacterial and fungal biomass from
the SOM pools they help create. The microbial option adds those missing links
while retaining the existing ELM litter and SOM hierarchy.

DOM, bacteria, and fungi are added as ordinary vertically resolved C-N-P
decomposition pools. Litter and SOM can solubilize into DOM, and they can be
assimilated by bacteria or fungi. The allocation between bacterial and fungal
uptake depends on substrate stoichiometry and the parameterized microbial C:N
traits. Carbon-use efficiency determines how much assimilated substrate
becomes microbial biomass and how much is respired. Bacterial and fungal
turnover then returns carbon to DOM, transfers microbial residues into SOM1-4,
or releases CO2. DOM can likewise be taken up by microbes or stabilized into
SOM. The result is a feedback loop in which plant litter and SOM supply
microbial substrates, microbial growth temporarily retains carbon and
nutrients, and microbial death supplies both relatively available DOM and
potentially persistent SOM residues.

These pathways make decomposition stoichiometry more explicit. Each transfer
attempts to move C, N, and P through the same generic cascade. If donor and
receiver C:N or C:P ratios differ, ELM's nutrient competition machinery
mineralizes excess organic nutrients or immobilizes NH4 and solution phosphate.
The DOM phosphorus pool is dissolved *organic* P and is therefore distinct
from `solutionp_vr`, which represents dissolved inorganic phosphate. DOM can
move with the organic pool's vertical-transport rule; inorganic phosphate
continues to use ELM's existing sorption, transport, and leaching processes.
This distinction permits rapid dissolved-organic cycling without conflating it
with plant-available mineral P.

The present cascade is microbial-explicit in pool topology, but it is not yet
an enzyme or population-control model. Litter and SOM turnover rates retain the
standard environmental controls and are not multiplied directly by total
bacterial or fungal biomass. Microbial biomass affects carbon and nutrient
trajectories through uptake allocation, carbon-use efficiency, its own
turnover, residue formation, and DOM recycling. Consequently, indirect
substrate and nutrient feedbacks are possible, but a strong biomass-driven
priming response to fresh carbon is not an explicit mechanism in this first
port. Adding enzyme production, biomass-dependent depolymerization, or a
formal priming formulation would be a separate science change with its own
parameter and validation requirements.

DOM provides the main bridge from the decomposition cascade to revised
methane chemistry. Under sufficiently moist, warm, and reducing conditions,
DOM-C is fermented to acetate, CO2, and H2. Acetate and H2/CO2 then support two
gross methane-production pathways: acetoclastic methanogenesis and
hydrogenotrophic methanogenesis. Their biomasses are represented by separate
functional guilds, distinct from the general bacterial and fungal pools.
Methane is consumed by aerobic methanotrophs where both CH4 and O2 are
available and by an anaerobic methanotroph guild under low-O2 conditions. All
four methane guilds grow and die explicitly; mortality carbon returns to DOM.
The initial parity model treats these guilds as carbon-only, so their growth
does not yet impose an additional N or P demand beyond the DOM transformations
handled by the C-N-P transaction.

Hydrology and oxygen determine which methane pathways can operate and whether
their products reach the atmosphere. Each layer contains saturated and
unsaturated subarea state. Changes in inundated fraction conservatively move
existing acetate, gases, and guild biomass between those partitions. CH4, O2,
CO2, and H2 state is stored as bulk-soil inventory density. Gas-phase diffusion
operates in air-filled pores above the water table, while Henry-law-equilibrated
aqueous diffusion operates in saturated layers. Gases can also exchange through
roots and aerenchyma, and methane above its solubility threshold can escape by
ebullition. Snow, ponded water, the top soil layer, and the atmospheric boundary
contribute resistances to surface exchange. Thus net CH4 emission is not equated
with production: it is the remainder after oxidation, temporary gas storage,
downward or upward redistribution, and the three transport pathways.

Oxygen is shared with the rest of ELM rather than treated as an isolated
methane control. In each saturated and unsaturated layer partition, potential
heterotrophic respiration, root respiration, nitrification, aerobic acetate
oxidation, and aerobic CH4 oxidation form one demand. A common stress scales
the demand when it exceeds the finite dissolved-O2 inventory. Nitrification
enters at two moles of O2 per mole of N. Following legacy `CH4Mod`, ELM carbon
and nitrogen fluxes have already used the preceding timestep's stress, while
the current stress scales microbial aerobic reactions and is published through
the established ELM anoxia interface for the next decomposition step. This
creates the intended lagged feedback: wetness restricts atmospheric O2 supply,
redox limitation alters decomposition and nitrogen transformations, DOM and
fermentation products accumulate, methane production increases, and oxidation
and transport determine the emitted fraction.

Carbon accounting differs deliberately from the legacy `CH4Mod` diagnostic
formulation. In the revised model, methanogenesis and methane oxidation are
internal conversions among represented carbon pools; neither is itself an
external source or sink. The carbon budget includes DOM through the standard
decomposition pools and adds acetate, guild biomass, bulk-soil CH4-C, and
bulk-soil CO2-C inventories as revised-backend storage. Only net CH4-C and CO2-C crossing
the land-atmosphere boundary are external fluxes. Gross diagnostics retain the
science needed to interpret that net exchange: `MM_CH4_PROD` integrates
acetoclastic plus hydrogenotrophic methanogenesis, while `MM_CH4_OXID`
integrates aerobic plus anaerobic methane consumption. They correspond in
purpose, though not in detailed formulation, to legacy `CH4PROD` and
`FCH4TOCO2`. The acetoclastic contribution to `MM_CH4_PROD` is the reaction
extent multiplied by the non-biomass fraction and the CH4 yield; counting the
whole process extent would incorrectly include biomass and CO2 carbon as CH4.
`MM_CH4_PROD_UNSAT`, `MM_CH4_PROD_SAT`, `MM_CH4_OXID_UNSAT`, and
`MM_CH4_OXID_SAT` are area-weighted contributions that sum to the bulk
diagnostics. `MM_CH4_OXID_AER` and `MM_CH4_OXID_AOM` separate aerobic
methanotrophy from anaerobic methane oxidation; their `_UNSAT` and `_SAT`
fields expose both area-weighted subarea contributions. The aerobic and AOM
bulk fields sum exactly to `MM_CH4_OXID`, and `MM_SAT_FRACTION` records the
area weighting used. `MM_CH4_SURF_DIFF`, `MM_CH4_SURF_AERE`, and
`MM_CH4_SURF_EBUL` expose the three surface pathways. Each has `_UNSAT` and
`_SAT` area-weighted contributions that sum to its bulk field; all use
`g C m-2 s-1` and are positive toward the atmosphere.

The expected land-carbon response is therefore not a uniform increase or
decrease in SOM. Carbon is diverted into rapidly cycling DOM and living
microbial pools, while microbial residues and DOM stabilization can feed slow
SOM. Early spinup may show lower litter or selected SOM stocks because of new
uptake and solubilization paths, together with a sizeable DOM pool if its
production exceeds microbial consumption and vertical export. Long-term SOM3
and SOM4 responses depend on residue routing, nutrient limitation, and the
unaccelerated turnover of the new pools. Accelerated decomposition spinup
remains useful for the inherited slow SOM pools, but bacteria, fungi, DOM, and
the revised methane state retain a spinup factor of one. Their values during an
AD run should therefore be interpreted as coupled fast-state adjustment around
accelerated SOM, not as an independently accelerated equilibrium. Final
spinup and transient tests are required before evaluating stocks, NEE, or CH4
emissions scientifically.

## 2. Source baselines inspected

| Tree | Revision inspected | Relevant source |
| --- | --- | --- |
| E3SM-Peatlands | branch `simplify/bog-zwt-no-perched`, commit `7b0703e07c7c933d53814f0c31856d3aed480d7b` | `components/elm` |
| CLM-SPRUCE | branch `master`, commit `a85800bad2c2a57abff77af36ccec319b202166a` | `models/lnd/clm/src/clm4_5` |
| CLM-Microbe | branch `master`, commit `9c2e0a048bb3799669d32b91e6cc76efb36d4b75` | `models/lnd/clm/src/clm4_5` and `inputdata/lnd/clm2/paramdata/microbepar_in` |

The important CLM-SPRUCE/CLM-Microbe source files are `microbeMod.F90`,
`initmicrobeMod.F90`, `microbeRestMod.F90`, `microbevarcon.F90`, and the
microbial changes in `CNDecompCascadeMod_BGC.F90`, `clm_varpar.F90`,
`pftvarcon.F90`, `CNEcosystemDynMod.F90`, and `clm_driver.F90`.

The important current ELM facilities are `CH4Mod.F90`, `CH4varcon.F90`, the
column C/N/P state and flux types, `CNDecompCascadeConType.F90`,
`DecompCascadeCNMod.F90`, `SoilLittDecompMod.F90`, `EcosystemDynMod.F90`,
`elm_driver.F90`, `elm_instMod.F90`, `readParamsMod.F90`, `restFileMod.F90`,
and `lnd2atmType.F90` / `lnd2atmMod.F90`.

## 3. Goals and non-goals

### Goals

- Add DOM, bacterial biomass, and fungal biomass as vertically resolved ELM
  decomposition pools when the option is enabled.
- Port the CLM-SPRUCE 11-pool, 47-transition microbial CTC cascade while
  preserving current ELM C, N, and P accounting.
- Port the revised methane reaction network and its acetate and functional-guild
  state.
- Connect the revised methane backend to current ELM hydrology, vegetation,
  history, restart, mass-balance, and land-atmosphere exchange interfaces.
- Generalize the source's two-column lateral gas diffusion to ELM's arbitrary
  multi-topounit connection graph without hard-coded hummock/hollow indices.
- Preserve exact current behavior when the new option is disabled.
- Make all scientific parameters named, unit-documented, validated NetCDF
  inputs.
- Establish process-level CLM-SPRUCE parity before calibrating the port or
  changing its equations.

### Non-goals for the initial integration

- Replaying or merging the old CLM-SPRUCE commit history.
- Replacing or cleaning up the current legacy `CH4Mod` path.
- Porting the empty `microbeCN` or `microben2o` stubs.
- Porting `HUM_HOL` hydrology or directly copying its ad hoc lateral BGC code.
  The new backend uses current ELM topounit topology for Fickian gas exchange;
  water-flux-driven DOM, nutrient, and acetate advection remains deferred.
- Supporting FATES, BeTR/sBeTR, PFLOTRAN, Alquimia, or EMI BGC in the first
  release.
- Enabling lake methane production in the revised backend initially.
- Recalibrating parameters during the mechanical/scientific parity port.
- Adding N or P stoichiometry to the four methane functional guilds. Those
  guilds are carbon-only state in CLM-SPRUCE and remain so for initial parity.

## 4. Findings that drive the design

### 4.1 Current ELM methane is an active, modern implementation

Current `CH4Mod` is not a placeholder. It owns typed column/gridcell state,
parameter reads, initialization, history, restart, production, oxidation,
aerenchyma, diffusion, ebullition, and atmosphere-facing outputs. It is called
from `elm_driver.F90` after `AnnualUpdate`, and current coupling exports net CH4
and applies the net-methane correction to CO2 when online coupling is active.

The revised backend must therefore plug into the current lifecycle and output
contract; it must not overwrite `CH4Mod.F90` with the older file.

### 4.2 The CLM-SPRUCE feature is two coupled models

The `MICROBE` CPP path contains both:

- a microbial decomposition cascade with DOM, bacteria, and fungi; and
- a methane reaction/transport model with acetate and four explicit functional
  guilds.

The bacteria and fungi in the decomposition cascade are not aliases for the
acetoclastic methanogens, hydrogenotrophic methanogens, aerobic methanotrophs,
or anaerobic methanotrophs. These are separate prognostic states and must remain
separate in the port.

### 4.3 Positional parameter input requires an explicit baseline

`microbevarcon.F90` declares `nummicrobepar = 88`, reads exactly 88 values from
`./microbepar_in`, ignores the names, and assigns by position. The original
CLM-SPRUCE snapshot inspected for the port has a 111-record file whose records
84-86 are `k_dom`, `k_bacteria`, and `k_fungi`; the reader instead assigns
those positions to `dom_diffus`, `m_Fick_ad`, and `m_dPlantTrans`. It therefore
cannot define the revised-methane runtime values safely. The local CLM-SPRUCE
test branch repairs the parameter file by moving `dom_diffus`, `m_Fick_ad`, and
`m_dPlantTrans` to records 84-86. The decomposition records remain later in the
file for provenance; their active values are read from the PFT parameter
NetCDF, not assigned by this methane reader.

The later CLM-Microbe snapshot has an 89-record file whose first 86 positions
align with the reader through `m_dPlantTrans`. Records 87 (`AOM`) and 88
(`H2maxCH4`) are read into the temporary array but never assigned, and record
89 (`AcemaxCH4`) is not read. For this integration, that file at commit
`9c2e0a048bb3799669d32b91e6cc76efb36d4b75` is the default upstream baseline
for the 38 active Phase 3 parameters it actually supplies. Phase 3 maps those
values by legacy identifier into named standard-ELM variables; it does not
reproduce positional reads at runtime. Before writing the ELM parameter file,
the converter translates the raw numbers from the units used by the active
legacy equations into ELM's declared units. Parameters absent from the file use
an audited declaration or active source literal with distinct provenance.

This resolves which checked-in values the integration tests use, but it does
not make them a production calibration. Archived-run provenance and site
suitability remain scientific-validation tasks. Many
PFT-level microbial-cascade parameters are separately declared in
`pftvarcon.F90` and are not supplied by this text file.

### 4.4 The old module directly mutates and duplicates state

The old methane routine copies DOM/DON out of decomposition arrays into
`cdocs`/`cdons`, advances separate saturated and unsaturated copies, then writes
DOM back to the decomposition arrays and resets DON from a fixed C:N ratio.
Functional biomass and gas arrays are updated in place. This obscures the mass
budget and is especially unsafe in current ELM's C-N-P model.

The port instead uses one authoritative state for every pool and an explicit
tendency/commit boundary.

### 4.5 The old source has ambiguous duplicate invocation sites

CLM-SPRUCE contains `microbech4` calls in both `CNEcosystemDynMod.F90` and
`clm_driver.F90`, and invokes the whole routine separately for bulk C, C13, and
C14 in one path. The port will have one timestep entry point and one bulk
reaction calculation. Isotope transfers will derive from that calculation.

### 4.6 CLM-Microbe structural delta audit

The later CLM-Microbe tree is an incremental evolution of the same model, not
a replacement formulation. Its core `microbeMod.F90` reaction equations,
`microbevarcon.F90` declarations, ecosystem driver, and Century cascade remain
substantially the same as CLM-SPRUCE. The following deltas affect the porting
decision:

- The seven litter/SOM-to-DOM transition respiration fractions changed from
  one to zero. This confirms that they are intended as conservative DOM
  transfers and matches the Phase 2 implementation.
- Dominant-PFT selection now explicitly handles bare-ground or nearly
  unvegetated columns. The ELM adapter independently implements the same safety
  requirement using current patch weights and PFT indices.
- The new fungi-allocation expressions add `(1 / 30)**0.6`. Because `1 / 30`
  is integer division in the source, the term evaluates to zero and does not
  change the old allocation. The port retains the effective equation; any
  intended nonzero fungi preference is a future science change.
- Inundation changed from a layer/water-head split toward an empirical
  `max(0.001, 1 - 0.04*zwt**2)` column fraction. Unsaturated surface exchange
  was also extended to every layer with a `1/j**2` weighting, frozen exchange
  was attenuated to one percent rather than disabled, and a top-five-layer
  atmospheric-CH4 oxidation diagnostic was subtracted from net CH4 flux. These
  are not copied literally: the port uses current ELM water state, conservative
  saturated/unsaturated repartition, and finite-volume transport. The added
  atmospheric-oxidation behavior remains a validation comparison target.
- Both inspected CLM-Microbe snapshots contain explicit vertical Fickian
  pathways: DOC and acetate use `dom_diffus`, while CH4/O2/CO2/H2 use the water
  diffusion coefficient (`Fick_D_w*m_Fick_ad`) throughout the soil profile,
  including unsaturated layers. `HUM_HOL` adds lateral concentration-gradient
  exchange. The later source partly masks weak atmospheric gas coupling by
  resetting a near-surface layer to atmospheric Henry equilibrium each
  timestep. The ELM port keeps DOM C/N/P on the standard SOM transport operator
  and gives acetate its own conservative `dom_diffus` solve. Revised gas
  transport now defaults to the source aqueous Fickian coefficients; the
  air-filled-pore/aqueous mapping and finite snow/pond/topsoil/boundary
  resistance remain available through a namelist option. Section 9.3.2 records
  both mappings and the diagnostics required to validate them.
- The source's lateral gas diffusion assumes columns 1 and 2, fixed 0.75/0.25
  areas, matching layer numbers, and a one-meter exchange distance. The port
  instead reads `regional_target_ti`, `lateral_dist`, topounit weights, and
  elevations from ELM. It exchanges all overlapping absolute-elevation layer
  segments simultaneously, applies one donor limiter across every outgoing
  edge, and conserves gridcell-area inventory. This supports the current
  three-topounit SPRUCE surface without a new neighbor data model.
- The later tree still owns state through global `clmtype` arrays and mutates
  duplicated saturated/unsaturated state in place. Its active calls remain in
  `CNEcosystemDynMod.F90` for bulk C, C13, and C14, while the driver call is
  commented. The ELM port therefore retains its single bulk call, typed state,
  transaction boundary, and derived future isotope transfers.
- The tree includes experiment-specific changes outside the microbial model,
  including fivefold N deposition, forced constant N fixation, and altered
  respiration summaries. Those changes are not module dependencies and are
  deliberately excluded.

No CLM-Microbe source examined here adds C:P parameters for DOM, bacteria, or
fungi, nor does its text file supply the Phase 2 microbial turnover and routing
inventory. Those ELM-CNP inputs remain explicit new or unresolved parameters.

## 5. Configuration contract

### 5.1 Namelist variables

Add to `elm_inparm`:

```fortran
logical :: use_microbe_methane = .false.
logical :: use_legacy_ch4_with_microbe = .false.
logical :: use_elm_microbe_methane_transport = .false.
logical :: use_microbe_nonbog_lateral_gas_transport = .true.
```

`use_legacy_ch4_with_microbe=.true.` selects the retained Phase 2 bridge:
DOM, bacteria, and fungi remain active in the decomposition cascade, while
legacy `CH4Mod` supplies methane, oxygen/anoxia, atmospheric exchange, restart,
and history behavior. Revised methane state and parameters are not allocated or
read in this mode. The option is intended for scientific attribution and
backward-compatible experiments; it defaults false, so the main enabled
configuration continues to use revised methane.

When revised methane is enabled, the transport option defaults false and uses
the CLM-Microbe aqueous Fickian mapping. Setting it true selects the ELM
multiphase air-filled-pore/aqueous mapping. It is invalid to enable the transport
option while `use_microbe_methane` is false.

The default-on `use_microbe_nonbog_lateral_gas_transport` switch preserves the
ported CLM-Microbe behavior on ELM's generalized topounit graph. Setting it
false excludes an edge when either endpoint is a non-bog topounit, while
retaining bog-to-bog gas exchange. The switch affects CH4, O2, CO2, and H2 gas
exchange only; it does not change lateral hydrology, vertical gas transport, or
aqueous DOM C/N/P and acetate transport. The default therefore preserves the
existing enabled-path result.

The module's scientific parameters will be variables in ELM's standard
parameter NetCDF selected by the existing `paramfile` machinery. There is no
second module-specific filename namelist. The added variables are required and
read only when `use_microbe_methane=.true.`; existing ELM parameter files remain
valid when it is false. A build-namelist convenience option such as
`-microbe_methane` may set both `use_lch4=.true.` and
`use_microbe_methane=.true.`, but the namelist remains the authoritative runtime
selection.

Keep `hist_wrtch4diag` for the established backend. Add
`hist_wrtmicrobediag=.false.` only if the full revised-model diagnostic set is
too large for default history. Scientific coefficients do not belong in the
namelist.

### 5.2 Selection truth table

| `use_lch4` | `use_microbe_methane` | `use_legacy_ch4_with_microbe` | Result |
| --- | --- | --- | --- |
| false | false | false | Current no-methane behavior, unchanged |
| true | false | false | Standard decomposition plus current `CH4Mod`, unchanged |
| false | true | either | Configuration error; the microbial cascade requires methane oxygen/anoxia coupling |
| either | false | true | Configuration error; the compatibility option requires the microbial cascade |
| true | true | false | Microbial decomposition plus revised methane backend |
| true | true | true | Microbial decomposition plus retained legacy `CH4Mod` (Phase 2 compatibility bridge) |

A string-valued `methane_model='legacy|microbial'` was considered but rejected
for this port. Adding a default-false boolean makes the backward-compatibility
boundary explicit and does not reinterpret existing user namelists or compsets.

### 5.3 Required configuration when enabled

The build-namelist layer and runtime initialization must both validate:

- `use_cn=.true.`;
- `use_vertsoilc=.true.`;
- `use_century_decomp=.false.` (the CTC cascade);
- `use_elm_microbe_methane_transport=.true.` only when
  `use_microbe_methane=.true.`;
- `use_legacy_ch4_with_microbe=.true.` only when both `use_lch4` and
  `use_microbe_methane` are true;
- the compatibility bridge cannot be combined with revised-methane saturation,
  gas-transport, DOM-relaxation, or physical aqueous-solute options;
- ELM-native BGC is active;
- `use_fates=.false.` and BeTR/sBeTR, PFLOTRAN, Alquimia, and EMI BGC are
  inactive;
- `allowlakeprod=.false.` for the revised backend; and
- the standard ELM parameter file contains the required microbial variables and
  compatible schema metadata.

The first production configuration is C-N-P with `nu_com='RD'`,
`suplnitro='NONE'`, and `suplphos='NONE'`, matching the main
`CNPRDCTCBC` E3SM compsets. C-N and nutrient-supplemented configurations are
useful for CLM-SPRUCE parity tests but are not production-supported until they
receive their own integration tests. `nu_com='ECA'` and `nu_com='MIC'` must fail
fast until explicitly validated.

Soil columns are supported. Lake columns may coexist in a grid, but revised
methane production and state are zero on them. Crop and glacier behavior is not
claimed until parameter coverage and tests exist.

## 6. Bit-for-bit disabled-mode contract

For a fixed machine, compiler, build options, task/thread layout, input data,
and initial condition, omitting the new option and setting it explicitly false
must produce bit-for-bit identical history and restart output to the pre-change
baseline. This applies independently to existing `use_lch4=.false.` and
`use_lch4=.true.` cases.

To make that achievable, disabled mode must satisfy all of the following:

- `ndecomp_pools`, transition counts, indices, allocation bounds, and loop
  ordering retain their current values.
- No microbial parameter variables are required or read; pre-feature standard
  ELM parameter files remain valid.
- No microbial state is allocated, initialized, registered with history, or
  written to restart.
- No old arithmetic expression or reduction order is changed merely to share
  code with the new backend.
- The existing `CH4()` call and its argument order remain the selected call.
- Existing methane parameter defaults and `ch4par_in` semantics remain
  unchanged.
- Branches selecting the new implementation are placed outside existing inner
  numerical loops.
- New pool indices are invalid/sentinel values while disabled and cannot be
  used accidentally.

Compilation of the new modules is acceptable; execution and data-layout changes
on the disabled path are not.

## 7. Proposed architecture

### 7.1 Module boundaries

Add these focused modules rather than one ported monolith:

| Module | Responsibility |
| --- | --- |
| `MicrobeDecompMod` | Define the 11-pool/47-transition microbial CTC topology, dynamic CUE/path fractions, turnover rates, and pool-role metadata |
| `MicrobeMethaneParamsMod` | Conditionally read and validate named variables from the standard ELM parameter file and record their schema/version |
| `MicrobeMethaneType` or `MicrobeMethaneMod` | Own revised methane prognostic state, fluxes, initialization, history, and restart |
| `MicrobeMethaneReactionMod` | Compute DOM-to-acetate/gas reactions, functional-guild growth/death, methanogenesis, and CH4 oxidation tendencies |
| `MicrobeGasTransportMod` | Repartition saturated/unsaturated state and compute diffusion, ebullition, aerenchyma, and surface exchange for CH4/O2/CO2/H2 |
| `MicrobeMethaneStateUpdateMod` | Limit tendencies, commit them once, and expose budget terms to ELM C/N/P balance checks |

The exact file split may be reduced during implementation, but the reaction,
transport, and commit interfaces must remain distinct and testable.

Longer term, common gas transport could be extracted from current `CH4Mod` and
shared. That is not the first implementation step because changing legacy
floating-point code conflicts with the disabled-mode B4B requirement.

### 7.2 Runtime flow

```text
namelist / build-namelist
        |
        +-- option false --> current pool topology and current CH4 call
        |
        `-- option true
              |
              +-- validate CTC + vertical C + native CNP configuration
              +-- initialize 11-pool microbial cascade
              +-- read/validate microbial variables from the standard ELM parameter file
              |
        ELM SoilLittDecomp and C/N/P state updates
              |
        AnnualUpdate (existing ordering constraint)
              |
        revised methane Step, exactly once
              +-- conservative saturated/unsaturated repartition
              +-- compute reaction tendencies
              +-- limit to available substrates/electron acceptors
              +-- transport gases
              +-- commit C/N/P and methane state once
              `-- publish CH4 and CO2-correction budget terms
              |
        existing atmosphere export, history, and balance checks
```

At the existing `elm_driver.F90` methane call site, dispatch outside the kernels:

```fortran
if (use_lch4 .and. .not. is_active_betr_bgc) then
   if (use_microbe_methane) then
      call microbe_methane_vars%Step(...)
   else
      call CH4(...)
   end if
end if
```

`Step` needs the current typed C/N/P state and flux objects in addition to the
forcing, hydrology, temperature, vegetation, and exchange objects already
provided to methane. There must be no second call from `EcosystemDynMod`.

### 7.3 Shared methane output contract

The implemented revised backend publishes distinct column CH4 and CO2 surface
fluxes from its own type. `lnd2atmMod` selects those fields when
`use_microbe_methane=.true.`, aggregates them into the established
`flux_ch4_grc` and `nem_grc` coupler fields, and otherwise reads the unchanged
legacy `ch4_vars%ch4_surf_flux_tot_col` path. Revised internal state therefore
does not masquerade as legacy CH4 state, while coupler field names remain
unchanged.

If this proves too confusing during implementation, introduce a small
backend-neutral `methane_exchange_type`; do not move legacy internal state into
the new type as part of the same change.

The column carbon budget must also receive explicit production, oxidation, and
surface-loss terms. Writing only an atmosphere diagnostic is insufficient for
ELM's carbon balance checks.

## 8. Decomposition pool design

### 8.1 Pool topology

Disabled mode retains current values: 8 pools and 9 transitions for CTC, or 7
pools and 10 transitions for Century. Enabled mode is valid only for CTC and
uses:

| Index | Pool | Authoritative ELM state | Initial CLM-SPRUCE C:N |
| --- | --- | --- | --- |
| 1-8 | LITR1, LITR2, LITR3, CWD, SOM1-SOM4 | Existing decomp arrays | Existing ELM values |
| 9 | bacteria | C/N/P decomp arrays | 5 |
| 10 | fungi | C/N/P decomp arrays | 15 |
| 11 | DOM | C/N/P decomp arrays | 10 |

The historical bacteria and fungi seed stock is `1e-5 g C m-3`; DOM starts at
zero. These become named parameters with validation rather than literals hidden
in topology code. Bacteria, fungi, and DOM also receive explicit C:P parameters
for CNP initialization and transfers. Those C:P values require science signoff
because CLM-SPRUCE did not define them.

`CNDecompCascadeConType` should gain pool-role metadata adequate to distinguish:

- litter;
- coarse woody debris;
- particulate SOM;
- dissolved organic matter; and
- microbial biomass.

At minimum this requires `is_dissolved` and `is_microbial_biomass`; overloading
the existing `is_soil` logical is not sufficient for fire, leaching, vertical
transport, and erosion routing.

### 8.2 Transition topology

The enabled cascade has 47 transitions, grouped as follows:

| Transition group | Count | Flows |
| --- | ---: | --- |
| CWD fragmentation | 2 | CWD to LITR2/LITR3 |
| Substrate uptake | 14 | LITR1-3 and SOM1-4 to bacteria/fungi |
| Direct stabilization | 6 | LITR1 to SOM1, LITR2 to SOM2, LITR3 to SOM3, and SOM1 to SOM2 to SOM3 to SOM4 |
| Microbial residues | 8 | bacteria/fungi to SOM1-4 |
| DOM stabilization | 4 | DOM to SOM1-4 |
| Microbial lysis to DOM | 2 | bacteria/fungi to DOM |
| DOM uptake | 2 | DOM to bacteria/fungi |
| Microbial respiration | 2 | bacteria/fungi to atmosphere |
| Solubilization | 7 | LITR1-3 and SOM1-4 to DOM |

For every donor, path fractions and respired fractions must be validated against
the equations, not merely forced to sum to one. The old cascade uses both path
fractions and pathway respiration/CUE, so a naive sum check can double-count the
respired part. Unit tests will encode the actual donor mass equation.

The original CLM-SPRUCE snapshot assigns a respiration fraction of one to the
seven litter/SOM-to-DOM solubilization paths. Under current ELM cascade
semantics that would respire all of their carbon and deliver none to DOM,
contradicting the named receivers and documented solubilization fractions.
Phase 2 therefore treats these as conservative transfers (`rf=0`). The later
CLM-Microbe snapshot independently changes all seven fractions to zero, so the
implemented behavior now agrees with the linked upstream tree as well as the
intended receiver semantics. It remains an expected difference from golden
vectors made with the older CLM-SPRUCE snapshot.

### 8.3 C-N-P behavior

- DOM, bacteria, and fungi use the standard vertically resolved C/N/P decomp
  arrays as their only bulk state.
- P in the DOM decomp pool is dissolved organic phosphorus (DOP). It is distinct
  from `solutionp_vr`, which is dissolved inorganic soil-solution phosphate;
  neither state aliases or duplicates the other.
- Transfer C, N, and P through standard cascade flux arrays and state-update
  infrastructure wherever possible.
- Cascade mineralization and immobilization provide the explicit exchange
  between DOM-P and solution phosphate required by donor/receiver C:P ratios.
  DOM transport carries its organic P with DOM, while solution phosphate keeps
  ELM's existing mineral-P transport, sorption, and leaching behavior.
- Bacterial/fungal C:N and C:P ratios and DOM C:N and C:P behavior are explicit
  named parameters.
- When methane chemistry consumes DOM-C, associated DOM-N and DOM-P changes are
  emitted as explicit organic-to-mineral fluxes needed to restore the selected
  stoichiometry. They are not assigned by overwriting pool values.
- Acetate and the four methane guilds remain carbon-only in the parity model and
  do not compete for mineral N or P. This limitation must be documented in
  output metadata and revisited as a separate science change.

### 8.4 Interaction with generic ELM processes

Every generic loop over `ndecomp_pools` must be audited. In particular:

- DOM can participate in dissolved vertical transport/leaching; bacteria and
  fungi should not automatically inherit DOM mobility.
- Fire, erosion, and land-use transitions need deliberate pool-role policies.
- Accelerated decomposition spinup must not accelerate living biomass or DOM
  without an explicit rule.
- Dynamic subgrid remapping and restart pool-name mapping must preserve the new
  pools.
- External model interfaces that assume 7 or 8 pools are rejected while this
  option is enabled rather than receiving a silently changed vector.

## 9. Revised methane state and algorithms

### 9.1 Prognostic state owned by the revised backend

For each active soil column and decomposition layer:

- acetate carbon;
- acetoclastic methanogen biomass C;
- hydrogenotrophic methanogen biomass C;
- aerobic methanotroph biomass C;
- anaerobic methanotroph biomass C; and
- bulk-soil CH4, O2, CO2, and H2 inventory densities.

State that genuinely differs between saturated and unsaturated subareas is
stored per partition. Bulk diagnostic fields are derived views, not additional
prognostic copies. DOM/DON/DOP and bacterial/fungal biomass are never duplicated
inside the methane type.

### 9.2 Reaction groups to port

- conversion of plant/SOM-derived DOM to acetate and H2/CO2;
- acetogenesis;
- acetoclastic methanogenesis;
- hydrogenotrophic methanogenesis;
- aerobic CH4 oxidation;
- anaerobic CH4 oxidation;
- growth and mortality of all four functional guilds;
- plant/root O2 demand used by the reaction constraints; and
- pH, temperature, redox, saturation, and substrate limitation functions that
  are active in the reference configuration.

Each process returns named tendencies. Substrate limiting is applied to the set
of competing tendencies before any state is changed. All pools must remain
nonnegative without post-hoc clipping that loses mass.

#### 9.2.1 Step 2 reaction-kernel contract

The Step 2 kernel is a pure, scalar, one-layer calculation. Inputs are immutable
and use `g C m-3 soil` for DOM, acetate, and functional guilds; gas inventory
densities use `mol gas m-3 bulk soil`. Outputs are named per-second tendencies in the corresponding
units plus carbon-basis process rates in `mol C m-3 s-1`. Hydrology supplies
bounded DOM-fermentation and aerobic-acetate-oxidation scalars; defining those
scalars and applying the kernel across saturated/unsaturated state belongs to
Steps 3-4.

For a process extent `R`, the implemented carbon and gas stoichiometry is:

- DOM fermentation: `-1.5 DOM-C + 1 acetate-C + 0.5 CO2 + 1/6 H2`;
- acetogenesis: `-1 CO2-C - 2 H2 + 1 acetate-C`;
- acetoclastic methanogenesis: acetate-C is divided between guild growth and
  the non-biomass remainder, with the named CH4 yield dividing that remainder
  between CH4 and CO2;
- hydrogenotrophic methanogenesis: `4 H2` produces one CH4-C, while an
  additional `y_H2` CO2-C supplies the explicitly represented biomass growth;
- aerobic methane oxidation: consumed CH4-C is divided between methanotroph
  growth and CO2, with the named O2:CH4 ratio applied to the full CH4 uptake;
- anaerobic methane oxidation: consumed CH4-C is divided between anaerobic
  methanotroph growth and CO2; and
- each guild's mortality carbon returns to DOM.

All four guild yields are interpreted consistently as biomass-C per
carbon-substrate process extent. For hydrogenotrophic methanogenesis, the
extent is CH4-C production, so applying `y_H2` makes the named growth-rate
parameter the maximum specific biomass growth rate. The kernel does not retain
CLM-SPRUCE's additional factor of four in hydrogenotroph biomass growth; that
source behavior is recorded as an explicit golden-vector review point.

Water and functional-guild H/O are not state variables, and CLM-SPRUCE does not
identify an electron acceptor for anaerobic methane oxidation. The kernel can
therefore close represented carbon and expose gas stoichiometric demands, but
cannot claim a complete H/O atom budget for biomass synthesis or AOM until
those missing constituents are specified.

The legacy equations contain several unambiguous accounting defects that are
not reproduced: full substrate product plus biomass growth creates carbon,
guild mortality deletes carbon, the above-water-table acetate-oxidation sign
creates O2, and acetogenesis limits on hydrogenotrophic methanogen biomass where
the parameter and documented equation identify CO2. Step 2 instead partitions
substrate carbon with the named yields, returns mortality carbon to DOM, consumes
O2 during aerobic acetate oxidation, and limits acetogenesis with dissolved CO2
rather than methanogen biomass. These changes mean the kernel intentionally
does not reproduce those carbon leaks in raw CLM-SPRUCE vectors.

Potential rates are computed first. A separate pure limiter then stages DOM
conversion; competing H2/CO2 consumers; competing acetate consumers; competing
CH4 consumers; and competing O2 consumers. Each competing group is scaled
proportionally to available substrate over the requested timestep. A final
pure assembly routine converts the limited process rates to tendencies. State
mutation, saturation repartition, transport, and ELM budget commits remain out
of scope for Step 2.

### 9.3 Gas transport

Port diffusion, ebullition, and aerenchyma terms for CH4/O2/CO2/H2 behind a
single transport interface. Reuse current ELM atmospheric forcing,
soil-temperature, porosity, water-filled pore space, root/aerenchyma, and water
table inputs. Do not port old global hydrology pointers or `HUM_HOL` branches.

The saturated area fraction can change between timesteps. Repartition before
reaction and transport using a conservative invariant for each layer and gas:

```text
M_total = V_layer * [(1 - f_sat) * C_unsat + f_sat * C_sat]
```

Changing `f_sat` must redistribute existing mass; it must not create or destroy
gas, acetate, DOM, or guild biomass. Handle the zero-area limits without
division by zero.

#### 9.3.1 Step 3 repartition and transport contract

Step 3 implements repartition and transport as independently testable kernels;
it does not dispatch the revised backend. For each partitioned acetate, guild,
or gas concentration, an area-transfer repartition moves the concentration of
the area that changed class into the receiving partition. If saturated area
increases from `f_old` to `f_new`, for example:

```text
C_sat,new = [f_old C_sat,old + (f_new - f_old) C_unsat,old] / f_new
C_unsat,new = C_unsat,old
```

The analogous expression mixes saturated concentration into the unsaturated
partition when saturated area decreases. This exactly preserves
`(1-f) C_unsat + f C_sat`, including transitions to and from `f=0` and `f=1`.
DOM remains the single authoritative bulk decomposition pool, so it needs no
partition remap; its Phase 2 vertical transport remains authoritative.

Vertical transport uses finite-volume, equal-and-opposite interface fluxes.
Acetate retains the explicit, closed-boundary kernel. Gas transport uses a
backward-Euler tridiagonal solve because gas-phase diffusion is too fast for an
ELM-length explicit timestep. For each layer the adapter supplies a transport
capacity `epsilon` and a diffusivity on a common gas-equivalent mobile-
concentration basis:

```text
C_bulk = epsilon * C_mobile
d(C_bulk)/dt = divergence(D_mobile * gradient(C_mobile))
```

In air-filled soil, `epsilon` includes gas-filled porosity plus Henry-weighted
liquid porosity and `D_mobile` is the effective gas diffusivity. In saturated
soil, `epsilon` is Henry-weighted liquid porosity and `D_mobile` is the aqueous
diffusivity multiplied by the same dimensionless Henry solubility. This change
of basis makes gas/aqueous interfaces continuous at phase equilibrium without
changing the authoritative state unit, which remains mol gas per m3 bulk soil.
The upper conductance includes the top half-layer, snow, ponded-water, and
boundary-layer resistances; the lower boundary is closed. The implicit solve is
nonnegative and exactly conservative to roundoff. A joint donor limiter then
bounds the explicit aerenchyma and ebullition losses against the post-diffusion
inventory.

Aerenchyma is represented as signed, first-order layer-to-atmosphere exchange.
This permits CH4, CO2, and H2 emission and O2 uptake through the same interface,
rather than suppressing reverse exchange. Ebullition removes a bounded fraction
of CH4 above a layer-specific solubility threshold. Its activation input
combines inundation, thaw state, and depth attenuation. Both kernels report a
positive-upward surface flux, and the layer-integrated tendency plus that flux
closes exactly.

The Step 4 ELM adapter constructs effective diffusivities, atmospheric boundary
concentrations, aerenchyma exchange rates, and ebullition activation from
temperature, atmospheric forcing, root state, and the named parameter-file
inputs. With `use_elm_microbe_methane_transport=.true.`, it additionally maps
current ELM hydrology into Henry capacities and phase-dependent diffusivities.
Keeping those model-state choices outside the pure transport kernel avoids old
global pointers and `HUM_HOL` branches. It also keeps the dimensional
interpretation of `m_dPlantTrans` inside the explicit adapter boundary: the
source equation requires `m s-1`, which becomes a first-order `s-1` exchange
after division by layer depth.

Two CLM-SPRUCE behaviors are intentionally not copied. Its inundation update
moves a fraction of each concentration without enforcing the area-weighted
inventory, and its vertical diffusion mutates adjacent layers sequentially
while mixing `1e-3` and `1e-4` metric factors. Step 3 instead uses the
conservative area-transfer equations and one simultaneous finite-volume flux
divergence. Generalized lateral exchange also multiplies the interface-area
flux by the smaller participating horizontal area fraction. The source's
hard-coded two-column equation has no corresponding geometry factor, but
omitting it in a weighted multi-topounit graph drives concentrations toward
infinity as a receiving saturated or unsaturated partition approaches zero
area. Golden-vector comparisons must classify these as accounting and
numerical repairs rather than port regressions.

An explicit parity test is available behind
`use_clm_microbe_dom_relaxation=.true.`. The executed source treats
`dom_diffus=1.8e-7` as an adjacent-layer relaxation rate in s-1, despite its
diffusivity metadata. The port applies that rate after reactions using a
closed-boundary backward-Euler finite-volume solve. The same transport matrix
is applied independently to DOM C, N, and P, so their column inventories are
conserved and a uniform C:N:P ratio remains uniform. This switch is off by
default and does not replace ELM's standard decomposition-pool transport. It
is a diagnostic for the large CLM-versus-ELM deep-DOM difference, particularly
because the ELM peat accumulation branch sets its base SOM diffusion to zero
and the normal 10x DOM multiplier therefore provides no diffusive mixing there.

For the scientific transport path, enable
`use_microbe_aqueous_transport=.true.`. This option is mutually exclusive with
the CLM relaxation diagnostic and is initially default-off while its transport
parameters and site behavior are validated. It replaces generic decomposition-
pool transport for dissolved pools only. DOM C, N, and organic P and the
adapter-owned acetate pool are transported; bacteria, fungi, litter, and SOM
remain attached to the soil matrix and continue to use their established ELM
behavior.

The aqueous operator is a conservative, backward-Euler finite-volume solve on
the actual ELM layer grid. It applies upwind advection using ELM's liquid-water
interface flux and centered diffusion/dispersion using porewater concentration:

```text
C_water = f_mobile C_bulk / theta_liq
J = q C_water - theta_liq D_eff d(C_water)/dz
D_eff = D_molecular S_liq^tortuosity f_temperature f_thaw
        + dispersivity |q| / theta_liq
```

Equivalently, the code combines `theta_liq*D_molecular` and
`dispersivity*|q|` into the conductivity multiplying the porewater gradient.
The same matrix is applied independently to DOM C, N, and P, preserving a
spatially uniform stoichiometric ratio while allowing an existing nonuniform
ratio to advect and diffuse conservatively. Clean infiltration carries zero
solute. Negative top soil-water flux associated with ground evaporation is
clipped to zero for solutes, because DOM and acetate cannot leave with water
vapor. The bottom interface is an open advective export when ELM diagnoses
downward water flow. Explicit coupling to later drainage and runoff terms is a
remaining hydrologic-integration task; neither is approximated as evaporation
or as an arbitrary relaxation sink. Carbon-isotope DOM transport is also not
yet implemented, so the namelist rejects this pathway with `use_c13` or
`use_c14` rather than allowing isotopic state to diverge silently.

The physical pathway deliberately differs from executed CLM-Microbe in four
ways: it uses metric grid spacing rather than equal-layer neighbor relaxation,
uses current liquid water and hydrologic flow, distinguishes molecular
diffusion from mechanical dispersion, and permits boundary leaching with
explicit C/N/P accounting. It publishes downward advective and diffusive DOM-C
interface fluxes, total aqueous C/N/P export, bottom DOM C/N/P export, and
per-timestep elemental residuals. Aqueous exports are included in column and
grid carbon balance, column N/P balance, and the monthly carbon budget.

The existing reconstructed `use_peatland_roots=.true.` control addresses a
separate coupling exposed by aqueous DOM-N transport. In peatland RD cases,
the standard uptake profile is the prescribed fine-root profile and does not
respond when transported DOM-N is mineralized below the surface root maximum.
The optional profile retains root presence as a hard constraint but weights
rooted, unsaturated layers by their current NH4 plus NO3 concentration:

```text
w_j = root_profile_j max(NH4_j + NO3_j, 0),  if theta_liq,j < porosity_j
w_j = 0,                                     otherwise
uptake_profile_j = w_j / sum_k(w_k dz_k)
```

The adaptive weighting is evaluated only on columns whose topounit has
positive peat depth and only for vascular PFTs. Uplands retain the ordinary
ELM uptake profile, while moss and other nonvascular PFTs retain their
prescribed shallow profiles. If no eligible mineral N exists on a peatland
column, the ordinary root profile is used. The
normalized profile only redistributes the existing column plant-N demand; the
standard RD competition code remains the sole plant-uptake flux and therefore
retains its N conservation and NH4/NO3 competition. This is intentionally a
general rooted-layer treatment rather than a hard-coded layer-6/7 rule. It is
part of the general peatland-root capability, defaults on with HUMHOL, and is
not owned by the methane module. It does not alter P uptake. A binary
below-porosity gate is used for the first experiment;
its sensitivity near saturation and the strong response to large deep-N
gradients require explicit validation before this becomes a recommended
configuration.

#### 9.3.2 Default Fickian and optional ELM gas-transport mappings

The revised backend defaults to a CLM-Microbe Fickian mapping and retains the
previous ELM multiphase mapping behind
`use_elm_microbe_methane_transport=.true.`. The following distinctions must be
retained in reviews and comparisons:

| Concern | Original CLM-Microbe | Default port | Optional ELM mapping |
| --- | --- | --- | --- |
| Prognostic gas state | Concentration semantics depend on the active equation and phase | Direct molar concentration-gradient basis with conservative area bookkeeping | Bulk-soil inventory mapped to one gas-equivalent mobile basis |
| Unsaturated diffusion | Aqueous `Fick_D_w*m_Fick_ad` in active layers | Same active coefficient and `T/298` response in both area partitions | ELM `CH4Mod` gas diffusivity scaled by air-filled porosity and soil structure |
| Saturated diffusion | Aqueous `Fick_D_w*m_Fick_ad` | Same active coefficient and temperature response | Aqueous diffusivity and Henry capacity on the common mobile basis |
| Atmosphere boundary | Timestep reset toward atmospheric Henry equilibrium | Finite top-half-layer Fickian conductance to Henry equilibrium | Boundary-layer, snow, pond-water, and top-half-soil resistance |
| Numerics | Sequential explicit layer updates with mixed metric factors | Conservative backward-Euler tridiagonal solve retaining the executed `1e-3` vertical coefficient | Same conservative backward-Euler solver |
| Lateral gas exchange with `use_humhol` | Hard-coded columns 1/2, 0.75/0.25 areas, matching layer indices, and `1e-4` coefficient | Arbitrary ELM `regional_target_ti` graph, actual weights/distances, shared horizontal-footprint scaling, absolute-elevation overlap, simultaneous donor limiting, and the source's factor-of-ten smaller lateral coefficient | Not applied; retains the earlier ELM-mapping behavior pending a separately validated multiphase lateral formulation |
| Reactions | Source microbial equations | Same ported equations | Same ported equations |

The default Fickian path also provides
`use_microbe_nonbog_lateral_gas_transport`. Its default value, true, uses every
configured edge and preserves the original generalized implementation. False
restricts exchange to edges for which both endpoint topounits are classified as
bog. This diagnostic control is useful for surfaces such as the three-topounit
SPRUCE configuration, where a bare non-bog fen/boardwalk column can otherwise
receive methane carbon from a productive hollow even though lateral dissolved
substrate transport has not yet been implemented.

The optional mapping is an intentional reuse of current ELM methane physics,
not a literal copy of all `CH4Mod` state or solver code. Above/below-water-table
phase selection and its existing standard-parameter controls are applied inside
the revised backend only when the option is true. Legacy `CH4Mod` remains
independently selectable and is not
modified by this path. Consequently, coupled results made with the earlier
aqueous-only Phase 3 adapter are useful reaction baselines but are not
scientifically comparable surface-transport baselines.

### 9.4 Units and conservation

Define units at module boundaries and convert only there:

- ELM decomp pools: `g C|N|P m-3 soil`;
- functional biomass and acetate: one documented carbon unit, preferably
  `g C m-3 soil`;
- gas state: `mol gas m-3 bulk soil` (inventory density); phase-equilibrated
  mobile concentration is an internal transport-solver quantity;
- internal reaction rates: matching state units per second;
- column surface CH4: `kg C m-2 s-1`, positive to atmosphere, for
  `flux_ch4_grc`; and
- net methane correction: `g C m-2 s-1` for `nem_grc`.

Central conversion functions must use ELM constants for atomic weights and
layer thickness. Every one-layer closed-box test must satisfy carbon atoms,
nitrogen, phosphorus, hydrogen, and oxygen to the tolerance appropriate for the
solver. The atmosphere boundary flux and any prescribed-atmosphere source are
included in those budgets.

### 9.5 Time integration

Reactions reproduce the CLM-Microbe explicit ordering through visible staged
tendencies. If a reaction can consume more substrate than available over an ELM
timestep, all competing consumers are scaled proportionally. The reaction
kernel remains explicit so one-layer reaction parity is directly testable.

Gas diffusion is the documented exception: the multiphase correction uses a
backward-Euler tridiagonal solve because gas-phase diffusivity makes an explicit
ELM timestep unstable. Aerenchyma and ebullition remain explicit and bounded.
Timestep-convergence tests must distinguish reaction sensitivity from the
unconditionally stable diffusion update; bounded reaction subcycling remains a
possible later change.

### 9.6 Step 4 transaction and ELM accounting contract

Step 4 is the boundary between the pure reaction/transport kernels and mutable
ELM state. It stages a complete candidate update before changing any ELM-owned
array. For each soil layer the transaction:

1. presents both saturated and unsaturated reaction kernels with the same
   authoritative bulk DOM concentration;
2. restricts reactive DOM C to the amount supported by DOM N and P;
3. area-weights the two reaction tendencies once, with mortality limited
   proportionally when its return to DOM would require more mineral N or P than
   is available;
4. produces paired DOM and mineral N/P tendencies and candidate methane state;
5. applies acetate and gas transport to candidate state; and
6. commits all candidate arrays together only after nonnegativity and C/N/P
   residual checks pass.

The transaction kernel is pure and independent of ELM data types. The ELM
adapter remains responsible for translating current hydrology, temperature,
pH, root distribution, atmospheric forcing, layer geometry, DOM C/N/P, and
mineral N/P into kernel units. It is also responsible for constructing the
transport boundary terms described in section 9.3.1. DOM continues to use the
Phase 2 decomposition-cascade transport path; Step 4 transports acetate and
the four gases and must not transport DOM a second time.

Carbon accounting uses one explicit ledger. DOM is excluded from the revised
methane storage addition because it is already included in ELM's standard
decomposition pools. The additional column storage is the area-weighted,
layer-integrated sum of acetate, the four functional-guild biomasses, bulk-soil
CH4-C, and bulk-soil CO2-C inventories. The matching external loss is the sum of
the net positive-upward CH4-C and CO2-C surface fluxes. Internal reaction,
partition exchange, and vertical transport do not enter the external budget.
The column, gridcell, and monthly C-budget interfaces will receive these
additional storage and flux terms only when `use_microbe_methane=.true.`;
disabled-mode calls and arithmetic remain unchanged.

The atmosphere-facing mapping is likewise explicit. Net CH4 surface exchange
is converted to `kg C m-2 s-1` for `flux_ch4_grc`. Net CO2 surface exchange is
converted to `g C m-2 s-1` and supplies the revised backend's methane-related
NEE correction. The adapter must not reuse legacy `CH4Mod`'s production-minus-
oxidation bookkeeping because those processes are internal transfers in the
revised carbon ledger. The legacy offline/online atmosphere policy remains an
outer ELM coupling decision rather than a reaction-kernel concern.

The implemented adapter uses `max(fsat_col, frac_h2osfc)` as the current
saturated-area fraction and conservatively remaps every partitioned state
before reaction. It takes layer temperature and geometry from the column
state, liquid/ice content and surface water from the column water state,
porosity, suction, water potential, and root fraction from `soilstate_vars`,
and atmospheric partial pressures from `atm2lnd_vars`. Native ELM currently
allocates `chemstate_vars%soil_pH` without populating it, while the external
chemistry modes that do populate it are unsupported by this option. The
adapter therefore uses a separate named soil-pH fallback as an explicit
temporary base pH; acetate feedback can still modify effective pH inside the
reaction kernel. The reference fallback is 7 and the initial SPRUCE test uses
4.5 without changing the response optimum (`ph_opt = 7`). Phase 4 must add a
spatially resolved native-ELM soil-pH input before
pH-response calibration. The named atmospheric mixing ratios are used only as
missing-forcing fallbacks. The source water-potential response is evaluated in
ELM's MPa units and combined with liquid saturation for the unsaturated
reaction scalar.

Because the legacy solver normally constructs `rootfr_col` internally, the
revised adapter instead aggregates the authoritative patch root profile to
columns before computing plant transport. Missing roots on unvegetated
columns become zero transport. In the default mapping, signed surface diffusion
relaxes the stored concentration toward atmospheric Henry equilibrium through
a top-half-layer Fickian conductance. In the optional ELM mapping, it relaxes
toward atmospheric concentration on the common gas-equivalent basis and Henry
capacity establishes the corresponding aqueous inventory. Plant exchange
for O2 and CO2 is also signed and bounded by the amount needed to reach its
temperature-dependent Henry equilibrium over one timestep. CH4 and H2
aerenchyma transport is a one-way emission path:
it is zero below the larger of atmospheric equilibrium and its named emission
threshold and cannot import atmospheric gas. CH4 ebullition retains the same
CH4 threshold. Keeping the emission threshold separate from the atmospheric
boundary prevents a cold-start soil from receiving an artificial CH4 source
while retaining exact surface-flux closure.

In the optional ELM mapping, gas surface conductance combines the patch-to-
column boundary conductance with top-half-layer, snow, and ponded-water
resistances. The first three gases reuse ELM's `d_con_w`, `d_con_g`, `c_h_inv`,
`kh_theta`, and `kh_tbase` constants. Gas-filled soil uses the same
organic/mineral Millington-Quirk/Moldrup blend as ELM `CH4Mod`; saturated soil
retains the CLM-Microbe `m_Fick_ad` multiplier and temperature exponent.
Aqueous diffusivity and capacity are both transformed by dimensionless Henry
solubility so the transport solve has one mobile basis.
Because those shared tables stop at CO2, H2 uses the corresponding non-tunable
CLM-SPRUCE physical constants (`4.5e-9 m2 s-1` in water, `6.11e-5 m2 s-1` in
air near 298 K, `1282.1 L atm mol-1`, and a `500 K` Henry temperature
coefficient). These are physical conversion constants, not additional
calibration parameters. Dimensional analysis of the source plant-flux equation
resolves `m_dPlantTrans` as `m s-1`; division by finite layer depth produces
the `s-1` exchange rate expected by the transport kernel.

The adapter publishes pathway-resolved CH4 surface exchange as
`MM_CH4_SURF_DIFF`, `MM_CH4_SURF_AERE`, and `MM_CH4_SURF_EBUL`, plus `_UNSAT`
and `_SAT` area-weighted contributions for each. The paired contributions must
sum to the bulk pathway field, and all three bulk pathways must sum to the net
CH4 surface exchange after unit conversion. These fields are required for the
upland atmospheric-uptake audit; net `FCH4` alone cannot distinguish a diffusion
defect from local production, oxidation, plant transport, or ebullition.

The adapter writes NH4 as the mineral-N counterpool and immediately refreshes
total mineral N as `NH4 + NO3`; solution P is the mineral-P counterpool. It
publishes additional methane-system storage, signed CH4+CO2 surface carbon
exchange, CH4-C exchange in the established atmosphere units, and the CO2 NEE
correction. Optional column, grid, and monthly carbon-budget arguments include
these terms only for the revised backend. Existing calls omit the arguments and
therefore preserve disabled-mode arithmetic.

The native RD phosphorus update can leave `solutionp_vr` slightly negative
when the solution pool is depleted. The revised transaction carries an
inherited negative value through its P ledger instead of clipping it and
creating phosphorus. Its mortality limiter treats that value as zero available
P, and the transaction may leave the inherited deficit unchanged or improve
it, but may not make it more negative. DOM-P remains subject to the strict
nonnegative-state check.

Nitrification participates in the revised oxygen budget. Before allocation,
its potential flux is calculated by the standard ELM nitrogen code; revised
methane no longer gives nitrification an independent cap against the complete
O2 inventory. After ELM resolves decomposition and NH4 competition, the
adapter reconstructs the CH4Mod-style potential demand from heterotrophic
respiration, root respiration, and `pot_f_nit_vr` (the latter at 2 mol O2 per
mol N). The reaction kernel adds aerobic acetate and methane oxidation and
applies one common stress when the total exceeds available O2. The accepted
non-microbial and microbial demands are removed together, so no process has
first access to the finite inventory and neither partition becomes negative.
The adapter then publishes revised O2 concentration, the common stress,
aerobic demand, and saturated fraction through the existing `ch4_vars`
oxygen/anoxia interface for the next standard decomposition and
nitrification/denitrification calculation. In revised mode that object is a
compatibility carrier; the legacy CH4 solver is not called. As in legacy
`CH4Mod`, nitrification's N flux is lag-coupled by ELM's existing call order;
it is not retroactively changed after nutrient allocation.

The root-respiration term must not be read from `col_cf%rr_vr` in revised
mode. That array is a legacy CH4Mod work product and remains at its fill value
when the legacy solver is bypassed. The implemented adapter instead aggregates
`veg_cf%rr * rootfr_patch * veg_pp%wtcol` over the soil-patch filter, matching
the legacy profile construction without dispatching legacy methane. This
prevents a fictitious near-infinite O2 demand and preserves the mutually
exclusive backend design.

## 10. Parameters and input data

All parameters used by the enabled feature will be migrated into ELM's standard
parameter NetCDF. There will not be a feature-specific parameter file. The
existing `paramfile` selection, staging, provenance, PFT dimension, and parallel
read machinery will be used. The new variables are conditionally required and
read only when `use_microbe_methane=.true.`, preserving compatibility with
existing parameter files and disabled-mode B4B behavior.

The tables below are the source-level inventory. They use the shared
CLM-SPRUCE/CLM-Microbe identifiers so that every equation can be traced during
the port. An
implementation may adopt clearer ELM-style NetCDF names, but each renamed
variable must carry a `legacy_name` attribute. Each new variable also needs
`units`, `long_name`, valid-range, source-revision, and provenance metadata.
The Phase 3 integration default uses the CLM-Microbe runtime file at commit
`9c2e0a048bb3799669d32b91e6cc76efb36d4b75` wherever it supplies an active
parameter. Runtime-file values and declaration defaults are not treated as
interchangeable. The runtime numbers are not copied blindly: inspection of the
active source shows that stored carbon is divided by 12 into mol C m-3 and that
reaction tendencies are multiplied by a timestep in seconds. The Phase 3
schema therefore applies these executable-equivalent conversions:

| Legacy quantity used by active code | ELM parameter units | Conversion |
| --- | --- | ---: |
| Concentration or half-saturation, mol m-3 | mmol m-3 | multiply by 1000 |
| Volumetric rate, mol m-3 s-1 | mmol m-3 d-1 | multiply by 86,400,000 |
| Specific growth, mortality, or oxidation rate, s-1 | d-1 | multiply by 86,400 |
| Functional biomass floor, mol C m-3 | g C m-3 | multiply by `catomw = 12.011` |

The acetate acidification coefficient is divided by 1000 because the ELM
equation supplies acetate in mmol C m-3 rather than mol C m-3. The same
concentration and rate conversions are applied to active source literals such
as the acetate feedback scale, gas-inhibition scales, H2 transport threshold,
and aerobic acetate oxidation coefficient. Each converted value carries a
provenance suffix in the standard parameter NetCDF. Production suitability
remains subject to the audit in section 10.6.

### 10.1 Microbial decomposition parameters

These are the parameters used by the 11-pool/47-transition cascade. Unless
noted otherwise, the source variables are PFT-indexed and will use the standard
ELM PFT dimension.

| Source identifier | Scope | Use in the enabled cascade |
| --- | --- | --- |
| `decomp_depth_efolding` | PFT | Depth e-folding scale for decomposition; reuse the existing standard ELM variable rather than duplicate it |
| `k_dom` | PFT | DOM daily turnover probability/rate |
| `k_bacteria` | PFT | Bacterial biomass daily turnover probability/rate |
| `k_fungi` | PFT | Fungal biomass daily turnover probability/rate |
| `m_rf_s1m`, `m_rf_s2m`, `m_rf_s3m`, `m_rf_s4m` | PFT | SOM1-SOM4 microbial carbon-use efficiencies (the retained fraction after pathway respiration) |
| `m_batm_f`, `m_fatm_f` | PFT | Bacterial and fungal turnover fractions respired to atmosphere |
| `m_bdom_f`, `m_fdom_f` | PFT | Bacterial and fungal lysis fractions routed to DOM |
| `m_bs1_f`, `m_bs2_f`, `m_bs3_f` | PFT | Bacterial residue fractions routed to SOM1-SOM3 |
| `m_fs1_f`, `m_fs2_f`, `m_fs3_f` | PFT | Fungal residue fractions routed to SOM1-SOM3 |
| `m_domb_f`, `m_domf_f` | PFT | DOM uptake fractions routed to bacteria and fungi |
| `m_doms1_f`, `m_doms2_f`, `m_doms3_f` | PFT | DOM stabilization fractions routed to SOM1-SOM3 |
| `cn_bacteria`, `cn_fungi` | PFT | Fixed bacterial and fungal biomass C:N ratios |
| `cn_dom` | Scalar | DOM C:N ratio; hard-coded as 10 in CLM-SPRUCE |
| `bacteria_pool_cn`, `fungi_pool_cn` | Scalar, new names | Initial/authoritative pool C:N values represented by literals 5 and 15 in CLM-SPRUCE; reconcile with `cn_bacteria`/`cn_fungi` before freezing the schema |
| `CUEmax` | Scalar | Maximum litter-to-microbial carbon-use efficiency; hard-coded as 0.8 |
| `microbe_cue_cn_target` | Scalar, new name | C:N target in the litter CUE equation; literal 8.0 in CLM-SPRUCE |
| `microbe_allocation_cn_exponent` | Scalar, new name | Bacteria/fungi allocation exponent; literal 0.6 in CLM-SPRUCE |
| `bacteria_initial_c`, `fungi_initial_c`, `dom_initial_c` | Scalar, new names | Cold-start C stocks; source literals are `1e-5`, `1e-5`, and 0 `g C m-3 soil` |
| `cp_bacteria`, `cp_fungi`, `cp_dom` | PFT or scalar, new names | Required C:P ratios for the target ELM CNP model; absent from CLM-SPRUCE and requiring science-owner values |
| `dom_som_diffusion_multiplier` | Scalar, new name | Multiplier on SOM-solver diffusivity for DOM; literal 10 in `CNSoilLittVertTranspMod.F90` |
| `microbe_som2_q10`, `microbe_som3_q10`, `microbe_som4_q10` | Scalar, new names | Pool-specific temperature responses; CLM-SPRUCE literals 1.5, 2.0, and 2.5 |
| `microbe_dom_q10` | Scalar, new name | DOM temperature response; CLM-SPRUCE literal 1.25 |

These four Q10 parameters use ELM's existing 25 degrees C reference and 10 K
temperature interval conventions; those shared mathematical constants are not
duplicated in the microbial parameter namespace. Litter, SOM1, bacteria, and
fungi continue to use ELM's existing `Q10_hr` and `froz_q10` parameters.

The following active cascade fractions are literals in CLM-SPRUCE and must be
promoted into the standard ELM parameter file rather than copied as literals:

| New parameter | CLM-SPRUCE value | Use |
| --- | ---: | --- |
| `l1dom_f`, `l2dom_f`, `l3dom_f` | 0.10, 0.08, 0.06 | Litter1-Litter3 solubilization to DOM |
| `s1dom_f`, `s2dom_f`, `s3dom_f`, `s4dom_f` | 0.18, 0.14, 0.10, 0.06 | SOM1-SOM4 solubilization to DOM |
| `l1s1_f`, `l2s2_f`, `l3s3_f` | 0.19, 0.21, 0.23 | Direct litter stabilization to corresponding SOM pools |
| `s1s2_f`, `s2s3_f`, `s3s4_f` | 0.14, 0.23, 0.27 | Direct SOM stabilization to the next pool |

`bs4_f`, `fs4_f`, and `doms4_f` are not independent inputs. They are residuals
computed respectively from the bacterial, fungal, and DOM path fractions. The
reader must reject a negative residual; the source's `max(0, residual)` behavior
silently loses the donor-fraction closure and will not be reproduced.

The Phase 2 test manifest records provenance per variable. Turnover and routing
values copied from `microbepar_in` remain explicitly *unverified* because the
legacy positional reader did not populate the correspondingly named PFT arrays.
The bacterial and fungal PFT C:N test values repeat the source pool literals,
while all three C:P values are new CNP test assumptions. None of those groups is
a production calibration until an archived successful-run parameter file or a
science-owner decision resolves it.

### 10.2 Revised methane reaction and isotope parameters

The following source variables are referenced by active reaction equations and
will be migrated as named scalar variables in the standard ELM parameter file.
The role descriptions are intentionally equation-oriented; final names and
units will be frozen only after the archived-run audit.

| Source identifiers | Role |
| --- | --- |
| `MFGbiomin` | Minimum biomass floor for each of the four methane functional guilds |
| `m_dKAce`, `m_dAceProdACmax`, `m_dKAceProdO2` | DOM/acetate production and acetate/O2 kinetic controls |
| `m_dH2ProdAcemax`, `m_dKH2ProdAce`, `m_dKCO2ProdAce`, `m_dH2AceProdQ10` | Acetogenesis maximum rate, H2 and CO2 half-saturation constants, and temperature response |
| `m_dGrowRH2Methanogens`, `m_dDeadRH2Methanogens`, `m_dYH2Methanogens` | Hydrogenotrophic methanogen growth, mortality, and yield |
| `m_dGrowRAceMethanogens`, `m_dDeadRAceMethanogens`, `m_dYAceMethanogens` | Acetoclastic methanogen growth, mortality, and yield |
| `m_dGrowRMethanotrophs`, `m_dDeadRMethanotrophs`, `m_dYMethanotrophs` | Aerobic methanotroph growth, mortality, and yield |
| `m_dGrowRAOMMethanotrophs`, `m_dDeadRAOMMethanotrophs`, `m_dYAOMMethanotrophs` | Anaerobic methanotroph growth, mortality, and yield |
| `m_dACMinQ10` | Temperature response for DOM-to-acetate mineralization |
| `m_dKH2ProdCH4`, `m_dKCO2ProdCH4`, `m_dH2CH4ProdQ10` | H2 and CO2 half-saturation constants and temperature response for hydrogenotrophic methanogenesis |
| `m_dKCH4ProdAce`, `m_dCH4ProdQ10`, `m_drCH4Prod` | Acetate half-saturation, temperature response, and CH4 yield for acetoclastic methanogenesis |
| `m_dKCH4OxidCH4`, `m_dKCH4OxidO2`, `m_dCH4OxidQ10`, `m_drCH4Oxid` | CH4 and O2 half-saturation constants, temperature response, and O2:CH4 ratio for aerobic oxidation |
| `m_dKAOMCH4OxidCH4`, `m_dAOMCH4OxidQ10` | CH4 half-saturation constant and temperature response for anaerobic oxidation |
| `m_drAer` | O2:C stoichiometric ratio used by aerobic decomposition |
| `m_dCH4min` | CH4 threshold used by ebullition and plant transport |
| `pHmin`, `pHmax`, `pHopt` | Minimum, maximum, and optimum pH response points |
| `frac_ace`, `frac_acch4`, `frac_hych4`, `frac_acetogenesis`, `frac_ch4ox`, `frac_ch4aom` | C13 fractionation for acetate production, acetoclastic and hydrogenotrophic methanogenesis, acetogenesis, aerobic oxidation, and anaerobic oxidation |

`frac_doc` is declared and read in CLM-SPRUCE, but its only inspected use is
commented out. It is retained in the audit table as a legacy candidate, not
treated as an active parameter. C13 parameters become required only when the
C13 implementation gate in section 12 is enabled.

### 10.3 Revised methane transport and boundary parameters

| Source identifier | Use and target treatment |
| --- | --- |
| `dom_diffus` | DOM and acetate liquid diffusion coefficient; migrate |
| `m_Fick_ad` | Multiplier on aqueous gas diffusion; migrate |
| `m_dPlantTrans` | Root/aerenchyma transport coefficient; migrate as m s-1 |
| `g_dMaxH2inWater` | Dissolved-H2 threshold for plant transport; migrate |
| `atmch4`, `atmo2`, `atmco2`, `atmh2` | Fallback atmospheric mixing ratios; migrate, but current ELM atmosphere forcing takes precedence when supplied |
| `Fick_D_w(1:4)` | Shared water diffusion coefficients for CH4, O2, CO2, and H2; reuse the current ELM equivalent rather than duplicate |
| `Henry_C_w(1:4)`, `Henry_kHpc_w(1:4)`, `kh_tbase` | Shared Henry-law coefficients and reference temperature; reuse current ELM equivalents |
| `rgasLatm` | Gas constant; use the current ELM physical constant, not a parameter-file duplicate |

The multiphase ELM adapter intentionally reuses the following parameters that
are already present in the standard ELM parameter file. They are scientific
dependencies of revised methane transport and must be included in parameter
audits, but must not be duplicated under `microbe_methane_*` names.

| Existing ELM parameter | Current test value | Revised-backend use |
| --- | ---: | --- |
| `f_sat` | 0.95 | Water-filled-porosity threshold selecting aqueous transport in an unsaturated-area layer |
| `satpow` | 2 | Porosity exponent for saturated aqueous diffusivity |
| `scale_factor_gasdiff` | 1 | Multiplier on effective air-phase gas diffusivity |
| `scale_factor_liqdiff` | 1 | Multiplier on effective aqueous diffusivity |
| `organic_max` | 130 kg m-3 | Organic-matter scale blending peat and mineral gas-diffusion structure |

The standard soil hydraulic exponent `bsw`, layer porosity, soil organic matter,
water/ice content, and layer geometry are spatial state/parameter inputs to the
same calculation. Molecular diffusion tables, Henry coefficients, the gas
constant, water/ice densities, and atomic weights are shared physical constants,
not tunable module parameters. All CLM-Microbe-specific scientific controls in
this section are migrated into the standard ELM parameter file; the table above
documents the additional standard-ELM parameters now consumed by the adapter.

The inspected equations also contain the following active science coefficients
as literals. They will either be promoted to clearly named variables in the
standard parameter file or, where indicated, replaced by a current ELM state or
unit conversion. This decision is made explicitly during the equation audit,
not by leaving unexplained literals in the port.

| Proposed name | Source value | Treatment |
| --- | ---: | --- |
| `saturation_reaction_threshold` | 0.99 | Map to current ELM saturated fraction/water-table state if equivalent; otherwise parameterize |
| `reaction_t_ref` | 286.65 K | Parameterize common Q10 reference temperature |
| `q10_temperature_interval` | 10 K | Named shared mathematical constant |
| `aom_t_ref` | 13.5 in a soil-temperature expression | Resolve the apparent Celsius/Kelvin inconsistency before parameterizing |
| `acetate_ph_trigger` | 5.5 | Parameterize threshold for acetate-induced pH adjustment |
| `acetate_acidification_coefficient` | `0.0042e-6` | Parameterize pH adjustment coefficient with explicit units |
| `acetate_feedback_half_saturation` | 0.1 | Parameterize acetate feedback term |
| `soil_water_potential_min` | -10 | Reuse an ELM moisture scalar if scientifically equivalent; otherwise parameterize with units |
| `soil_suction_unit_conversion` | `-9.8e-6` | Implement as a documented units conversion, not a tunable parameter |
| `hydrogenotrophic_co2_inhibition_scale` | 9.2 | Parameterize CO2 limitation scale after verifying the source comment/equation |
| `aom_o2_inhibition_scale` | 4.6 | Parameterize O2 inhibition scale for AOM |
| `aerobic_acetate_oxidation_rate` | 0.05 | Parameterize coefficient in the active alternate unsaturated branch |
| `ch4_h2_root_efolding_depth` | 0.25 m | Parameterize CH4/H2 root-transport depth scale |
| `o2_root_efolding_depth` | 0.15 m | Parameterize O2 root-transport depth scale |
| `ebullition_efolding_depth` | 0.35 m | Parameterize ebullition depth scale |
| `transport_thaw_threshold` | -0.1 in source temperature units | Map to current ELM frozen-state logic or parameterize after resolving units |
| `plant_o2_consumption_fraction` | 0.001 | Parameterize plant O2 consumption multiplier |
| `plant_co2_flux_fraction` | 0.001 | Parameterize plant CO2 flux multiplier |
| `aqueous_diffusion_t_ref` | 298 K | Reuse shared diffusion reference temperature if available; otherwise parameterize |
| `aqueous_diffusion_t_exponent` | 1.87 | Parameterize temperature exponent |

Cold-start acetate and all four bulk-soil gas inventory densities are zero in the
source. Each functional guild starts at `MFGbiomin`. These are documented
initialization rules rather than additional independent parameters; if science
owners require nonzero configurable initial concentrations, the corresponding
`acetate_initial_c`, `ch4_initial_concentration`, `o2_initial_concentration`,
`co2_initial_concentration`, and `h2_initial_concentration` variables will also
be added to the standard ELM parameter file. Restart values always take
precedence.

The source also substitutes hard-coded atmospheric concentrations
`8.2e-5`, `10.16`, `0.019`, and `2.7e-5 mol m-3` when atmospheric temperature
is zero. The port will reject that invalid forcing rather than promote these
fallbacks as science parameters.

Reaction stoichiometry such as 2 O2 per CH4, the four-H2 hydrogenotrophic
reaction, atomic-weight conversions, seconds-per-day, and exact metric unit
factors are documented constants, not calibration parameters. Site-specific
`HUM_HOL` layer geometry, fixed area fractions, and hard-coded SPRUCE grid
indices are deliberately excluded. The generalized gas operator uses the
standard ELM topounit graph, weights, elevations, and lateral distances. The
source's executed factor-of-ten distinction between vertical (`1e-3`) and
lateral (`1e-4`) gas coefficients is retained for parity and must be
scientifically validated.

#### 10.3.1 Phase 3 Step 1 standard parameter-file names

Step 1 freezes 64 active bulk-model inputs under a
`microbe_methane_` namespace and reads them through ELM's existing `paramfile`
path only when `use_microbe_methane=.true.`. This avoids collisions with the
legacy `CH4Mod` variables and allows both backends to remain supported. The
exact inventory is:

- initialization and acetate production:
  `microbe_methane_mfg_biomass_min`, `microbe_methane_k_acetate`,
  `microbe_methane_acetate_prod_max`,
  `microbe_methane_k_acetate_prod_o2`,
  `microbe_methane_dom_to_acetate_q10`;
- acetogenesis:
  `microbe_methane_acetogenesis_max`,
  `microbe_methane_k_acetogenesis_h2`,
  `microbe_methane_k_acetogenesis_co2`,
  `microbe_methane_acetogenesis_q10`;
- hydrogenotrophic methanogens:
  `microbe_methane_h2_methanogen_growth_rate`,
  `microbe_methane_h2_methanogen_death_rate`,
  `microbe_methane_h2_methanogen_yield`,
  `microbe_methane_k_h2_methanogenesis_h2`,
  `microbe_methane_k_h2_methanogenesis_co2`,
  `microbe_methane_h2_methanogenesis_q10`,
  `microbe_methane_h2_methanogenesis_co2_inhibition_scale`;
- acetoclastic methanogens:
  `microbe_methane_acetate_methanogen_growth_rate`,
  `microbe_methane_acetate_methanogen_death_rate`,
  `microbe_methane_acetate_methanogen_yield`,
  `microbe_methane_k_acetoclastic_methanogenesis_acetate`,
  `microbe_methane_acetoclastic_methanogenesis_q10`,
  `microbe_methane_acetoclastic_methanogenesis_ch4_yield`;
- aerobic methanotrophs and oxidation:
  `microbe_methane_aerobic_methanotroph_growth_rate`,
  `microbe_methane_aerobic_methanotroph_death_rate`,
  `microbe_methane_aerobic_methanotroph_yield`,
  `microbe_methane_k_aerobic_oxidation_ch4`,
  `microbe_methane_k_aerobic_oxidation_o2`,
  `microbe_methane_aerobic_oxidation_q10`,
  `microbe_methane_aerobic_oxidation_o2_ch4_ratio`,
  `microbe_methane_aerobic_decomp_o2_c_ratio`,
  `microbe_methane_aerobic_acetate_oxidation_rate`;
- anaerobic methanotrophs and oxidation:
  `microbe_methane_anaerobic_methanotroph_growth_rate`,
  `microbe_methane_anaerobic_methanotroph_death_rate`,
  `microbe_methane_anaerobic_methanotroph_yield`,
  `microbe_methane_k_anaerobic_oxidation_ch4`,
  `microbe_methane_anaerobic_oxidation_q10`,
  `microbe_methane_aom_o2_inhibition_scale`;
- pH, saturation, temperature, and moisture controls:
  `microbe_methane_ph_min`, `microbe_methane_ph_max`,
  `microbe_methane_ph_opt`, `microbe_methane_acetate_ph_trigger`,
  `microbe_methane_acidification_coefficient`,
  `microbe_methane_acetate_feedback_half_saturation`,
  `microbe_methane_soil_water_potential_min`,
  `microbe_methane_saturation_reaction_threshold`,
  `microbe_methane_reaction_t_ref`, `microbe_methane_aom_t_ref`;
- diffusion, plant transport, and ebullition:
  `microbe_methane_dom_diffusivity`,
  `microbe_methane_dom_relaxation_rate`,
  `microbe_methane_aqueous_dom_molecular_diffusivity`,
  `microbe_methane_aqueous_acetate_molecular_diffusivity`,
  `microbe_methane_aqueous_dom_mobile_fraction`,
  `microbe_methane_aqueous_acetate_mobile_fraction`,
  `microbe_methane_aqueous_solute_dispersivity`,
  `microbe_methane_aqueous_solute_tortuosity_exponent`,
  `microbe_methane_aqueous_solute_min_liquid_fraction`,
  `microbe_methane_aqueous_gas_diffusion_multiplier`,
  `microbe_methane_plant_transport_coefficient`,
  `microbe_methane_h2_plant_transport_threshold`,
  `microbe_methane_ch4_transport_threshold`,
  `microbe_methane_ch4_h2_root_efold_depth`,
  `microbe_methane_o2_root_efold_depth`,
  `microbe_methane_ebullition_efold_depth`,
  `microbe_methane_transport_thaw_threshold`,
  `microbe_methane_plant_o2_consumption_fraction`,
  `microbe_methane_plant_co2_flux_fraction`,
  `microbe_methane_aqueous_diffusion_t_ref`,
  `microbe_methane_aqueous_diffusion_temperature_exponent`;
- atmospheric fallback boundary values:
  `microbe_methane_atmospheric_ch4_mixing_ratio`,
  `microbe_methane_atmospheric_o2_mixing_ratio`,
  `microbe_methane_atmospheric_co2_mixing_ratio`, and
  `microbe_methane_atmospheric_h2_mixing_ratio`.

The versioned reference manifest records the legacy identifier, value, units,
and provenance of every name. It is an injection and integration-test fixture,
not a production calibration. Manifest schema
`elm_microbe_methane_phase3_aqueous_transport_v1` contains 73 parameters. It
uses 38 values copied from the linked CLM-Microbe runtime file; 28 more are
declarations, promoted active literals, or the two explicit unit corrections
in the same inspected source tree. The seven aqueous-solute parameters are
new scientific hypotheses and numerical controls with explicit validation
provenance, not values inferred from CLM-Microbe's relaxation equation. The
runtime baseline differs materially from the original Step 1
declaration-based fixture: for example, AOM net unlimited growth decreases from
`0.022 d-1` to `0.002 d-1`, the aerobic decomposition O2:C coefficient changes
from 2 to 0.002, and `dom_diffus` changes to `1.8e-7`. These are intentional
baseline changes, not new calibrations, and require sensitivity testing.

The active concentration half-saturation use sites have now been audited: the
legacy arithmetic uses mol m-3 (numerically mmol L-1), so the ELM mmol m-3
values require the documented factor of 1,000. This establishes execution
parity, including `m_dKAce=16` mapping to 16,000 mmol C m-3, but it does not
establish that the source comments or values are scientifically appropriate.
Kinetic rate intent and transport units still need the archived-run audit. In
particular, the linked DOM diffusion equation is not dimensionally equivalent
to the port's finite-volume operator, and `m_drAer=0.002` is scientifically
suspicious for a quantity represented in the port as mol O2 per mol C. The
runtime value is retained for reproducibility and flagged for validation rather
than silently replaced. Step 1 also resolves two unambiguous source defects:
`aom_t_ref=286.65 K` replaces a 13.5-Celsius literal used against Kelvin state,
and `transport_thaw_threshold=273.05 K` replaces a -0.1-Celsius literal used
against Kelvin state. C13 fractionation inputs remain deferred to Phase 5.

### 10.4 Parameters already owned by ELM

The enabled modules also consume ordinary ELM parameter/state inputs such as
soil pH, porosity and suction parameters, litter and SOM decay rates, CWD
cellulose/lignin fractions, `ksomfac`, root fraction, PFT identity, and layer
geometry. Those remain owned by their current ELM structures and standard
parameter variables. They must be recorded as interface dependencies in the
implementation, but must not be copied into a second microbial namespace.

Similarly, `ch4offline`, `allowlakeprod`, history controls,
`use_microbe_methane`, `use_elm_microbe_methane_transport`, and
`use_peatland_roots` are runtime controls, not scientific parameter-file
variables.

The inspected CLM-SPRUCE branch also hard-codes a 1.35 above-freezing base Q10,
a 0.3 moisture-response floor, and a 0.5 moisture exponent outside its
`MICROBE` conditional. Those are branch-wide SPRUCE tuning changes, not
microbial-module parameters, and are not imported by this option. Current ELM
temperature and moisture controls remain authoritative for litter, SOM1,
bacteria, and fungi. The SOM2-SOM4 and DOM Q10 literals that occur inside the
`MICROBE` path are migrated explicitly in section 10.1.

### 10.5 Declared or supplied legacy values not active in the inspected path

For completeness, the following names occur in `microbevarcon.F90` or one of
the inspected 111-record CLM-SPRUCE and 89-record CLM-Microbe parameter files,
but are not referenced by active calculations in the inspected
revised-methane/microbial-cascade path. They will be reported by the conversion
audit and will not be added to the standard ELM parameter file unless an
archived reference executable demonstrates an active use:

- older methane formulation: `q10ch4base`, `q10ch4`, `vmax_ch4_oxid`, `k_m`,
  `q10_ch4oxid`, `smp_crit`, `aereoxid`, `mino2lim`, `rootlitfrac`,
  `scale_factor_aere`, `vgc_max`, `organic_max`, `satpow`, `cnscalefactor`,
  `f_ch4`, `k_m_o2`, `nongrassporosratio`, `usephfact`, `k_m_unsat`,
  `vmax_oxid_unsat`, `scale_factor_gasdiff`, `scale_factor_liqdiff`, `redoxlag`,
  `usefrootc`, `redoxlag_vertical`, and `fin_use_fsat`;
- hard-wired/unused logicals: `transpirationloss`, `ch4rmcnlim`,
  `anoxicmicrosites`, and `ch4frzout`;
- unused DOM production controls: `plant2doc`/text-file `lit2avc` and
  `som2doc`/text-file `som2avc`;
- imported but inactive microbial decomposition constants: `micbiocn`,
  `micbioMR`, `dock`, `CUEref`, `CUEt`, `Tcueref`, `Tsref`, `Tmref`, `Tsmin`,
  `Tmmin`, `Msmin`, `Msmax`, and `Mmmin`;
- declared/read reaction values whose equations are commented or inactive:
  `m_dAceProdQ10`, `m_dACProdQ10`, `m_dAceH2min`, `m_dCH4H2min`,
  `m_dKCH4ProdO2`, `m_dKAerO2`, `m_dAerDecomQ10`, `m_dKe`, `m_dAirCH4`,
  `m_dAirH2`, `m_dAirO2`, `m_dAirCO2`, and `frac_doc`; and
- orphan CLM-Microbe text-file records: `AOM` and `H2maxCH4` are read but never
  assigned, and `AcemaxCH4` lies beyond the 88-record read boundary.

The existing `ch4offline` control is active at the coupling level and remains a
namelist/runtime control as described in section 12.2; its presence in the old
parameter reader does not make it a new standard parameter-file variable.

### 10.6 Migration and provenance procedure

Write a conversion/audit utility that reads the historical text file, an
archived run copy, and the old physiology NetCDF, then augments a copy of the
standard ELM parameter file and emits a report of:

- active variables migrated and their final NetCDF names;
- existing standard ELM variables reused;
- hard-coded coefficients promoted;
- ignored/inactive/orphan legacy values;
- missing values requiring investigator input;
- duplicate or conflicting definitions; and
- all unit conversions.

The utility must never modify the canonical input file in place. It writes a
new versioned ELM parameter file with global attributes recording the source
ELM parameter-file checksum, CLM-SPRUCE revision, conversion-tool revision, and
microbe/methane schema version.

The original CLM-SPRUCE reader/file pair cannot define a canonical mapping:
file records 84-86 are named `k_dom`, `k_bacteria`, and `k_fungi`, but the
reader assigns those positions to `dom_diffus`, `m_Fick_ad`, and
`m_dPlantTrans`; records after the fixed read boundary cannot supply the PFT
cascade values. The linked CLM-Microbe 89-record file corrects those first 86
positions and is therefore the Phase 3 default baseline. The conversion utility
must still report its three trailing orphan records and must never use it as a
source for absent Phase 2 or new C:P parameters. An archived successful run
directory or an investigator-approved table remains required before freezing
production values.

Initialization validates finite values, units/schema version, PFT coverage,
positive half-saturation and rate constants, yields and path fractions, pH
ordering, nonnegative stocks, and exact donor-path closure before the first
timestep.

## 11. Initialization, history, and restart

### 11.1 Lifecycle

Add one instance to `elm_instMod` whose lifecycle mirrors current ELM types:

- `InitAllocate` only when enabled;
- `InitCold` for cold/arb initialization;
- `InitHistory` only when enabled;
- `Restart(flag='define'|'write'|'read')`; and
- `ReadParams` conditionally from the standard ELM parameter file.

Step 5 makes the two executing backends exclusive without removing the legacy
implementation or its namelist capability. The existing `ch4_vars` allocation,
cold initialization, and restart lifecycle temporarily remain active whenever
`use_lch4` is true because standard decomposition and
nitrification/denitrification still consume its oxygen/anoxia fields. Legacy
methane history fields are registered only for the legacy backend. In revised
mode the adapter overwrites the shared oxygen subset from revised O2 state each
timestep; no legacy methane production, oxidation, or transport routine
executes. Moving the shared oxygen interface into a backend-neutral type is a
later structural refactor, not part of the Step 5 science change.

### 11.2 Restart compatibility

The revised restart schema includes:

- the three new C/N/P decomp pools through existing named-pool machinery;
- acetate, four functional-guild biomasses, and gas state by partition/layer;
- lagged redox/saturation state that affects future evolution; and
- a microbial-methane restart schema version.

Supported restart operations are:

- revised mode to the same revised mode: exact restart;
- cold start in revised mode: documented seeds/defaults; and
- legacy/disabled mode to itself: current behavior unchanged.

Switching an existing legacy restart into revised mode, or a revised restart
into legacy mode, changes the state vector and must fail with an actionable
message unless the user supplies an explicit offline conversion product. A
converter may seed the new pools and document the mass redistribution, but the
model must not do that silently.

### 11.3 History output

Default revised-mode output should be sufficient to close budgets without
creating an impractically large history file:

- DOM, bacteria, and fungi C/N/P;
- acetate and four functional-guild biomass C;
- CH4/O2/CO2/H2 profiles or column totals;
- acetoclastic and hydrogenotrophic production;
- aerobic and anaerobic oxidation;
- diffusive, ebullitive, and aerenchyma CH4 fluxes;
- total surface CH4 flux and CO2 correction; and
- C/N/P and gas-reaction residuals.

Full saturated/unsaturated layer diagnostics are enabled by
`hist_wrtmicrobediag`. Field names should be new and explicit; do not reuse a
legacy `CH4Mod` history name for a quantity with different semantics.

## 12. Isotopes and atmospheric coupling

### 12.1 C13 and C14

The old code executes its reaction model independently for bulk, C13, and C14.
That can select different limiting rates and violate isotope-to-bulk
consistency. The new design computes bulk rates once, then transfers isotopes
along the accepted bulk fluxes with process-specific fractionation factors and
zero-mass guards.

Implementation order:

1. bulk C with a runtime error if `use_c13` or `use_c14` is requested;
2. C13 transfers and emitted-CH4 diagnostics, validated against the latest
   CLM-SPRUCE reference; and
3. C14 only after confirming that the placeholder portions of the old code are
   scientifically required and defining restart/history expectations.

C13 is required before claiming full parity with the latest CLM-SPRUCE methane
work. C14 is not on the initial production critical path unless a target case
requires it.

### 12.2 Offline and online methane

Preserve the existing `ch4offline` meaning:

- offline: use prescribed atmospheric CH4, report diagnostic surface flux, and
  do not alter atmospheric CH4 or NEE; and
- online: consume atmosphere-provided CH4, export net CH4, and apply the
  carbon-conserving CO2 correction.

The revised backend should implement the current exchange interface, but the
first production qualification remains offline because current ELM documents
online mode as not yet functional. Enabling revised online coupling requires a
coupled carbon-budget test and atmosphere acceptance test; it is not implied by
merely populating `Fall_methane`.

## 13. Expected implementation surface

Likely existing files to change are grouped below. This is a dependency map,
not a commitment to edit every file.

| Area | Likely files/modules |
| --- | --- |
| Namelist/configuration | `namelist_definition.xml`, `namelist_defaults.xml`, `ELMBuildNamelist.pm`, `elm_varctl.F90`, `controlMod.F90` |
| Dimensions/topology | `elm_varpar.F90`, `CNDecompCascadeConType.F90`, `DecompCascadeCNMod.F90` or a new selected `MicrobeDecompMod` |
| Decomp calculations | `SoilLittDecompMod.F90`, C/N/P state and flux types/update modules, vertical transport, summaries, balance checks |
| Parameters/PFT data | `readParamsMod.F90`, potentially `pftvarcon.F90`, the standard ELM parameter NetCDF, plus the new parameter module and audit/conversion utility |
| Instance lifecycle | `elm_instMod.F90`, `elm_initializeMod.F90`, `restFileMod.F90` |
| Timestep dispatch | `elm_driver.F90`; avoid a second invocation from `EcosystemDynMod.F90` |
| Exchange/history | `lnd2atmType.F90`, `lnd2atmMod.F90`, history registration in the new type; preserve coupler field names |
| New science | the modules listed in section 7.1 |

All current generic pool loops and every allocation dimensioned by
`ndecomp_pools` or `ndecomp_cascade_transitions` require review even if no edit
is ultimately needed.

## 14. Implementation sequence and gates

Do not claim the combined revised microbial-methane backend until the
end-to-end path passes its gate. Phase 2 site tests run the microbial cascade
with the established methane backend fully initialized so oxygen limitation
and nitrification/denitrification remain active. This is an explicitly
temporary integration bridge, not the Phase 3 methane implementation.

### Phase 0: freeze references

- Archive disabled-mode E3SM history/restart baselines for methane off and
  current methane on.
- Recover the exact parameter files and configuration of a successful
  CLM-SPRUCE microbial run.
- Generate one-timestep and short-run golden vectors for pool states, reaction
  rates, and surface flux components.
- Record active CPP flags and whether the duplicated old call sites were both
  reached in that executable.

Gate: reference outputs and parameter provenance are reproducible.

### Phase 1: configuration and dormant type

- Add the option, validation, type skeleton, and conditional lifecycle.
- Keep the option false in all defaults/compsets.
- Add disabled-mode B4B tests before changing the pool topology.

Gate: all disabled-mode comparisons are exact and invalid combinations fail.

### Phase 2: microbial decomposition cascade

- Add conditional pool dimensions/indices and pool roles.
- Implement 11 pools and 47 transitions with named parameters.
- Route C, N, and P and audit fire, transport, erosion, subgrid, spinup, and
  restart assumptions.
- Preserve CLM-SPRUCE's 10x DOM diffusivity in the ordinary SOM vertical
  transport solver as a named parameter. DOM remains on this shared C/N/P path
  in Phase 3; acetate and gas use the selectable revised-methane transport
  mapping without duplicating movement of the authoritative DOM pool.
- Validate the cascade independently of revised methane reactions.

Gate: closed decomp tests conserve C/N/P, all pools stay nonnegative, and the
target CNP site case completes with balance checks enabled.

Implementation status: the Phase 2 branch contains the conditional topology,
rates, named parameter reader/injector, pool-role policies, and independent
closed-box tests. For site integration, the established methane backend is
fully initialized and continues to own its state, history, restart, gas
transport, and oxygen/anoxia coupling while the microbial cascade supplies the
generic decomposition respiration inputs. Phase 3 replaces this bridge with
the revised backend only when `use_microbe_methane=.true.`; the established
`CH4Mod` path remains a permanent supported capability when the option is
false. Site runs with the test parameter manifest are integration
evidence only; unverified PFT values and CNP-only assumptions still require
recovery or science-owner approval before production use.

### Phase 3: bulk revised methane

Implement Phase 3 as separately reviewable steps:

1. Parameter and state foundation: conditionally read the 64 named variables
   from the standard ELM parameter file; allocate saturated/unsaturated acetate,
   four functional-guild, and four bulk-soil gas-inventory states; add cold-start,
   inactive-by-default history fields, and complete restart I/O. Keep legacy
   `CH4Mod` dispatch unchanged. **Implemented on the Phase 3 Step 1 branch.**
2. Standalone reaction kernel: compute named tendencies without mutating state,
   with unit and closed-box tests. **Implemented on the Phase 3 feature
   branch.**
3. Conservative saturated/unsaturated repartition and gas transport.
   **Implemented on the Phase 3 feature branch.**
4. State limiting/commit, ELM C-budget integration, and atmosphere-facing CH4
   and NEE terms. **Implemented on the Phase 3 feature branch; the opt-in
   interfaces are activated by Step 5 dispatch.**
5. Explicit dispatch at the existing methane call site: select the revised
   backend when `use_microbe_methane=.true.` and the legacy backend otherwise.
   Both implementations remain in the model; neither may execute twice.
   **Implemented on the Phase 3 feature branch, including revised carbon
   accounting, atmosphere export, and nitrification O2 consumption.**

Gate: golden-vector process parity, closed-box conservation, exact restart, and
short SPRUCE integration tests pass.

### Phase 4: diagnostics, conversion, and site validation

- Complete history/restart schema and parameter audit/conversion tooling.
- Validate at US-SPR against CLM-SPRUCE and available observations.
- Document expected differences caused by current hydrology, soil layering,
  forcing, and numerical precision.

Gate: agreed site metrics and budgets pass without unexplained residuals.

### Phase 5: isotopes and extended configurations

- Implement C13 from accepted bulk fluxes.
- Decide C14 scope.
- Qualify online atmosphere coupling, topounit variants, crops, and additional
  nutrient competition methods independently.
- Add OpenACC/GPU optimization after CPU numerical behavior is frozen.

## 15. Verification plan

### 15.1 Disabled-mode regression

For both omitted and explicit-false new option:

- methane off, current CTC CNP;
- current methane on, CTC CNP (`I1850CNPRDCTCBC` family);
- current methane on, Century decomposition; and
- restart/continue for representative site and global configurations.

Acceptance: bit-for-bit identical history and restart output to the Phase 0
baseline. No new input file is required and no new restart/history field appears.

### 15.2 Unit and component tests

- configuration truth table and incompatible-mode errors;
- 11-pool index/name/role initialization and 47-transition mapping;
- donor mass equations and CUE/path-fraction bounds;
- one-layer legacy-to-ELM potential-rate parity after explicit unit conversion;
- one-layer closed-box reaction atom balance;
- DOM C/N/P consumption and mineralization balance;
- nonnegative limiting with two or more competing consumers;
- saturated/unsaturated repartition at `f_sat=0`, `f_sat=1`, and changing
  intermediate values;
- diffusion against an analytic profile and zero-gradient equilibrium;
- ebullition threshold and aerenchyma sign/unit tests;
- cold initialization and missing/invalid parameter errors; and
- isotope/bulk consistency once C13 is implemented.

The Phase 3 focused suite now includes the one-layer equation-level parity
vector. Ten directly comparable rates agree at floating-point tolerance. It
classifies the legacy AOM Celsius/Kelvin error as an intentional divergence and
checks its exact Q10 ratio. This does not replace the Phase 0 requirement for a
vector emitted by an archived CLM-SPRUCE executable.

### 15.3 Integration tests

Use `elm_olmt` to create fresh site cases by default, including one-day and
one-year smoke tests. Do not use `create_clone` as the routine site workflow:
OLMT setup generates case-specific domain, surface, and parameter inputs in the
run directory, and a clone can omit those files while retaining absolute paths
to its parent. A clone that reuses an executable is permitted only for a
deliberately limited smoke check after every generated input and absolute path
has been audited. Fresh OLMT cases may explicitly share a compatible executable.
The configured prefix/date and full case name identify the run directory;
generated inputs are written there directly. A shared `elm-olmt/temp` staging
file is prohibited because concurrent setup can corrupt the NetCDF output.

- exact restart (`ERS`) with revised mode;
- multi-instance/processor-layout reproducibility appropriate to ELM's current
  standard;
- one-day, one-year, and spinup segments at US-SPR;
- C/N/P balance checks with diagnostic residuals enabled;
- no NaN, infinity, or negative prognostic state;
- timestep sensitivity; and
- topounit water-table/inundation forcing without old `HUM_HOL` code.

### 15.4 Scientific parity and validation

Compare the new model to a frozen CLM-SPRUCE reference at the process level:

- DOM, bacteria, fungi, and SOM trajectories;
- acetate and functional-guild biomasses;
- acetoclastic versus hydrogenotrophic CH4 production;
- aerobic versus anaerobic CH4 oxidation;
- gas concentration profiles;
- diffusion, ebullition, aerenchyma, and total CH4 flux; and
- C13 signatures when enabled.

Bit-for-bit parity with the ten-year-old executable is not expected after
adapting hydrology and state infrastructure. Each material difference must be
classified as an intended infrastructure difference, an identified old-code
defect, or an unresolved port error before calibration begins.

Run the repository's ELM developer regression suite after component/site tests
on a supported machine, with pre-change baselines generated from the recorded
E3SM revision.

## 16. Major risks and mitigations

| Risk | Mitigation |
| --- | --- |
| Uncertain parameter values and units despite an aligned CLM-Microbe test baseline | Preserve per-value provenance, require archived-run evidence, and generate an audit report before freezing production inputs |
| Hidden C/N/P imbalance from direct DOM mutation | One authoritative DOM state, explicit tendencies, single commit, and residual history fields |
| New pool count changes unrelated modes | Conditional counts and allocations, false-path B4B tests before science work |
| Old C-only science does not define microbial phosphorus | Parameterize DOM/bacteria/fungi C:P, close P budgets, and require science signoff; keep methane guilds explicitly carbon-only |
| Double invocation of revised methane | One driver dispatch and no ecosystem-dynamics call |
| Hydrology divergence from CLM-SPRUCE | Use current ELM/topounit inputs and diagnose hydrologic versus reaction differences separately |
| Restart mass loss when switching modes | Reject cross-mode restarts unless an explicit offline conversion file is supplied |
| Lakes or external BGC receive incomplete state | Fail unsupported configurations and keep lake revised state/flux zero |
| Isotopes select inconsistent rates | Derive isotope transfers from one accepted bulk reaction calculation |
| Memory and GPU cost of partitioned layer state | Allocate only when enabled; establish CPU parity before accelerator optimization |

## 17. Decisions requiring science-owner confirmation

These decisions do not block the structural implementation, but they must be
closed before a production parameter set or full validation claim:

1. Production parameter values and effective units beyond the reproducible
   CLM-Microbe test baseline, including conflicts with declarations, the old
   CLM-SPRUCE file, PFT inputs, and archived runs.
2. C:P ratios and fixed-versus-floating stoichiometry for DOM, bacteria, and
   fungi in the target CNP model.
3. Whether methane functional guilds should remain carbon-only after parity or
   be extended to explicit N/P limitation in a later science change.
4. Required C13 and C14 scope for the first production release.
5. Whether online atmosphere coupling is a release requirement or a later
   qualification target.
6. Desired water-flux-driven lateral transport of DOM C/N/P and acetate, and
   whether the optional ELM multiphase gas mapping should also receive a
   lateral operator. Default Fickian gas inventories now use a generalized,
   conservative multi-topounit graph operator; the old advection/remapping code
   is still not suitable for direct porting.
7. The archived US-SPR case, forcing, initial/restart data, and diagnostic list
   that define the CLM-SPRUCE reference.

## 18. Acceptance criteria

The integration is complete only when:

- the new option defaults false and disabled ELM runs are B4B for both legacy
  methane states;
- enabled configuration is explicitly validated and unsupported combinations
  fail before allocation or stepping;
- DOM, bacteria, and fungi are standard, conservative C/N/P decomp pools;
- the revised backend runs exactly once per land timestep;
- default Fickian gas exchange conserves area-weighted inventory for two- and
  multi-topounit graphs with unequal weights and vertically offset soil grids;
- bulk revised methane closes all required budgets and restarts exactly;
- parameter provenance is resolved and the complete active inventory is stored
  as named, versioned variables in the standard ELM parameter NetCDF;
- current atmosphere field names and sign/unit contracts are honored;
- US-SPR process differences from CLM-SPRUCE are explained and accepted; and
- the applicable ELM regression suite passes against the frozen baseline.

## 19. Remaining scientific validation and parameter work

The implemented equations and interfaces are not yet a production-calibrated
model. The current Phase 2 and Phase 3 parameter files are integration-test
fixtures: they make every required input explicit and allow conservation,
restart, and long-run tests, but they do not establish that all values are the
ones used by a successful CLM-SPRUCE experiment or that they are appropriate
for ELM's current C-N-P formulation. Phase 3 schema v3 now makes the linked
CLM-Microbe runtime file the reproducible default for every active value it
contains and records executable-equivalent unit conversions; this removes
ambiguity in the integration-test fixture, not the need for scientific
validation. The remaining work must distinguish three classes of information:
recoverable legacy values, equation choices that require interpretation, and
genuinely new parameters introduced by the ELM integration.

### 19.1 Recover the historical scientific reference

The first priority is to locate an archived CLM-SPRUCE build and run directory
with its exact `microbepar_in`, physiology parameter file, initial conditions,
compiler flags, CPP options, and output. The 88-value positional reader paired
with the original 111-record input file is not sufficient evidence of the
values actually used. The aligned 89-record CLM-Microbe file is now the test
baseline, but an archived run is still needed to establish whether that file
was staged unchanged and which source revision produced accepted results. For
each active parameter, the audit must record the runtime value, effective units
after all source conversions, its equation and call site, and whether it came
from a text record, PFT parameter, source default, or literal. The audit must
also establish which historical `microbech4` call sites executed. Agreement
with a checked-in source file alone is not a scientific parity result.

The archived executable should be used to produce one-layer and one-timestep
vectors for DOM consumption, acetate and H2/CO2 production, both
methanogenesis pathways, both oxidation pathways, guild growth and mortality,
and diffusion, ebullition, and aerenchyma fluxes. Differences caused by the
intentional carbon-accounting and units repairs in sections 9.2 and 9.3 must be
quantified rather than tuned away. These vectors define which remaining
differences are porting errors and which are accepted corrections to the old
model.

### 19.2 Resolve and calibrate microbial decomposition parameters

The litter/SOM-to-DOM fractions, direct-stabilization fractions, bacterial and
fungal allocation traits, pool turnover rates, carbon-use-efficiency controls,
microbial-residue routing, DOM stabilization, and DOM vertical diffusivity all
need traceable values. Several currently used PFT turnover and routing values
are marked unverified because the legacy positional reader could not have
filled the named PFT arrays as implied by the source. They should first be
recovered from an archived run or investigator table. Only parameters that
cannot be recovered should be calibrated.

Calibration must use multiple observables rather than matching total soil C
alone. Useful constraints include litter mass loss, heterotrophic respiration,
dissolved organic carbon concentration or export, microbial biomass C,
bacterial-to-fungal allocation where data exist, and the vertical and
fractional distribution of SOM. SOM3 and SOM4 must be evaluated after final
spinup because their accelerated-spinup values are deliberately rescaled,
whereas DOM and microbial biomass are not accelerated. Parameter combinations
that reproduce total soil C while producing implausible DOM or microbial
stocks must be rejected.

The initial paired 50-year US-MOz test provides a concrete calibration target,
not a validated result. Relative to the legacy-CH4 case, the enabled case ended
with `+216%` litter C, `+180%` SOM C, `+29%` total column C, and `-10%`
vegetation C. Its DOM stock reached `759.94 gC m-2`, with bacterial and fungal
stocks of `24.40` and `40.27 gC m-2`. These large shifts are consistent with a
substantially altered decomposition cascade but are too large to accept
without observational constraints and an unaccelerated final-spinup test.
Routing fractions, DOM stabilization and transport, microbial turnover, and
nutrient stoichiometry should be evaluated jointly rather than calibrating
only total soil C.

### 19.3 Define the new phosphorus science

CLM-SPRUCE did not define phosphorus stoichiometry for the added DOM,
bacterial, or fungal pools. The following are therefore new ELM science
parameters, not values that can be recovered from the legacy methane module:

| Parameter | Required scientific decision and evidence |
| --- | --- |
| `cp_bacteria` | Bacterial biomass C:P, including whether it is fixed, PFT-dependent, soil-dependent, or allowed to acclimate; constrain with microbial biomass C and P measurements or an accepted synthesis |
| `cp_fungi` | Fungal biomass C:P and its variability; constrain separately from bacteria where fungal:bacterial composition data exist |
| `cp_dom` | Effective C:P of the reactive dissolved-organic pool; constrain with paired DOC and DOP measurements and recognize that bulk extractable DOM may not equal the modeled reactive fraction |

The cold-start P assigned to bacteria, fungi, and DOM follows these C:P ratios
and must be included in initialization budgets. We must decide whether fixed
ratios are adequate or whether flexible microbial and DOM stoichiometry is
needed. That decision affects immobilization, mineralization, plant-microbe P
competition, and the amount of DOM-C available to methane chemistry. It cannot
be made by simply selecting large C:P values to make P limitation disappear.

The use of `solutionp_vr` as the inorganic P counterpool also requires tests
over representative soil P regimes. Those tests must verify conservation when
DOM is produced, consumed, stabilized, and transported; competition with plant
uptake; interaction with sorption and secondary-mineral pools; and behavior
when solution P approaches zero. Sensitivity experiments should bracket
literature-supported bacterial, fungal, and DOM C:P ranges. The selected values
and uncertainty ranges need science-owner approval before the parameter file
is promoted beyond testing status.

### 19.4 Resolve revised-methane kinetics and redox assumptions

The growth, mortality, yield, half-saturation, Q10, pH, moisture, inhibition,
and biomass-floor parameters need a dimensional audit and, after legacy
recovery, calibration against process observations. Particular attention is
required for the interpretation of growth rates versus substrate-consumption
rates, the hydrogenotrophic yield after removal of the old extra factor of
four, the acetoclastic CH4/CO2 split, and the O2:CH4 ratio for aerobic
oxidation. Gross `MM_CH4_PROD` and `MM_CH4_OXID` should be evaluated separately;
matching net FCH4 can otherwise hide compensating errors in production and
oxidation.

Native ELM also needs a spatial soil-pH data path. The current adapter uses a
separate `microbe_methane_soil_ph_fallback` parameter because native
`chemstate_vars%soil_pH` is not populated. Phase 4 should identify an
appropriate soil-profile or surface-data product, define interpolation and
depth behavior, and test whether pH is prescribed or evolves. Methane kinetic
calibration must not use either the response optimum or the temporary fallback
to absorb a missing site-pH constraint.

Anaerobic methane oxidation remains scientifically incomplete because the old
model does not identify or budget its electron acceptor. The initial port can
retain the carbon-only AOM parameterization for parity, but a production claim
must explicitly state this limitation. Adding sulfate, nitrate, ferric iron,
or another acceptor would require new state, stoichiometry, parameterization,
and site data; it must not be inferred from the existing O2-inhibition scalar.

The first 50-year US-MOz AD-spinup integration test demonstrates that this was
an active calibration blocker rather than only a conceptual caveat. With the
former declaration-default fixture, additional methane-system C rose from
`0.38 gC m-2` after model year 25 to `316.51 gC m-2` after year 50, dominated
by anaerobic-methanotroph biomass. At year 50, gross CH4 oxidation was
`1.31e-6 gC m-2 s-1`, versus gross production of only
`4.13e-10 gC m-2 s-1`, and the atmospheric boundary supplied a net CH4 influx.
The run remained numerically stable and carbon-conservative, so this result
specifically flags the unconstrained AOM science and former parameters. The
linked-file baseline lowers AOM growth from 0.024 to 0.004 d-1 and yield from
0.40 to 0.15 while retaining 0.002 d-1 mortality; its unlimited net growth is
therefore 0.002 rather than 0.022 d-1. The corrected follow-up (`P3S`) gives
year-50 gross production of 21.842, oxidation of 1.270, and net emission of
20.557 g C m-2 yr-1, with 10.588 g C m-2 in additional methane-system state.
An earlier reported production value of 51.82 was a history-diagnostic error:
it counted the complete acetoclastic extent instead of its CH4-C yield. Before
scientific use, evaluate an explicit electron-acceptor limitation or a
documented AOM-off configuration and constrain AOM growth, death, yield, and
CH4 half-saturation against observations. Do not tune transport merely to hide
biomass growth.

The shared oxygen budget needs evaluation using observed or credible modeled
soil O2/redox profiles, nitrification and denitrification rates, and wetting or
water-table transitions. The former sequential coupling gave standard
respiration, roots, and nitrification first access to O2. A 50-year 2-by-2
sensitivity using
`K_CH4={1000,1}` and `K_O2={4000,4}` mmol m-3 produced identical output in all
four cases. A runtime probe confirmed that the low values were read, but the
standard-demand bridge had reduced pre-reaction top-layer O2 to exactly zero;
gas transport restored O2 only after the reaction call. The first common-
stress implementation then exposed a second coupling defect: it read the
legacy-only `col_cf%rr_vr` profile while that field still held its fill value,
creating a fictitious root O2 demand and stresses near `1e-40`. The corrected
adapter reconstructs the profile from patch respiration, root fraction, and
patch weight. All potential aerobic demands now share one O2 stress, microbial
aerobic rates use it immediately, and ELM decomposition uses it on the next
timestep.

A clean paired 50-year retest confirms that the Monod constants are active.
With `K_CH4/K_O2=1000/4000 mmol m-3`, aerobic biomass remained at its
functional floor and year-50 aerobic oxidation was `2.83e-11 g C m-2 yr-1`.
With `1/4 mmol m-3`, maximum saturated aerobic biomass reached
`0.01168 g C m-3`, aerobic oxidation reached `5.679 g C m-2 yr-1`, and net
`FCH4` fell from `5.453` to `0.014 g C m-2 yr-1`. Hydrology and maximum O2 were
effectively identical. Compare this lagged shared-stress treatment with bounded
subcycling if rapid water-table changes produce timestep sensitivity. Annual
post-transport O2 history alone remains insufficient to validate the reaction
environment.

There is also a factor-of-1,000 contradiction between the legacy documentation
and executable. Its parameter declaration labels the raw `1.0` and `4.0`
constants as mmol m-3, while its active arithmetic compares them directly with
concentrations that are numerically mol m-3 (or mmol L-1). Consequently,
1,000/4,000 mmol m-3 in ELM unambiguously reproduces the legacy executable,
whereas 1/4 mmol m-3 follows the stated parameter units. Retain the former as
an execution-parity baseline and the latter as a candidate science setting,
but treat the large paired response as an identifiability constraint rather
than evidence that either is calibrated.

### 19.5 Calibrate transport and hydrologic controls

DOM/acetate diffusivity, the aqueous-gas diffusion multiplier, plant transport
coefficient, CH4 and H2 thresholds, root e-folding depths, ebullition depth
scale, saturation threshold, thaw threshold, and plant O2/CO2 factors need to
be checked against their legacy units and current ELM hydrology. The reused ELM
`f_sat`, `satpow`, `scale_factor_gasdiff`, `scale_factor_liqdiff`, and
`organic_max` controls also require sensitivity tests in the revised backend;
their legacy-CH4 defaults are not automatically calibrated for the microbial
reaction model. Surface CH4 must be decomposed into diffusion, ebullition, and
aerenchyma components, with saturated and unsaturated contributions reported
separately. Water-table position, inundated fraction, air-filled porosity, ice
state, snow and surface-water resistance, rooting depth, and plant functional
type should be validated before modifying reaction kinetics to fix a flux
timing error that is actually hydrologic or transport-driven.

The new solute-transport parameters require independent research and
calibration: molecular diffusivities for DOM and acetate, their mobile
fractions, longitudinal dispersivity, the saturation/tortuosity exponent, and
the minimum liquid fraction used as a numerical mobility threshold. None of
these seven values exists as a validated parameter in CLM-Microbe. Use soil-
solution profiles, lysimeter or porewater data, tracer experiments, and runoff
or drainage DOC/DON/DOP fluxes to constrain them. Evaluate whether DOM mobility
should vary with substrate class, mineral versus peat soil, pH, ionic strength,
and adsorption; a single equilibrium mobile fraction is only the first
hypothesis. The existing new microbial P controls and DOM C:P assumptions must
be calibrated jointly with organic-P transport and solution-P competition,
because plausible carbon profiles alone do not validate the phosphorus budget.

Hydrologic coupling still needs explicit tests for bottom drainage, saturation-
excess exfiltration, ponded surface water, and runoff routing. Ground
evaporation must retain solute, as implemented, but true liquid exfiltration
should transfer mass to a surface aqueous pool rather than delete it. Likewise,
subsurface drainage that is diagnosed after the current reaction call must be
connected at the correct operator-split point before interpreting leaching
fluxes. Compare diffusion-only, advection-only, and combined cases under dry,
steady saturated, pulsed infiltration, and fluctuating-water-table conditions;
require grid and timestep convergence and exact C/N/P closure in each case.

The CLM-Microbe Fickian and optional ELM multiphase mappings must be compared as
separate transport hypotheses. Existing coupled runs remain useful for reaction
and O2-competition comparisons only when their selected mapping is recorded;
their net or pathway surface fluxes must not be mixed across modes. Repeat the
execution-parity 50-year US-MOz and SPRUCE cases in both modes with the
`MM_CH4_SURF_*` fields. A reaction-off column
test should first recover atmospheric equilibrium across a gas/aqueous
interface, conserve inventory under a moving water table, and converge across
timestep. Then require a plausible unsaturated diffusion sink in an upland
case without introducing a separate high-affinity oxidation formulation.

The linked CLM-Microbe tree provides an additional structural sensitivity
target: compare its continuous empirical inundated fraction, all-layer
unsaturated gas exchange, one-percent frozen exchange, and atmospheric-CH4
oxidation diagnostic against the ELM adapter's hydrology and conservative
top-boundary transport. These behaviors should be tested separately; they
should not be bundled into a kinetic calibration or copied without closing the
gas and carbon inventories.

US-SPR remains the key scientific comparison because it motivated the source
model, but US-MOz is the simpler integration and parameter-sensitivity site.
The recommended sequence is: use US-MOz to establish numerical stability,
budgets, and interpretable production/oxidation behavior; reproduce the
archived US-SPR configuration; then add independent wetland and upland sites
that span water-table, temperature, vegetation, and nutrient regimes. A
site-specific fit at SPRUCE alone is not evidence of transferable parameter
values.

### 19.6 Required validation experiments and release evidence

Before calibration, the implementation must pass exact restart, timestep
convergence, C/N/P and revised-carbon closure, nonnegative-state, and long-run
stability tests with the option both off and on. The disabled path must retain
the recorded B4B behavior. Enabled runs should include unaccelerated segments
after AD and final spinup so slow SOM, fast microbial/DOM state, vegetation,
and nutrients can adjust on compatible clocks.

Scientific evaluation should report at least DOM, bacteria, fungi, litter C,
SOM1-4, mineral N and solution P, GPP, NPP, HR, NEE, gross CH4 production,
gross CH4 oxidation, bulk-soil CH4 and O2 inventory profiles, phase-equivalent
mobile concentrations where needed for interpretation, and net CH4 exchange.
For the revised backend, production must be split into acetoclastic and
hydrogenotrophic components and oxidation into aerobic and anaerobic
components in detailed diagnostic runs, even if routine history stores only
the gross totals. Seasonal cycles and responses to water-table and temperature
changes are at least as important as annual means.

Calibration should use parameter priors and an identifiability analysis rather
than adjusting all microbial and methane coefficients simultaneously. A
sensible order is: hydrology and transport; microbial decomposition and DOM;
P stoichiometry and nutrient competition; gross methane production; methane
oxidation; then net surface exchange. Part of the available site record must be
withheld for validation, and parameter uncertainty should be propagated to
soil-C, NEE, and CH4 predictions. The production parameter file must record the
calibration data, objective functions, priors, posterior or selected values,
software revision, and validation results. Until these tasks are complete, the
reference parameter file remains explicitly labeled for testing only.

### 19.7 CLM-SPRUCE SOM initialization sensitivity

The source CLM-Microbe reaction loop raises saturated acetoclastic methanogen
biomass to `1.e-5 mol C m-3` before every reaction call when it falls below
that value. This is distinct from its cold initializer, which assigns every
guild `1.e-15 mol C m-3`, and it introduces carbon without a donor pool. A
100-year ELM sensitivity found only a 4.17% year-50 production increase from
the repeated reset, so that experimental capability was removed.

No high-biomass initializer or floor is retained. All methane guilds use the
CLM cold-start value `MFGbiomin` and subsequently evolve prognostically.

The relevant initialization sensitivity is instead the independently spun-up,
unaccelerated SOM state in CLM-SPRUCE's
`SPRUCE-finalspinup-peatland-carbon-initial.nc`. For a controlled comparison,
CLM hollow SOM1--SOM4 C and N map to the ELM boardwalk/fen and hollow
topounits, while CLM hummock maps to ELM hummock. ELM P is rescaled at fixed
C:P because the source CLM case does not prognose P. All non-SOM pools and
physical state remain copied from one common normal-mode ELM restart. This
isolates the effect of substrate initialization on CH4 production, oxidation,
and surface flux without confounding it with methane-guild seeding or AD
exit-spinup scaling.

The completed 100-year sensitivity rules out SOM initialization as the main
CLM/ELM production discrepancy. Copying CLM SOM increased grid-mean production
because it filled ELM's nearly empty boardwalk/fen topounit, but year-50
production in the shared hollow changed from 0.9917 to 0.8875 and the shared
hummock from 0.4454 to 0.4456 g C m-2 yr-1. The matched-SOM ELM values remain
far below the CLM-SPRUCE reference values of 14.18 and 13.18 g C m-2 yr-1.

The next parity target is saturated-area semantics. CLM-SPRUCE's `HUM_HOL`
path hard-codes `finundated` and `micfinundated` to 0.99, making both reference
columns approximately 99% saturated for methane state and rate aggregation.
ELM uses native `max(FSAT, frac_h2osfc)`; the matched-SOM year-50 values are
0.313 in the hollow and 0.0526 in the hummock. This is an intentional ELM
science change, but it must be isolated with a parity-only namelist option and
a one-layer reaction comparison before the remaining production difference
can be assigned to the reaction kernel, substrate coupling, or units.
