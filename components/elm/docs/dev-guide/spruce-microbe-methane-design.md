# CLM-SPRUCE microbial decomposition and methane integration design

Status: proposed design; no implementation has begun
Date: 2026-09-17

## 1. Executive summary

This change should be implemented as a selective port and refactor, not as a Git
merge or rebase. The CLM-SPRUCE source predates the current ELM architecture and
uses global `clmtype` storage, CPP selection, positional text parameters, and
direct state mutation. Current ELM already contains a substantially evolved
`CH4Mod`, modern typed state, restart and history infrastructure, C-N-P
decomposition, topounit hydrology, and atmosphere coupling. Replacing those
facilities with the older CLM-SPRUCE implementation would regress current ELM
capabilities.

The proposed design adds one opt-in science switch:

```fortran
use_microbe_methane = .false.
```

`use_lch4` remains the existing umbrella switch for methane. When
`use_microbe_methane` is omitted or false, ELM follows its present call graph,
pool topology, parameter reads, restart schema, history registration, and
floating-point operations. When both switches are true, the new switch selects:

1. the CLM-SPRUCE-style microbial CTC cascade, adding DOM, bacteria, and fungi
   to the standard decomposition pools; and
2. a revised methane backend with acetate, explicit methanogenic and
   methanotrophic guilds, dissolved CH4/O2/CO2/H2, and gas transport.

The existing `CH4Mod` remains intact and remains the default methane backend.
The two methane backends must never run on the same soil column and timestep.

The primary production target is the E3SM-Peatlands C-N-P, relative-demand,
vertical CTC configuration. C:P behavior for the three new decomposition pools
is therefore part of the initial integration, not a deferred enhancement.

## 2. Source baselines inspected

| Tree | Revision inspected | Relevant source |
| --- | --- | --- |
| E3SM-Peatlands | branch `simplify/bog-zwt-no-perched`, commit `7b0703e07c7c933d53814f0c31856d3aed480d7b` | `components/elm` |
| CLM-SPRUCE | branch `master`, commit `a85800bad2c2a57abff77af36ccec319b202166a` | `models/lnd/clm/src/clm4_5` |

The important CLM-SPRUCE source files are `microbeMod.F90`,
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
- Preserve exact current behavior when the new option is disabled.
- Make all scientific parameters named, unit-documented, validated NetCDF
  inputs.
- Establish process-level CLM-SPRUCE parity before calibrating the port or
  changing its equations.

### Non-goals for the initial integration

- Replaying or merging the old CLM-SPRUCE commit history.
- Replacing or cleaning up the current legacy `CH4Mod` path.
- Porting the empty `microbeCN` or `microben2o` stubs.
- Porting `HUM_HOL` hydrology or its ad hoc lateral BGC code. The new backend
  consumes current ELM/topounit hydrologic state.
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

### 4.3 The old parameter input is not safe to reproduce

`microbevarcon.F90` declares `nummicrobepar = 88`, reads exactly 88 values from
`./microbepar_in`, ignores the names, and assigns by position. The inspected
`microbepar_in` contains 111 records. Its records 84-86 are `k_dom`,
`k_bacteria`, and `k_fungi`, while the reader assigns those positions to
`dom_diffus`, `m_Fick_ad`, and `m_dPlantTrans`. Later records are never read by
that routine. Many PFT-level cascade parameters are separately declared in
`pftvarcon.F90`; most were not present by name in the inspected SPRUCE physiology
NetCDF file.

Consequently, the checked-in text file, declared defaults, PFT reader, and the
values actually used by historical runs must be audited against archived run
directories before declaring a canonical parameter set. The port must not
encode the positional behavior accidentally.

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

## 5. Configuration contract

### 5.1 Namelist variables

Add to `elm_inparm`:

```fortran
logical :: use_microbe_methane = .false.
```

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

| `use_lch4` | `use_microbe_methane` | Result |
| --- | --- | --- |
| false | false | Current no-methane behavior, unchanged |
| true | false | Current `CH4Mod` behavior, unchanged |
| false | true | Configuration error with an actionable message |
| true | true | Microbial decomposition plus revised methane backend |

A string-valued `methane_model='legacy|microbial'` was considered but rejected
for this port. Adding a default-false boolean makes the backward-compatibility
boundary explicit and does not reinterpret existing user namelists or compsets.

### 5.3 Required configuration when enabled

The build-namelist layer and runtime initialization must both validate:

- `use_cn=.true.`;
- `use_vertsoilc=.true.`;
- `use_century_decomp=.false.` (the CTC cascade);
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

For the lowest-risk first integration, the revised backend may publish the net
column CH4 flux through the already allocated
`ch4_vars%ch4_surf_flux_tot_col` exchange field and publish the CO2 correction
through `lnd2atm_vars%nem_grc`. This leaves the existing `lnd2atmMod` aggregation
and coupler field names unchanged. All revised-backend internal state remains in
its own type.

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

### 8.3 C-N-P behavior

- DOM, bacteria, and fungi use the standard vertically resolved C/N/P decomp
  arrays as their only bulk state.
- Transfer C, N, and P through standard cascade flux arrays and state-update
  infrastructure wherever possible.
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
- dissolved CH4, O2, CO2, and H2.

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

### 9.4 Units and conservation

Define units at module boundaries and convert only there:

- ELM decomp pools: `g C|N|P m-3 soil`;
- functional biomass and acetate: one documented carbon unit, preferably
  `g C m-3 soil`;
- dissolved gases: `mol gas m-3` of the explicitly documented phase;
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

The first port should reproduce the CLM-SPRUCE explicit ordering but make the
ordering visible as staged tendencies. If a reaction can consume more substrate
than available over an ELM timestep, scale all competing consumers
proportionally. Do not add an implicit solver or subcycling until parity vectors
exist; either change would combine a numerical-method change with the port.

After parity, timestep-convergence tests may justify bounded subcycling as a
separate change.

## 10. Parameters and input data

All parameters used by the enabled feature will be migrated into ELM's standard
parameter NetCDF. There will not be a feature-specific parameter file. The
existing `paramfile` selection, staging, provenance, PFT dimension, and parallel
read machinery will be used. The new variables are conditionally required and
read only when `use_microbe_methane=.true.`, preserving compatibility with
existing parameter files and disabled-mode B4B behavior.

The tables below are the source-level inventory. They use the CLM-SPRUCE
identifiers so that every equation can be traced during the port. An
implementation may adopt clearer ELM-style NetCDF names, but each renamed
variable must carry a `legacy_name` attribute. Each new variable also needs
`units`, `long_name`, valid-range, source-revision, and provenance metadata.
The authoritative values and units remain subject to the audit in section 10.6;
source declaration defaults and values in the checked-in positional text file
must not be treated as interchangeable.

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
| `m_rf_s1m`, `m_rf_s2m`, `m_rf_s3m`, `m_rf_s4m` | PFT | SOM1-SOM4 fractions routed to microbial uptake |
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
| `m_dPlantTrans` | Root/aerenchyma transport coefficient; migrate |
| `g_dMaxH2inWater` | Dissolved-H2 threshold for plant transport; migrate |
| `atmch4`, `atmo2`, `atmco2`, `atmh2` | Fallback atmospheric mixing ratios; migrate, but current ELM atmosphere forcing takes precedence when supplied |
| `Fick_D_w(1:4)` | Shared water diffusion coefficients for CH4, O2, CO2, and H2; reuse the current ELM equivalent rather than duplicate |
| `Henry_C_w(1:4)`, `Henry_kHpc_w(1:4)`, `kh_tbase` | Shared Henry-law coefficients and reference temperature; reuse current ELM equivalents |
| `rgasLatm` | Gas constant; use the current ELM physical constant, not a parameter-file duplicate |

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

Cold-start acetate and all four dissolved-gas concentrations are zero in the
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
`HUM_HOL` layer geometry, lateral-exchange literals, and hard-coded SPRUCE grid
indices are deliberately excluded because that hydrology is not being ported.

### 10.4 Parameters already owned by ELM

The enabled modules also consume ordinary ELM parameter/state inputs such as
soil pH, porosity and suction parameters, litter and SOM decay rates, CWD
cellulose/lignin fractions, `ksomfac`, root fraction, PFT identity, and layer
geometry. Those remain owned by their current ELM structures and standard
parameter variables. They must be recorded as interface dependencies in the
implementation, but must not be copied into a second microbial namespace.

Similarly, `ch4offline`, `allowlakeprod`, history controls, and
`use_microbe_methane` are runtime controls, not scientific parameter-file
variables.

### 10.5 Declared or supplied legacy values not active in the inspected path

For completeness, the following names occur in `microbevarcon.F90` or the
111-record `microbepar_in`, but are not referenced by active calculations in
the inspected revised-methane/microbial-cascade path. They will be reported by
the conversion audit and will not be added to the standard ELM parameter file
unless an archived reference executable demonstrates an active use:

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
- orphan text-file records not declared as active parameters: `AOM`,
  `H2maxCH4`, and `AcemaxCH4`.

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

Because the inspected reader consumes names only as comments, declares 88
records, and is paired with a 111-record file, it cannot define the canonical
mapping. In particular, file records 84-86 are named `k_dom`, `k_bacteria`, and
`k_fungi`, but the reader assigns those positions to `dom_diffus`, `m_Fick_ad`,
and `m_dPlantTrans`; records after the fixed read boundary cannot supply the PFT
cascade values. The conversion utility must flag this mismatch rather than
reproduce it. An archived successful SPRUCE run directory or an
investigator-approved table is required before freezing production values.

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

Legacy `ch4_vars%InitHistory`, `InitCold`, and `Restart` remain selected for
`use_lch4 .and. .not. use_microbe_methane`. Its allocation behavior need not be
refactored during this work.

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

Do not expose `use_microbe_methane=.true.` as usable until the end-to-end path
passes its gate. Intermediate branches may contain the option, but must fail
early with a development-only message rather than run a half-connected model.

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
- Validate the cascade independently of revised methane reactions.

Gate: closed decomp tests conserve C/N/P, all pools stay nonnegative, and the
target CNP site case completes with balance checks enabled.

### Phase 3: bulk revised methane

- Add acetate, functional guilds, reaction tendencies, partitioning, and gas
  transport.
- Call the backend once at the existing methane dispatch location.
- Publish current ELM methane/NEE exchange and balance terms.

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
- one-layer closed-box reaction atom balance;
- DOM C/N/P consumption and mineralization balance;
- nonnegative limiting with two or more competing consumers;
- saturated/unsaturated repartition at `f_sat=0`, `f_sat=1`, and changing
  intermediate values;
- diffusion against an analytic profile and zero-gradient equilibrium;
- ebullition threshold and aerenchyma sign/unit tests;
- cold initialization and missing/invalid parameter errors; and
- isotope/bulk consistency once C13 is implemented.

### 15.3 Integration tests

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
| Uncertain parameter values caused by the 88/111-record mismatch | Require archived-run provenance and generate an audit report before freezing inputs |
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

1. Canonical parameter values and units, given the mismatch among the old text
   file, reader, source defaults, PFT reader, and archived runs.
2. C:P ratios and fixed-versus-floating stoichiometry for DOM, bacteria, and
   fungi in the target CNP model.
3. Whether methane functional guilds should remain carbon-only after parity or
   be extended to explicit N/P limitation in a later science change.
4. Required C13 and C14 scope for the first production release.
5. Whether online atmosphere coupling is a release requirement or a later
   qualification target.
6. Desired lateral transport of DOM and dissolved gases in topounit/hillslope
   cases; the old `HUM_HOL` implementation is not suitable for direct porting.
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
- bulk revised methane closes all required budgets and restarts exactly;
- parameter provenance is resolved and the complete active inventory is stored
  as named, versioned variables in the standard ELM parameter NetCDF;
- current atmosphere field names and sign/unit contracts are honored;
- US-SPR process differences from CLM-SPRUCE are explained and accepted; and
- the applicable ELM regression suite passes against the frozen baseline.
