# Phase 3 revised-methane tests

This directory verifies the revised-methane parameter/state foundation, the
side-effect-free Step 2 reaction kernel, the conservative Step 3
repartition/transport kernels, the Step 4 ELM adapter and accounting
interfaces, and the Step 5 mutually exclusive backend dispatch. It supplements
CIME; it does not replace ELM build, run, or exact-restart tests.

## Docker architecture

`elmv3:latest` is the native ARM64 image. Use it for new ARM64 builds and
runs of ARM64 executables. `elmv3:amd64` is the x86-64 image. The shared
Phase 3 P3E executable under `phase3-artifacts/us-moz-20260918` is x86-64, so
P3I/P3J cases that reuse that executable must run with `elmv3:amd64` and
`--platform linux/amd64`. Running the P3E executable in `elmv3:latest` fails
with a Rosetta `/lib64/ld-linux-x86-64.so.2` loader error.

## Coupled site-case workflow

Create coupled site tests with `elm_olmt` by default, including short smoke
tests. Do not default to CIME `create_clone`: OLMT cases can reference generated
`domain.nc`, surface data, and modified parameter files in the original run
directory, while the clone does not reproduce all of those run-local inputs
and can silently retain absolute paths to the parent case. A clone with
`--keepexe` is acceptable only as a deliberately limited executable smoke test
after its input paths and generated files have been audited. Fresh `elm_olmt`
cases may still share a compatible executable explicitly.

For Docker cases, always run CIME case setup with `--disable-git`. OLMT
enables this by default for Docker; `disable_git = True` under
`[case_options]` makes the intent explicit. CIME's case-local Git bookkeeping
is unnecessary for bind-mounted development cases and can fail on
host/container ownership differences.

Generated `domain.nc`, `surfdata.nc`, parameter files, and optional dynamic-
landuse/FATES inputs are written directly into each case's run directory. The
run path already includes the configured prefix/date and full case name. Do
not use a shared `elm-olmt/temp` file for these inputs: concurrent case setup
can otherwise collide in the NetCDF/HDF5 writer.

When running independent one-rank OLMT cases concurrently in Docker, do not
assume that `taskset` on `case.submit` pins the actual model rank. OpenMPI may
discard the wrapper's intended placement and bind every spawned `e3sm.exe` to
core 0. Disable OpenMPI rebinding while retaining the wrapper affinity, for
example:

```console
OMPI_MCA_hwloc_base_binding_policy=none taskset -c 3 ./case.submit --no-batch
```

After launch, verify the affinity of each live `e3sm.exe`, not merely the
`case.submit`, shell, or `mpiexec` processes (`taskset -pc PID` on Linux). If
the ranks are already running on the same core, they can be moved safely with
`taskset -pc CORE PID`. Record the case-to-PID mapping from `/proc/PID/cwd`
before changing affinity.

The 85 scalar values in `phase3_reference_parameters.json` are traceable test
inputs, not a production calibration. The default reaction/transport baseline
is CLM-Microbe commit `9c2e0a048bb3799669d32b91e6cc76efb36d4b75`:
38 raw numbers come from its runtime `microbepar_in`, while values absent from
that file come from audited declarations or active literals in the same source
tree. The raw numbers are converted from the units in which the legacy
executable actually evaluates them to the explicit units of the ELM kernel:

- concentration controls: mol m-3 to mmol m-3, multiply by 1000;
- volumetric rates: mol m-3 s-1 to mmol m-3 d-1, multiply by 86,400,000;
- specific rates: s-1 to d-1, multiply by 86,400; and
- the functional-guild floor: mol C m-3 to g C m-3, multiply by ELM's
  `catomw = 12.011 g mol-1`.

The legacy source comments often label the unconverted numbers as daily or
millimolar, but the active code divides stored g C m-3 by 12 to obtain mol C
m-3 and multiplies tendencies by a timestep in seconds. Copying those numbers
unchanged into the unit-explicit ELM equations therefore does not reproduce the
legacy executable. Concentration-dependent source literals are converted by
the same rule. Active literals used by the port have been promoted to named
inputs. The two Kelvin-correction values still repair comparisons that used
Celsius literals with Kelvin soil temperature.

The CLM-Microbe positional reader consumes 88 of that file's 89 records and
assigns only the first 86, so `AOM` and `H2maxCH4` are read but unused and
`AcemaxCH4` is not read. The local CLM-SPRUCE parameter file contains inserted
decomposition records; it has been reordered so `dom_diffus`, `m_Fick_ad`, and
`m_dPlantTrans` again occupy the three positions consumed by the methane
reader. The decomposition values remain governed by the PFT parameter NetCDF.

The linked runtime values supersede the declaration defaults used by the
original Step 1 fixture. Their provenance begins with
`clm_microbe_9c2e0a_runtime_file` and records any executable-unit conversion.
This choice gives integration tests a reproducible upstream baseline; it is
not evidence that the values are scientifically calibrated for ELM-Peatlands.

Run focused tests:

```console
python3 -m unittest -v test_phase3_state.py
python3 -m unittest -v test_phase3_reactions.py
python3 -m unittest -v test_phase3_transport.py
python3 -m unittest -v test_phase3_state_update.py
python3 -m unittest -v test_phase3_adapter.py
```

`reaction_oracle.py` is an independent scalar implementation used for unit and
closed-box tests. The tests exercise every reaction and mortality path, require
represented-carbon closure, and confirm that proportional substrate limiting
cannot drive a state negative over the requested timestep. A one-layer parity
vector evaluates the raw CLM-SPRUCE molar/per-second equations and the
unit-converted ELM equations at the same physical state. Ten directly
comparable production, oxidation, and mortality rates must agree to floating
point tolerance. Legacy AOM is classified separately because its active source
subtracts 13.5 degrees C from a Kelvin soil temperature; the test requires the
resulting difference to equal that exact erroneous Q10 factor. This is an
equation/oracle parity test, not yet an archived CLM-SPRUCE executable vector.
CIME compilation remains the authoritative Fortran interface and dependency
check.

`transport_oracle.py` independently exercises saturated/unsaturated
area-transfer repartition, the default CLM-Microbe Fickian coefficients,
explicit closed-boundary acetate diffusion, capacity-aware implicit gas
diffusion, conservative two- and multi-topounit lateral gas diffusion,
an opt-in conservative DOM C/N/P adjacent-layer relaxation operator and the
grid-aware physical aqueous advection/diffusion operator,
bidirectional O2/CO2 aerenchyma
exchange, one-way threshold CH4/H2 plant emission, and threshold ebullition.
Tests cover the zero- and one-area limits, phase-equilibrated gas/aqueous
profiles, zero-gradient diffusion, nonnegative updates under stiff gas
transport, atmospheric uptake by surface diffusion, rejection of artificial
CH4 plant influx, surface-flux units/signs, and exact inventory closure.

Revised methane defaults to the CLM-Microbe aqueous Fickian mapping. Both area
partitions use the active source's `Fick_D_w*m_Fick_ad*(T/298)^1.87`
coefficient and direct concentration-gradient basis. The implementation keeps
the conservative backward-Euler solver and replaces the source's timestep
surface reset with a finite top-half-layer Fickian boundary to atmospheric
Henry equilibrium.

With `use_humhol=.true.`, this default path exchanges CH4, O2, CO2, and H2 on
ELM's arbitrary `regional_target_ti` graph. It uses actual topounit/column
weights and lateral distances, overlaps layers by absolute elevation, and
scales every pair flux by the smaller participating horizontal footprint before
limiting simultaneous outgoing transfers against each donor inventory. That
area factor prevents concentration spikes as one area partition approaches
zero. The lateral coefficient remains 0.1 times the vertical coefficient for
executed CLM-Microbe parity. `MM_LATERAL_C_FLUX` is positive into a column and
is passed to the column carbon-balance check; its gridcell-area sum closes to zero.
Water-flux-driven lateral DOM C/N/P and acetate transport remains deferred.

`use_microbe_nonbog_lateral_gas_transport=.true.` is the default and preserves
gas exchange on every configured topounit edge. Set it false to skip an edge
when either endpoint is non-bog while retaining bog-to-bog exchange. This is a
gas-only diagnostic switch: it does not disable lateral hydrology or alter
vertical gas and aqueous transport.

Set `use_elm_microbe_methane_transport=.true.` to select the alternative ELM
multiphase mapping. That path uses gas-phase diffusion in air-filled pores,
aqueous diffusion in saturated media, Henry capacity across phase interfaces,
and finite snow/pond/topsoil/boundary resistance. It reuses standard ELM
`f_sat`, `satpow`, `scale_factor_gasdiff`, `scale_factor_liqdiff`, and
`organic_max` parameters. Both modes retain the same reactions, conservative
state accounting, aerenchyma and ebullition kernels, and split history fields.
Gas-transport selection is independent of DOM transport. By default DOM C/N/P
remains on the Phase 2 standard ELM transport operator; the gas switch changes
acetate and gas transport only.

When `use_microbe_aqueous_transport=.true.`, the same conservative
backward-Euler water-flux operator now transports standard mineral NH4 and NO3
after the reaction update. NO3 is fully mobile. NH4 uses equilibrium linear
sorption: only `theta/(theta + rho_b Kd)` of the authoritative total NH4 pool
is represented in porewater and transported. The sorbed remainder stays in
the same layer; no duplicate dissolved-N state is created. Bottom NH4 and NO3
exports are added to `MM_AQUEOUS_N_EXPORT` so the existing column N balance
accounts for them. Molecular diffusivities and `Kd` are named standard ELM
parameter-file inputs in the Phase 3 manifest and are explicitly uncalibrated
science hypotheses.

For the CLM-source diagnostic, set
`use_clm_microbe_dom_relaxation=.true.`. This separate, default-off switch
applies the source's executed adjacent-layer relaxation timescale to DOM C/N/P
and the saturated and unsaturated acetate stores after the
reaction update to the shared DOM C, N, and P profiles. The port uses one
backward-Euler finite-volume matrix for all three elements, closed vertical
boundaries, and element-specific inventory checks. Thus it preserves the
source timescale for equal layers without copying the source's sequential,
non-conservative layer mutation. This test is especially relevant to peat
topounits, where ELM's standard peat branch sets the underlying SOM diffusion
coefficient to zero, so the nominal 10x DOM multiplier also multiplies zero.
The first paired smoke test is `20260920MM27` (enabled) versus `20260920MM28`
(disabled), both one-year SPRUCE cases initialized from the CLM-SOM restart
with forced 0.99 saturation. MM28 reproduces the earlier MM26 annual numeric
history fields exactly; MM27 completes with machine-precision carbon balance
and a strong downward redistribution of DOM C/N/P.
The 50-year continuation is `20260920MM29`. At year 50, relaxation raises
hollow/hummock gross production to 9.69/9.04 g C m-2 yr-1, about 68% of the
corrected pH-4.5 CLM-SPRUCE reference. Oxidation remains only 42--45% of CLM,
and the resulting surface flux is about 3.6 times CLM. This isolates gas-profile
transport and oxidation, plus unmatched initial DOM, as remaining parity
issues rather than supporting a larger DOM relaxation rate by itself.

For the physical scientific hypothesis, set
`use_microbe_aqueous_transport=.true.` and leave
`use_clm_microbe_dom_relaxation=.false.`. This default-off path removes DOM
from generic decomposition-pool transport and applies one conservative,
backward-Euler operator to DOM C/N/P and acetate. It uses actual ELM layer
thickness, liquid water, and downward-positive interlayer water flux; molecular
diffusion, saturation/tortuosity, temperature, thaw state, and longitudinal
dispersion determine the diffusive conductivity. Bacteria and fungi remain
soil-attached. Clean infiltration contains no solute, and negative top flux is
clipped so ground evaporation cannot remove DOM. Bottom downward flow exports
solute and is included in the C/N/P balance hooks.

The seven added solute parameters are initial test hypotheses, not values taken
from CLM-Microbe: DOM and acetate molecular diffusivity, their mobile fractions,
longitudinal dispersivity, saturation/tortuosity exponent, and the minimum
liquid fraction. The first coupled smoke case is OLMT case `20260920MM30`; it
demonstrated machine-precision operator and whole-column balance but exposed
the top evaporative-boundary issue, which was then corrected. Re-test the
corrected boundary and add explicit post-hydrology drainage/runoff coupling
before using aqueous export as a scientific result. Corrected OLMT case
`20260920MM31` completed with zero boundary export in its first year,
`CMASS_BALANCE_ERROR=-1.36e-14 g C m-2`, and DOM C/N/P transport residuals
near `1e-17`--`1e-19 g element m-2`.

Root access is supplied by the reconstructed peatland-root implementation,
`use_peatland_roots`. It resolves N demand by PFT, excludes nonvascular PFTs
from adaptive deep access, retains the prescribed root profile as the physical
constraint, and excludes saturated layers. It is enabled by default for
HUMHOL cases and is not a methane-specific control.

`state_update_oracle.py` independently exercises the Step 4 transaction
boundary. Tests require one shared DOM pool across the two area partitions,
pre-commit C/N/P and mineral-nutrient limiting, preservation of the functional
biomass floor without post-hoc clipping, and explicit conversions for stored
carbon, surface CH4, and the CO2 correction. The transaction routines return
candidate state. The Step 4 adapter gathers ELM decomposition, nutrient,
hydrology, temperature, root, atmosphere, snow, and pond state; validates the
complete column; and commits state and budget diagnostics once. Its interfaces
are called by `elm_driver.F90` only when `use_microbe_methane=.true.` and
`use_legacy_ch4_with_microbe=.false.`. Setting the latter option true restores
the Phase 2 attribution mode: the DOM-bacteria-fungi cascade stays active while
the retained legacy `CH4Mod` supplies methane. Revised CH4 and CO2 surface
fluxes feed the established atmosphere fields, and revised storage plus surface
carbon exchange feed the opt-in balance and monthly-budget arguments. Potential
nitrification, standard heterotrophic and root respiration, aerobic acetate
oxidation, and aerobic methane oxidation now share one finite-O2 limiter in
each layer partition. The common stress scales microbial aerobic rates in the
current reaction call and feeds ELM's established `o_scalar` interface on the
next timestep, following the lagged structure of legacy `CH4Mod`.
Pathway histories `MM_CH4_SURF_DIFF`, `MM_CH4_SURF_AERE`, and
`MM_CH4_SURF_EBUL`, each with `_UNSAT` and `_SAT` area-weighted contributions,
separate atmospheric diffusion from plant transport and bubbling.
The adapter constructs column root fractions directly from patch state because
the retained legacy solver is not dispatched, and all inactive revised state
is initialized to finite zero so standard single-precision history output is
safe.

The pre-shared-limiter US-MOz `K_CH4`/`K_O2` integration sensitivity exposed
an ordering issue, not a parameter-file problem. Four 50-year cases using `(1000,4000)`,
`(1,4000)`, `(1000,4)`, and `(1,4)` mmol m-3 were numerically identical. A
runtime check logged the intended low values but found zero O2 immediately
before methane reactions, after the standard-demand bridge and before gas
transport. Consequently aerobic oxidation was zero and the two aerobic Monod
parameters could not affect the run. Preserve this result as a coupling test;
do not treat it as evidence that either unit interpretation is scientifically
equivalent.

A subsequent shared-stress diagnostic found a second, concrete adapter defect:
`col_cf%rr_vr` is filled only by legacy `CH4Mod` and therefore retained its
initial fill value when revised mode bypassed that solver. The resulting
fictitious root O2 demand reduced stress to about `1e-40`. The adapter now
reconstructs vertically resolved root respiration directly from patch
respiration, root fraction, and patch weight. In clean paired 50-year runs,
the execution-parity `(1000,4000)` values retained aerobic biomass at its floor
and gave effectively zero aerobic oxidation, while the stated-unit `(1,4)`
values produced `5.679 g C m-2 yr-1` of aerobic oxidation and reduced `FCH4`
from `5.453` to `0.014 g C m-2 yr-1` at year 50. This confirms that the Monod
parameters are active after the coupling fix; it does not establish that the
low values are calibrated. See `microbe-methane-parameter-validation.md` for
the full stock, flux, water-table, and balance comparison.

Starting with the Phase 2 parameter file, add the revised methane variables in
the Docker environment:

```console
python3 phase3.py inject-parameters \
  --input /output/clm-parameters-phase2.nc \
  --output /output/clm-parameters-phase3-reference.nc
python3 phase3.py validate-parameter-file \
  --file /output/clm-parameters-phase3-reference.nc
```

The input is always copied; this tool refuses to edit the standard source file
in place. The output is still an ELM standard parameter NetCDF, with a schema,
source revision, legacy identifier, units, and provenance attached.
