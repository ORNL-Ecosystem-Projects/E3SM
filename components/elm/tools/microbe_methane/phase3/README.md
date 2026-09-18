# Phase 3 revised-methane tests

This directory verifies the revised-methane parameter/state foundation, the
side-effect-free Step 2 reaction kernel, the conservative Step 3
repartition/transport kernels, and the Step 4 ELM adapter and accounting
interfaces. It supplements CIME; it does not replace ELM build, run, or
exact-restart tests.

The 64 scalar values in `phase3_reference_parameters.json` are traceable test
inputs, not a production calibration. Values tagged `clm_spruce_declared_default`
come from the declarations in CLM-SPRUCE `microbevarcon.F90`; because its
positional file reader is inconsistent with the checked-in input, those values
are not assumed to reproduce an archived SPRUCE run. Active literals have been
promoted to named inputs. The two Kelvin values tagged
`corrected_celsius_kelvin_source_bug` correct comparisons that used Celsius
literals with Kelvin soil temperature.

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
cannot drive a state negative over the requested timestep. The Fortran kernel
is not dispatched by ELM in Step 2; CIME compilation remains the authoritative
Fortran interface and dependency check.

`transport_oracle.py` independently exercises saturated/unsaturated
area-transfer repartition, finite-volume vertical diffusion, signed
aerenchyma exchange, and threshold ebullition. Tests cover the zero- and
one-area limits, analytic and zero-gradient diffusion, nonnegative donor
limiting, surface-flux units/signs, and exact inventory closure.

`state_update_oracle.py` independently exercises the Step 4 transaction
boundary. Tests require one shared DOM pool across the two area partitions,
pre-commit C/N/P and mineral-nutrient limiting, preservation of the functional
biomass floor without post-hoc clipping, and explicit conversions for stored
carbon, surface CH4, and the CO2 correction. The transaction routines return
candidate state. The Step 4 adapter gathers ELM decomposition, nutrient,
hydrology, temperature, root, atmosphere, snow, and pond state; validates the
complete column; and commits state and budget diagnostics once. Its interfaces
are not called by `elm_driver.F90` yet: Step 5 will select revised methane or
the retained legacy `CH4Mod` backend from the namelist option.

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
