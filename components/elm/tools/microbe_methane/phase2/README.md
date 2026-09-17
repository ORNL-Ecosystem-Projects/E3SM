# Phase 2 microbial decomposition tests

This directory provides focused tests and parameter-file tooling for the
11-pool, 47-transition DOM/bacteria/fungi cascade. It supplements rather than
replaces CIME:

- `test_phase2.py` checks the parameter inventory, exact transition topology,
  donor-path closure, feedback connections, nonnegative closed-box evolution,
  and C/N/P conservation.
- `phase2.py` copies a standard ELM parameter NetCDF and injects the named
  Phase 2 variables. It never edits the input file in place.
- the Phase 0 CIME runner remains the source of disabled-mode build, run, and
  bit-for-bit regression coverage.

Run the focused tests from this directory:

```console
python3 -m unittest -v test_phase2.py
```

Create a development parameter file inside the ELM Docker environment:

```console
python3 phase2.py inject-parameters \
  --input /inputdata/path/to/current-parameters.nc \
  --output /output/elm-parameters-phase2-test.nc
python3 phase2.py validate-parameter-file \
  --file /output/elm-parameters-phase2-test.nc
```

The committed parameter values are test fixtures, not an approved production
calibration. The manifest assigns every parameter to one provenance class:

- active CLM-SPRUCE source literal;
- value named in `microbepar_in` but not proven to have reached the historical
  PFT arrays because of the legacy input mismatch;
- a test assumption that reuses the bacterial/fungal pool C:N literals; or
- a new CNP-only assumption with no CLM-SPRUCE analogue.

The injector writes that class as each NetCDF variable's
`source_provenance` attribute. Microbial C:P ratios require science-owner
review, and the unverified text values require an archived successful-run
parameter file before they can be called production values.

## DOM transport boundary

CLM-SPRUCE gave DOM the standard SOM advection path and ten times the SOM
diffusivity. Phase 2 preserves that distinction through the named
`dom_som_diffusion_multiplier`. CLM-SPRUCE also copied DOM into separate
saturated and unsaturated DOC/DON state and applied `dom_diffus` there. That
duplicated-state path is deliberately deferred to Phase 3, where DOM, acetate,
and dissolved gases can use one conservative, water-coupled transport design.

DOM-P is dissolved organic phosphorus in the decomposition-pool arrays; it is
not ELM's inorganic `solutionp_vr` pool. Cascade mineralization and
immobilization exchange P between those two states. Phase 2 transports DOM-P
with DOM and leaves solution-phosphate transport, sorption, and leaching on the
existing mineral-P path.

## Site integration mode

Phase 2 site tests select both `use_microbe_methane=.true.` and
`use_lch4=.true.`. The established methane backend remains fully initialized
during this phase so its oxygen state and nitrification/denitrification
coupling are preserved; only the decomposition cascade is replaced. This is a
temporary bridge, not the revised methane implementation planned for Phase 3.
Run site tests with `do_budgets=.true.`; accelerated-decomposition CN spinup
checks C and N, while a subsequent ordinary CNP segment is required to exercise
P conservation.
