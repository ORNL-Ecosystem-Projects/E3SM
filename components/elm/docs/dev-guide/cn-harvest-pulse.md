# CN harvest pulse

## Purpose

ELM normally interprets the annual harvest values in the land-use time series
as fractional mortality rates and applies them continuously through the model
year. For an annual input fraction `a`, the default instantaneous rate is

```text
m = a / seconds_per_year.
```

That behavior is appropriate for spatially aggregated land-use data, where an
annual harvested area can represent many events distributed through the year.
It is less appropriate for a site history that represents one known stand
disturbance. In particular, continuously applying `a = 0.99` leaves roughly
`exp(-0.99)` of the starting pool in the absence of other fluxes; it does not
remove exactly 99 percent at one event.

The optional CN harvest pulse treats the annual input as an event fraction. On
the configured date, ELM sets

```text
m = a / dtime
```

for one timestep. The ordinary state update then removes exactly fraction `a`
of each affected pool over that timestep. A pulse fraction outside `[0, 1]`
is rejected to prevent negative vegetation pools.

This capability reconstructs the historical SPRUCE stand-establishment
treatment without tying it to the hummock/hollow hydrology switch or a C
preprocessor definition. It is implemented in the generic dynamic-subgrid CN
harvest path and can therefore be used for other site histories when an event
interpretation is scientifically justified.

## Namelist controls

The controls are in the `dynamic_subgrid` group:

```text
use_cn_harvest_pulse = .false.
cn_harvest_pulse_month = 1
cn_harvest_pulse_day = 31
cn_harvest_pulse_tod = 0
```

`cn_harvest_pulse_tod` is seconds after midnight. The defaults reproduce the
legacy SPRUCE timing: January 31 at 00:00. The trigger detects the timestep
that crosses the requested time, so the event is not missed when the requested
time lies between timestep endpoints.

The option defaults to false and only affects non-FATES CN runs with
`do_harvest = .true.`. It may remain enabled in shared spinup executables or
phases where `do_harvest = .false.`; it has no effect there. FATES harvest
timing is managed by FATES and is intentionally unchanged.

## Carbon, nitrogen, and phosphorus accounting

The option changes the time distribution and event interpretation of the
existing annual harvest fraction, not the downstream routing. Existing ELM
harvest fluxes still move:

- leaves, fine roots, live stems, coarse roots, storage pools, and transfer
  pools to their established litter destinations;
- dead stem carbon and associated nutrients to the 10- and 100-year wood
  product pools according to the PFT product fraction; and
- excess maintenance-respiration carbon to the atmosphere.

The same mortality fraction is applied consistently to the corresponding C,
N, and P vegetation pools. Existing state updates, product pools, litter
routing, balance checks, `NEE`, and `NBP` accounting remain in use. Useful
diagnostics include `WOOD_HARVESTC`, `WOOD_HARVESTN`, `WOOD_HARVESTP`, the
vegetation carbon pools, `NEE`, and `NBP`. Subdaily output is required to see
the instantaneous pulse directly; monthly output preserves its time-integrated
flux but smooths the event in the reported mean.

## SPRUCE configuration and limitations

OLMT's `SPRUCE.cfg` enables the pulse explicitly and specifies January 31. The
annual harvest value still comes from `flanduse_timeseries`; the namelist does
not create or alter a harvest event. The historical US-SPR point-data recipe
places a 0.99 primary-forest harvest fraction at the 1974 land-use endpoint,
so this setting is a site-history choice rather than a general peatland
process. ELM's dynamic-harvest reader uses file year `Y+1` for model year `Y`;
therefore a value stored at the 1974 endpoint is applied during model year
1973. This follows the existing dynamic-land-use convention and is independent
of pulse timing within the model year.

The date is recurring: any modeled year with a nonzero annual harvest input
will pulse on the configured month and day. Month and day must exist in the
run calendar. The current interface provides one date for all grid cells and
harvest categories. Multiple event dates, gridcell-specific dates, and FATES
events would require a richer event input format and are outside this option.
