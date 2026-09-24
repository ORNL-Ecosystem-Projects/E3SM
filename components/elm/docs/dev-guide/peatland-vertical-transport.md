# Peatland vertical transport

`use_peatland_vertical_transport` selects peat-specific vertical transport for
the litter and soil-organic-matter pools handled by `SoilLittVertTranspMod`.
It defaults to true when `use_humhol` is enabled and false otherwise. An
explicit false setting keeps standard ELM bioturbation and SOM advection in a
peatland hydrology case, which is useful for attribution tests.

The option is runtime-gated and requires `use_humhol`. With the option off,
the existing assignments of `som_adv_flux` and `som_diffus` are unchanged.
With the option on, each soil column obtains its topounit through
`col_pp%topounit`. A positive `TopounitPeatDepth` selects the peat coefficients;
a zero depth retains standard ELM coefficients. This avoids the historical
assumption that a particular column index is always the hummock or hollow and
works with any number or ordering of topounits.

## Coefficients

When `use_peatland_compaction_profile` is false, the compatibility advection
coefficient in topounit *t* is

```text
v_peat(t) = peat_som_adv_flux
            * max(TopounitPeatDepth(t), 0)
            / peat_adv_reference_depth
```

When `use_peatland_compaction_profile` is true, the density parameters instead
define a depth-dependent solid-carbon storage capacity:

```text
rho_target(z) = peat_compaction_surface_density
                + (peat_compaction_deep_density
                   - peat_compaction_surface_density)
                  * (1 - exp(-z / peat_compaction_efolding_depth))

F_excess(j) = max(rho_physical(j) - rho_target(j), 0)
              * dz(j) / peat_compaction_timescale_years

v(j+1/2) = F_excess(j) / rho_physical(j)
```

The modeled physical-equivalent density is the sum of all solid litter, CWD,
SOM, bacterial, and fungal C pools; DOM is excluded and remains under the
aqueous transport operator. Carbon above a layer's target capacity is buried
toward the next layer over the configured timescale. The surface and bottom
boundaries are closed, so this operation redistributes material without adding
or removing it. C, N, P, C-13, and C-14 move through the same conservative
matrix and retain pool stoichiometry and isotope composition.

During accelerated decomposition spinup, `rho_physical` is reconstructed with
the same pool acceleration factors (including the accumulated scalar after
year 40) used by vertical transport. Thus reduced raw AD slow-pool stocks do
not falsely suppress compaction. The variable-velocity transport operator also
includes the interface-flux divergence term that was unnecessary for ELM's
historical constant-velocity formulation.

The physical velocity remains common to all transported solid pools. During
accelerated decomposition spinup, the existing transport solver subsequently
multiplies each pool's velocity by its decomposition acceleration factor (and,
after the scalar accumulation period, divides by `scalaravg_col`). Thus the
different numerical AD velocities accelerate convergence; they do not assert
that SOM pools have different physical burial velocities. Normal and final
spinup use the common physical profile directly.

The initial capacity-overflow defaults are:

- `peat_som_diffus = 0 m2 s-1`;
- `use_peatland_compaction_profile = true` when peatland vertical transport is
  enabled;
- `peat_compaction_surface_density = 25 kg C m-3`;
- `peat_compaction_deep_density = 100 kg C m-3`;
- `peat_compaction_efolding_depth = 0.25 m`;
- `peat_compaction_timescale_years = 1 year`.

The four compaction scalars are read from the standard ELM parameter file;
they are not namelist options. This keeps the process switch in the namelist
while placing its scientific coefficients with the rest of ELM's calibrated
parameters.

`peat_som_adv_flux` and `peat_adv_reference_depth` now apply only to the
compatibility path. They do not affect capacity overflow. The fen/boardwalk
topounit receives peat transport because its mapped peat depth is positive;
transport is determined by peat presence rather than vegetation or the
`TopounitIsBog` label.

Set `use_peatland_compaction_profile = .false.` to recover the historical
constant-with-depth velocity and total-peat-depth scaling for controlled
comparisons.

Dynamic peat columns bury the full solid peat matrix, including CWD and living
microbial biomass. Standard columns retain ELM's historical exclusions. The
option does not add or remove material and requires no new balance-check term;
it only changes vertical redistribution within a column.

The following inactive-by-default history fields support validation:

- `PEAT_C_DENSITY`: physical-equivalent modeled solid-C density (g C m-3);
- `PEAT_C_TARGET_DENSITY`: prescribed capacity curve (g C m-3);
- `PEAT_C_BURIAL_FLUX`: downward excess-C flux leaving each layer
  (g C m-2 s-1);
- `SOM_ADV_COEF`: corresponding interface velocity (m s-1).

## Scientific status

The target densities, compaction e-folding depth, one-year overflow timescale,
and zero diffusivity are not yet independently validated general peatland
parameters. Validate them against peat bulk-density profiles, accumulation
ages, and vertical C/N/P stocks. The target density is a carbon-density target,
not dry bulk density; observational comparisons must account for peat carbon
fraction. Nonzero peat diffusion can be tested through `peat_som_diffus`, but
should be constrained rather than enabled solely to smooth profiles.

The current logic retains ELM's pre-existing cryoturbation selection: when a
column qualifies for the frozen-soil cryoturbation branch, that branch takes
precedence over peat coefficients. This ordering preserves historical model
behavior and should be reviewed separately for permafrost peatlands.
