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

The peat advection coefficient in topounit *t* is

```text
v_peat(t) = peat_som_adv_flux
            * max(TopounitPeatDepth(t), 0)
            / peat_adv_reference_depth
```

The compatibility defaults reproduce the final pre-microbial Peatlands branch:

- `peat_som_adv_flux = 0.0004 m yr-1` (stored in the namelist in m s-1);
- `peat_adv_reference_depth = 3 m`;
- `peat_som_diffus = 0 m2 s-1`.

For the current SPRUCE surface depths of 2.95, 3.00, and 3.15 m, these defaults
give burial velocities of approximately 0.393, 0.400, and 0.420 mm yr-1. The
fen/boardwalk topounit receives peat transport because its mapped peat depth is
positive; transport is determined by peat presence rather than vegetation or
the `TopounitIsBog` label.

All non-coarse-woody-debris decomposition pools use the same coefficient
matrix. Consequently the existing generic transport machinery applies the
peat coefficients consistently to C, N, P, C-13, and C-14 pools. The option
does not add or remove material and requires no new balance-check term; it only
changes vertical redistribution within a column.

## Scientific status

The 0.4 mm yr-1 burial rate, linear peat-depth scaling, 3 m reference depth,
and zero diffusivity are compatibility assumptions inherited from the SPRUCE
development branch, not independently validated general peatland parameters.
In particular, present-day peat thickness is not necessarily a mechanistic
predictor of contemporary accumulation rate. Before using this option outside
its calibration context, validate the burial rate and profile evolution
against peat accumulation ages, bulk-density profiles, and vertical C/N/P
stocks. Nonzero peat diffusion can be tested through `peat_som_diffus`, but
should be constrained rather than enabled solely to smooth profiles.

The current logic retains ELM's pre-existing cryoturbation selection: when a
column qualifies for the frozen-soil cryoturbation branch, that branch takes
precedence over peat coefficients. This ordering preserves historical model
behavior and should be reviewed separately for permafrost peatlands.
