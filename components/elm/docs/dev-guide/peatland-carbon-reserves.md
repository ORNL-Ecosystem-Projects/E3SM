# Peatland carbon reserves

## Biomass-scaled nonstructural-carbon target

On topounits with positive `peat_depth`, ELM treats `cpool` as a
biomass-scaled nonstructural-carbon (NSC) reserve. Its target is

```text
cpool_target = cpool_target_wood_frac
             * (livestemc + deadstemc + livecrootc + deadcrootc)
             + cpool_target_leafroot_frac * (leafc + frootc).
```

The initial parameters are `cpool_target_wood_frac = 0.03` and
`cpool_target_leafroot_frac = 0.10`, both in gC reserve per gC tissue. The
wood term uses all attached stem and coarse-root carbon because the empirical
whole-organ concentrations integrate radial declines from sapwood to
heartwood. This is an *effective* capacity: the model's single `cpool` remains
fully accessible, whereas some measured inner-wood NSC may not be.

The 3% wood fraction is anchored to an ecosystem-scaled mixed temperate forest
estimate of 3.6% NSC dry mass relative to biomass, which converts to roughly
3% on a C:C basis. An organ-resolved whole-tree budget for evergreen white
pine provides a lower reference: 2.6 kg of sugars and starch in 103 kg of
woody dry mass converts to about 2.1% NSC-C per unit woody C when assuming 40%
carbon in NSC and 48% carbon in dry woody tissue. Boreal conifer observations
also show large seasonal redistribution and
especially dynamic root reserves, so these initial fractions should be
treated as calibration priors rather than fixed truths. The 10%
leaf/fine-root fraction represents the higher soluble-sugar and starch
concentrations expected in metabolically active tissues. It is more weakly
constrained because the available whole-tree study found foliage to be a
small mass contribution and did not resolve fine roots.

For peatland vegetation, excess respiration is

```text
cpool_relative = min(cpool / cpool_target, cpool_xr_scale_max)
XR = cpool * br_xr * cpool_relative * temperature_scalar.
```

Thus `br_xr` retains its units of s-1 and is the fractional turnover
coefficient when `cpool` equals its target. At twice the target the fractional
coefficient doubles; because the substrate pool is also twice as large, total
XR is four times its value at the target. This quadratic total-flux response is
intentional: merely scaling an absolute XR flux linearly with `cpool` would
algebraically reproduce the old `XR = br_xr * cpool` equation and the target
would have no effect. Below target, turnover declines smoothly rather than
switching off at a hard threshold. The initial `cpool_xr_scale_max` is 10, so
the quadratic response transitions to a linear response above ten times the
target. This prevents a collapsing target from producing an unbounded rate
when nonwoody vegetation senesces and its leaf and fine-root carbon disappear.
If the target reaches zero while `cpool` remains positive, ELM uses the same
10-times maximum; if both are zero, the multiplier is immaterial and is set
to one.

Leaf and fine-root storage and transfer pools are deliberately excluded from
the target. They represent carbon already committed to future structural
growth, not physical tissue that sets current NSC capacity. If the cap is not
sufficient for seasonal grasses, a future alternative should use a lagged or
annual-maximum tissue biomass rather than these committed-growth pools.

Explicit graminoid rhizomes are deferred. They are absent from upstream ELM
and from this reconstruction; the older ELM-Peatlands implementation reused
`livecroot` for rhizome C/N/P and changed allocation, respiration, phenology,
turnover, litter, and isotope pathways. Until that behavior is ported as a
separate sedge-focused change, nonwoody target capacity uses only current leaf
and fine-root carbon.

Non-peatland topounits retain the standard ELM equation exactly. The two
target fractions and maximum XR multiplier are optional during parameter-file
migration, with the values above used as fallbacks, but should be carried in
the standard ELM parameter file. They require evaluation against SPRUCE tissue NSC observations and
seasonal cpool dynamics; in particular, separate branch, bole, coarse-root,
leaf, and fine-root fractions would be preferable when sufficiently complete
data become available.

Primary empirical sources for these initial choices are:

- Furze et al. (2019), [Whole-tree nonstructural carbohydrate storage and
  seasonal dynamics in five temperate
  species](https://doi.org/10.1111/nph.15462).
- Richardson et al. (2015), [Distribution and mixing of old and new
  nonstructural carbon in two temperate
  trees](https://doi.org/10.1111/nph.13273).
- Schoonmaker et al. (2021), [Seasonal dynamics of non-structural carbon pools
  and their relationship to growth in two boreal conifer tree
  species](https://doi.org/10.1093/treephys/tpab013).

## Maintenance-respiration debt recovery

ELM's `xsmrpool` is a carbon-only maintenance-respiration reserve and debt
account. When maintenance respiration exceeds current photosynthesis, the pool
can become negative. Standard ELM repays that deficit from newly available
carbon over the parameterized `dayscrecover` timescale, giving repayment
priority over new growth.

On a topounit with positive `peat_depth`, existing nonstructural carbon in
`cpool` can also fund the repayment. The requested recovery rate remains

```text
-xsmrpool / dayscrecover.
```

At each timestep, no more than
`cpool_xsmr_recovery_max_frac * cpool` can come from the existing reserve. Any
unfunded remainder is taken from newly available carbon using the standard ELM
priority rule. The default maximum fraction is 0.1. Non-peatland topounits use
the standard ELM calculation exactly.

The transfer uses the existing `cpool_to_xsmrpool` flux, so carbon and carbon
isotopes remain conservative and the ecosystem carbon balance requires no new
source or sink. The maximum reserve-withdrawal fraction and the existing
`dayscrecover` parameter require scientific validation against nonstructural
carbon observations. `cpool_xsmr_recovery_max_frac` should be maintained in
the standard ELM parameter file; it is optional during migration so older
files use the documented default.
