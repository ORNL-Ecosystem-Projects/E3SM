# Peatland moss physiology

The peatland moss implementation is active only when `use_humhol` is true and
the PFT has `nonvascular = 1`. Standard ELM behavior is retained otherwise.
The implementation is generalized from the historical ELM-Peatlands/SPRUCE
code: it identifies mosses from PFT properties rather than a hard-coded PFT
number, and it uses the model's atmospheric N-deposition forcing rather than a
site-specific constant.

## Water, photosynthesis, and respiration

Moss internal water is diagnosed from mean liquid volumetric water content in
soil layers 3 and 4. The diagnosed volumetric water content is limited to
0--0.25 before applying the historical polynomial relation. Intercepted canopy
water supplies an external store. `H2O_MOSS_INTER` and `H2O_MOSS_WC` expose the
two diagnosed quantities on history files.

Nonvascular photosynthesis uses the historical moss surface-conductance
polynomial, tissue-water limitation, and VPD limitation instead of vascular
root `btran`. The minimum VPD scalar is 0.1 rather than zero. Moss also uses its
own photosynthetic temperature-acclimation coefficients. PFT-specific
maintenance respiration and woody respiration acclimation are documented and
implemented separately from the core moss physiology.

Snow and surface water reduce the exposed fraction of nonvascular vegetation.
Nonvascular PFTs with liquid water available at the surface bypass the
vascular root-resistance calculation; this does not give them adaptive deep-N
access.

## Nitrogen deposition

A Beer-law interception fraction is diagnosed from column-weighted exposed
nonvascular LAI using a fixed coefficient of 1.2. The intercepted fraction of
the actual `forc_ndep` flux is sent directly to the nonvascular plant N pool;
the remainder enters soil mineral N through the standard pathway. Allocation
among multiple moss PFTs is LAI weighted and exactly conservative after patch
weights are applied. `NDEP_TO_NPOOL` reports the patch-level intercepted flux.

## Parameters

The implementation reads the following variables from the ELM PFT parameter
file. They should be migrated into, documented in, and maintained with the
standard ELM parameter file rather than kept in a site-private file.

| Parameter | Scope | Units | Default or requirement | Purpose |
| --- | --- | --- | --- | --- |
| `vpd_min_moss` | scalar | Pa | 1000 | VPD below which moss is unstressed |
| `vpd_max_moss` | scalar | Pa | 1500 | VPD at which the moss VPD scalar reaches 0.1 |
| `vwc_moss_offset` | scalar | m3 m-3 | 0 | Offset subtracted from the layer-3/4 moss water proxy |
| `blower_u0` | scalar | m s-1 | 0 (disabled) | Surface wind contribution from the SPRUCE blower treatment |
| `blower_lambda` | scalar | m | 2 | E-folding height of the blower wind contribution |

`vpd_max_moss` must exceed `vpd_min_moss`, and `blower_lambda` must be
positive. The blower contribution is additionally gated by `use_humhol`, a
positive `blower_u0`, and simulation year 2015 or later.

## Intentional differences and deferred work

- Atmospheric N deposition uses `forc_ndep`; the historical hard-coded
  0.57 g N m-2 yr-1 value is not retained.
- PFT roles use `nonvascular`, `crop`, `iscft`, and PFT names rather than fixed
  numeric indices.
- Moss-water polynomial input is bounded below at zero to prevent unphysical
  extrapolation.
- Shrub-over-moss radiation shading is intentionally deferred to a separate
  commit and is not part of this implementation.
- Satellite-phenology modifications from the historical branch are excluded.
- The empirical SPRUCE tree/shrub phenology and transfer-pool conservation
  changes are kept in a separate commit.

Scientific validation should compare moss water content, surface conductance,
GPP, respiration, and N interception against observations. The VPD thresholds,
water-content offset, layer-3/4 proxy, fixed N-interception coefficient, and
blower profile remain calibration targets. The hydraulic-stress
photosynthesis path currently receives the moss temperature response but not
the empirical moss water/VPD scalar; that path needs a dedicated formulation
and test before it should be used for moss.
