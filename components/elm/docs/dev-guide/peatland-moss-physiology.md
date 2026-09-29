# Peatland moss physiology

The peatland moss implementation is active only when `use_humhol` is true and
the PFT has `nonvascular = 1`. Standard ELM behavior is retained otherwise.
The implementation is generalized from the historical ELM-Peatlands/SPRUCE
code: it identifies mosses from PFT properties rather than a hard-coded PFT
number, and it uses the model's atmospheric N-deposition forcing rather than a
site-specific constant.

## Water, photosynthesis, and respiration

`use_prognostic_moss_water` selects a conserved living-moss water store and
defaults true when `use_humhol` is true. Water capacity scales with live moss
dry mass rather than a fixed soil pore volume. Live carbon is converted to dry
mass using the parameterized carbon fraction. A Clapp--Hornberger retention
relation maps the dry-mass water ratio to effective saturation and matric
potential, while a separately parameterized power law maps saturation to
hydraulic conductivity. The living-layer thickness is only the hydraulic path
length for exchange with the top peat layer. Atmospheric moss water loss is
removed from this store and is not also imposed as a vascular root sink.
Below the structural holding capacity, soil--moss exchange is one-way: soil
water may recharge the moss, but tightly held internal water does not drain
back into the soil in response to a downward matric-potential gradient. Only
water exceeding the drainage threshold is returned conservatively to the
upper peat.

CanopyFluxes receives the stored water divided by live moss dry mass in
g H2O g-1 dry mass, matching the units of the historical Sphagnum water
response. The cold start uses 10 gC m-2 of live moss, a carbon fraction of
0.5 gC g-1 dry mass, and 5 g H2O g-1 dry mass, giving 0.10 kg H2O m-2. At
200 gC m-2, the default maximum of 8 g g-1 permits 3.2 kg H2O m-2.
`H2O_MOSS_STORAGE`, `MOSS_WATER_POTENTIAL`, and `QFLX_MOSS_SOIL` expose the
state and positive soil-to-moss exchange.

Setting `use_prognostic_moss_water=.false.` retains the historical diagnostic
mapping from mean liquid volumetric water content in soil layers 3 and 4. Its
volumetric proxy remains limited to 0--0.25 before applying the historical
polynomial relation. This branch is retained for controlled attribution.

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
| `moss_water_layer_thickness` | scalar | m | 0.05 | Hydraulic path length between living moss and upper peat; does not set capacity |
| `moss_water_saturated_suction` | scalar | mm | 10 | Magnitude of matric potential at saturation |
| `moss_water_clapp_hornberger_b` | scalar | 1 | 3.5 | Water-retention exponent |
| `moss_water_hydraulic_conductivity_sat` | scalar | mm s-1 | 1.5 | Saturated living-moss hydraulic conductivity |
| `moss_water_conductivity_exponent` | scalar | 1 | 5 | Effective-saturation exponent for conductivity |
| `moss_initial_carbon` | scalar | g C m-2 | 10 | Cold-start live Sphagnum seed carbon |
| `moss_carbon_fraction_dry_mass` | scalar | g C g-1 dry mass | 0.5 | Converts prognostic live moss carbon to dry biomass |
| `moss_water_initial_content` | scalar | g H2O g-1 dry mass | 5 | Cold-start dry-mass water ratio; gives 0.10 kg H2O m-2 with the default seed |
| `moss_water_content_min` | scalar | g H2O g-1 dry mass | 0.5 | Dry-mass water ratio mapped to zero effective saturation |
| `moss_water_content_max` | scalar | g H2O g-1 dry mass | 8 | Maximum dry-mass water ratio; gives 3.2 kg H2O m-2 at 200 gC m-2 |
| `moss_water_drainage_threshold` | scalar | g H2O g-1 dry mass | 8.5 | Structural holding capacity; only water above this ratio drains from moss to upper peat |

`vpd_max_moss` must exceed `vpd_min_moss`, and `blower_lambda` must be
positive. The blower contribution is additionally gated by `use_humhol`, a
positive `blower_u0`, and simulation year 2015 or later.

## Intentional differences and deferred work

- Atmospheric N deposition uses `forc_ndep`; the historical hard-coded
  0.57 g N m-2 yr-1 value is not retained.
- PFT roles use `nonvascular`, `crop`, `iscft`, and PFT names rather than fixed
  numeric indices.
- Moss-water polynomial input is bounded below at zero to prevent unphysical
  extrapolation in the retained legacy branch.
- The retention and conductivity exponents are intentionally independent.
  Forcing the mineral-soil `2b+3` conductivity relation would make the fitted
  living-moss conductivity collapse too rapidly as it dries.
- Shrub-over-moss radiation shading is intentionally deferred to a separate
  commit and is not part of this implementation.
- Satellite-phenology modifications from the historical branch are excluded.
- The empirical SPRUCE tree/shrub phenology and transfer-pool conservation
  changes are kept in a separate commit.

Scientific validation should compare moss water content, matric potential,
soil-to-moss flux, surface conductance, GPP, respiration, and N interception
against observations. The initial hydraulic curve is a shared Sphagnum
hypothesis; hummock and hollow retention/conductivity curves, hydraulic path
length, initial dry-mass water ratio, and the physiological g/g mapping all
require calibration.
The VPD thresholds, legacy water-content offset and layer-3/4 proxy, fixed
N-interception coefficient, moss carbon fraction, and blower profile remain calibration targets. The hydraulic-stress
photosynthesis path currently receives the moss temperature response but not
the empirical moss water/VPD scalar; that path needs a dedicated formulation
and test before it should be used for moss.
