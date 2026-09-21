# Peatland phenology

The peatland phenology extensions preserve standard ELM behavior outside
peatland topounits. The SPRUCE-derived seasonal rules are used only when
`use_humhol` is enabled, the patch belongs to a topounit with positive peat
depth, and its PFT name identifies either the peatland needleleaf deciduous
boreal tree or peatland broadleaf deciduous boreal shrub. PFT names are resolved
at initialization rather than assuming fixed numeric PFT indices.

For those two PFTs, spring onset uses 2-m air temperature and accumulated
chilling. Autumn offset uses a combined air-temperature and daylength index.
Other seasonal-deciduous PFTs retain the standard soil-temperature and critical
daylength formulation. The generic onset threshold remains parameterized by
`crit_gdd1` and `crit_gdd2`; absent fields default to the standard ELM values
4.8 and 0.13, respectively. These parameters should ultimately be included in
the standard ELM parameter file for every PFT.

The chilling and autumn indices are initialized for new patches, written to
restart files, and exposed as `ONSET_CHIL` and `DAYL_TEMP` history fields.

## Transfer-pool conservation

Phenology counters can overshoot zero when their duration is not an exact
multiple of the model timestep. Cleanup therefore tests for counters less than
or equal to zero. At the final onset timestep, residual leaf, fine-root,
coarse-root, stem C, N, and P transfer pools are moved into displayed biomass
instead of being reset and lost. A warning reports nontrivial residuals so that
large timestep sensitivity remains visible. The same behavior is applied to
seasonal- and stress-deciduous phenology.

Crop detection in the BeTR phenology path uses the PFT `crop` and `iscft`
properties instead of assuming that all crop indices form a fixed numeric
range.

Satellite phenology changes are intentionally excluded from this commit.
