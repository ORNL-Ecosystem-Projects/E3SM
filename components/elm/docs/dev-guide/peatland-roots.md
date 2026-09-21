# Peatland roots

`use_peatland_roots` is the model-level switch for peatland root behavior. It
defaults to true when `use_humhol` is enabled and false otherwise. Setting it
explicitly to false permits an attribution experiment with peatland hydrology
and standard ELM roots.

## Peatland root profile

When `use_peatland_roots` is enabled, vegetated topounits with positive peat depth use
a linear cumulative root profile originally calibrated from SPRUCE field
observations. The code calls this the **peatland root profile**; SPRUCE is its
calibration provenance, not a restriction on where it may be applied. Other
topounits retain ELM's standard Zeng (2001) profile.

For the peatland profile, `roota_par` is interpreted as the cumulative-profile
slope in 1/m and `rootb_par` as its dimensionless intercept. These meanings
differ from their meanings in the Zeng equation, so a peatland parameter file
must supply the corresponding values for peatland PFTs. The generated layer
fractions are checked during initialization. ELM aborts with the patch, PFT,
topounit, and complete profile if any fraction is negative or NaN, or if their
sum exceeds one. This is intended to expose a misapplied parameter set before
it produces opaque hydrology or carbon-balance failures.

`ROOTFR` is the profile used by root water uptake and plant-mediated methane
transport. ELM's vertically resolved biogeochemistry separately constructs
`FROOT_PROF`; the adaptive-N option below uses that BGC profile as its ordinary
root-access constraint. `FROOT_PROF` now also has an explicit nonnegative/NaN
check in addition to its existing unit-integral check.

## Adaptive mineral-N access

With native ELM vegetation and the RD nutrient competition model,
`use_peatland_roots` also enables experimental adaptive mineral-N access. It
only changes columns whose topounit has positive peat depth. Other nutrient
competition schemes and FATES still receive the peatland biophysical root
profile, but do not use this adaptive-N extension. The option is independent
of the methane backend: aqueous carbon and nutrient transport can motivate its
use, but plant nutrient acquisition remains a separate process.

For each vascular PFT, the ordinary fine-root profile remains a hard access
constraint. Demand within that profile is redistributed toward mineral N in
layers that are not completely water saturated:

```text
w(p,j) = root_profile(p,j) max(NH4(j) + NO3(j), 0),  if theta_liq(j) < porosity(j)
w(p,j) = 0,                                           otherwise
```

The weights are normalized over layer thickness, so the option redistributes
the PFT's existing N demand without creating demand or a second uptake flux.
The PFT-resolved demands are summed before the standard RD competition step.
The realized column uptake is then attributed back to PFTs using the same
layer fulfillment fractions. This construction conserves the column N sink.

Nonvascular PFTs are explicitly excluded from the adaptive weighting. Mosses
and lichens retain their prescribed fine-root-profile demand and therefore do
not gain access to deep mineral-N hotspots exposed by the vascular option.
Non-peat topounits retain standard ELM behavior.

Shrubs and trees currently use the same adaptive rule. Field behavior suggests
that shrubs may sustain nutrient acquisition from wet peat better than trees.
This lifeform contrast should be revisited after suitable observations are
identified; it should not be represented by an uncalibrated hard-coded PFT
multiplier. Candidate future formulations include lifeform-specific moisture
response curves, effective rooting-depth limits, and acclimation time scales.

The present binary saturation gate and direct mineral-N concentration
weighting should be treated as an upper-bound experiment. Before recommending
the option for production simulations, validate PFT-level uptake profiles,
soil mineral-N depletion, vegetation productivity, and sensitivity near full
saturation. The empty-eligible-layer fallback is the unchanged prescribed root
profile.
