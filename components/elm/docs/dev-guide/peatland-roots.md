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
`FROOT_PROF`; vascular mineral-N demand follows that BGC profile. `FROOT_PROF`
now also has an explicit nonnegative/NaN
check in addition to its existing unit-integral check.

## Vascular mineral-N access

With native ELM vegetation and the RD nutrient competition model,
`use_peatland_roots` resolves vascular demand using each PFT's prescribed
`FROOT_PROF`. It only changes columns whose topounit has positive peat depth.
Other nutrient competition schemes and FATES still receive the peatland
biophysical root profile, but do not use this PFT-resolved RD pathway. The
behavior is independent of the methane backend.

For each vascular PFT, the ordinary fine-root profile and the diagnosed main
water table are hard access constraints. Demand is zero in layers entirely
below ZWT. The layer intersecting ZWT is weighted by its geometric fraction
above the water table, and the existing local-saturation check continues to
exclude saturated or perched water within the nominally unsaturated interval:

```text
f_above_zwt(j) = fraction of layer j geometrically above ZWT
w(p,j) = root_profile(p,j) f_above_zwt(j),
          if theta_total(j) < porosity(j)
w(p,j) = 0, otherwise
```

The weights are normalized over layer thickness. They do not contain an NH4
or NO3 concentration factor, so vascular demand is not adaptively redirected
toward deep mineral-N hotspots. The PFT-resolved demands are summed before the
standard RD competition step. The realized column uptake is then attributed
back to PFTs using the same layer fulfillment fractions. This construction
conserves the column N sink.

Nonvascular PFTs are explicitly excluded from the vascular pathway. Mosses
and lichens retain their prescribed profile unless the separate capillary
connectivity option below is enabled. Non-peat topounits retain standard ELM
behavior.

## Experimental moss capillary nutrient access

`use_moss_capillary_nutrients` is a separate sensitivity option for native ELM
vegetation with RD nutrient competition. It defaults on with `use_humhol` and
can be disabled explicitly for controlled comparisons. When enabled, a
nonvascular PFT's ordinary mineral-N demand follows the peatland biophysical
`ROOTFR` profile rather than `FROOT_PROF`. For the current SPRUCE moss
parameterization, this confines ordinary access to approximately the upper
10 cm. It does not alter the profile used to place fine-root carbon turnover.

The option can redistribute part of that existing demand into layers between
the bottom of `ROOTFR` and the layer intersecting the main water table. Layers
below the water-table layer are not donors. Eligible-layer weights combine
NH4 plus NO3 concentration with a hydraulic-connectivity factor calculated
from the complete surface-to-layer path. Layer conductivity uses the same
Clapp-Hornberger liquid-saturation exponent and frozen-soil impedance used by
peatland hydrology:

```text
K(j) = HKSAT(j) Sliq(j)^(2 b(j) + 3) 10^(-e_ice f_ice(j))
tpath(j) = sum[k=1..j] dz(k) / K(k)
C(j) = exp(-tpath(j) / tau_cap)
```

The fraction assigned to the capillary profile is
`moss_capillary_max_demand_fraction` times the connectivity to the water-table
layer. The default maximum is 0.10 and the default connectivity time scale is
30 days. These two scalars are read from the standard ELM parameter file as
`moss_capillary_max_demand_fraction` and
`moss_capillary_connectivity_timescale_days`; they are not namelist options.
Both are explicitly experimental sensitivity values, not calibrated SPRUCE
parameters. A dry or frozen intervening layer increases series
resistance and suppresses access to every layer below it. Total PFT demand is
renormalized, so the option redistributes demand without creating additional
demand or bypassing standard plant-microbe competition.

This first implementation affects mineral N only. A mechanistic successor
should transport NH4, NO3, and solution P upward with capillary water and let
moss retain strictly surface-localized uptake.

Shrubs and trees currently use the same water-table clipping rule. Field
behavior suggests that shrubs may sustain nutrient acquisition from wet peat
better than trees. This lifeform contrast should be revisited after suitable
observations are identified; it should not be represented by an uncalibrated
hard-coded PFT multiplier. Candidate future formulations include
lifeform-specific moisture response curves, effective rooting-depth limits,
and acclimation time scales.
