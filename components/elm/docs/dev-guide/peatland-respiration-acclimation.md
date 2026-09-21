# Peatland respiration acclimation

Maintenance respiration retains the standard ELM formulation outside peatland
topounits. When `use_humhol` is enabled and a patch belongs to a topounit with
positive peat depth, roots and live wood use PFT-specific base respiration and
temperature-sensitivity parameters. C3 leaf respiration uses the same
PFT-specific base rate and a Q10 temperature response only on those peatland
topounits; non-peatland patches continue to use the standard ELM leaf
respiration temperature function.

Woody maintenance respiration can acclimate to warming relative to an annual
mean reference temperature collected during spinup. The reference is the
running mean of up to `mr_acclim_spinup_years` annual temperatures. After the
reference is established, the effective woody base rate is adjusted by the
configured fraction of the modeled warming. Fine-root respiration uses the
PFT-specific parameters but is not acclimated by this term.

## Parameters

All parameters below are read from the ELM parameter file. The new parameters
will be migrated into the standard ELM parameter file so that custom peatland
files are not required solely to define respiration behavior.

| Parameter | Scope | Units | Missing-field default | Purpose |
| --- | --- | --- | --- | --- |
| `br_mr_pft` | PFT | g C g N-1 s-1 | 2.525e-6 | Maintenance-respiration base rate at 20 C |
| `q10_mr_pft` | PFT | dimensionless | 1.5 | Multiplicative respiration response per 10 K |
| `mr_acclim_warming_frac` | scalar | dimensionless | 0 | Fraction of warming compensated by woody acclimation; zero disables acclimation |
| `mr_acclim_spinup_years` | scalar | years | 0 | Maximum spinup years used to construct the baseline; zero disables baseline collection |

The fallback values reproduce the standard ELM maintenance-respiration
constants. The PFT-specific values and acclimation controls require scientific
validation for each peatland PFT before production simulations.

`SPINUP_T` exposes the stored reference temperature in history output. The
reference temperature and the number of contributing years are restart state.
