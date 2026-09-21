# Peatland carbon reserves

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
