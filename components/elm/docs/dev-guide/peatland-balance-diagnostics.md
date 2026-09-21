# Peatland balance diagnostics

The standard ELM water and C/N/P conservation checks remain authoritative for
peatland simulations. Peatland lateral aquifer exchange is included in the
column water balance as a signed flux: positive values add water to a column
and negative values remove it. Surface routing and the existing generic
lateral-flow terms remain separate entries in the same balance.

When a fatal water-balance error occurs with `use_humhol` enabled, ELM reports
the failing topounit's weight, elevation, peat depth, lateral distance, regional
target, bog flag, and uphill water store. It also reports surface, snow,
aquifer, and layer-resolved soil water states together with a reconstruction
of the complete column flux sum. These values are diagnostic only; the balance
threshold and model state are not changed by this feature.

The carbon check snapshots vegetation, coarse woody debris, litter, and SOM at
the beginning and end of each timestep. On failure it reports the component
changes, phenology litter source and sink, and patch-level carbon stocks,
fluxes, and root profiles. The existing nitrogen and phosphorus checks already
report their complete external inputs, outputs, and storage changes.

This implementation intentionally does not retain the historical Peatlands
branch's relaxed water-balance abort threshold. Peatland simulations therefore
have the same numerical conservation tolerance as standard ELM. It also does
not move or truncate any lateral-routing state from inside a diagnostic path.

