# Peatland shading options

`use_shrub_moss_shading` applies a Beer-Lambert attenuation to the PAR received
by moss PFTs. Mosses and shrubs are identified from the `nonvascular` and
`woody` PFT metadata rather than fixed PFT numbers. Shrub LAI is weighted by
the actual patch fraction and only shrubs in the moss patch's own topounit and
column contribute, so radiation is not transferred between hummock, hollow,
and fen columns.

The option is false by default and therefore does not alter standard ELM or
peatland simulations unless explicitly enabled. The extinction coefficient is
0.5, matching the legacy ELM-Peatlands implementation. That coefficient still
requires scientific validation for the SPRUCE shrub and moss communities.

## Possible moss understory overlay

The current subgrid structure treats each PFT as a horizontally tiled patch.
Consequently, patch canopy fluxes are already multiplied by the moss patch
fraction when they are aggregated to the column. Simply replacing moss LAI by
`moss LAI * moss patch fraction` inside the existing canopy calculation would
apply that fraction a second time.

A future distributed-moss option should instead treat moss as an understory
overlay. It would calculate moss canopy exchange over the topounit using
column-area moss LAI and overstory-transmitted radiation, then convert the
result back to moss-patch-area fluxes before the standard patch-to-column
aggregation. Photosynthesis, canopy evaporation and energy exchange, and the
flux-to-carbon/water state updates must be handled together to keep C, water,
and energy conservation exact. This overlay is not implemented by
`use_shrub_moss_shading`.
