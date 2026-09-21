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

## Surface-structure ground shading

`use_surface_structure_shading` applies optional, topounit-resolved shading by
boardwalks or other surface structures. Two optional variables configure each
topounit in the surface dataset:

* `TopounitStructureShadeFrac`: fraction of the topounit covered by the
  structure, from zero to one.
* `TopounitStructureLightTrans`: shortwave transmissivity through the covered
  area, from zero to one.

If either variable is absent, its no-effect default is used: zero covered
fraction and unit transmission. The effective multiplier on ground shortwave
absorption is

```
1 - TopounitStructureShadeFrac * (1 - TopounitStructureLightTrans)
```

The reduction is applied consistently to soil, snow, snow-layer, and aerosol
forcing absorption diagnostics. The removed ground absorption is assigned to
reflected shortwave so the modeled land-surface energy budget remains closed.
The option is false by default, so adding the optional surface properties does
not activate the physics by itself.

OLMT generates both properties with their no-effect defaults for topounit
cases. An active configuration can override individual values through
`[surface_data]`; each comma-separated triplet is a zero-based topounit index,
zero-based gridcell index, and value:

```ini
[case_options]
use_surface_structure_shading = .true.

[surface_data]
TopounitStructureShadeFrac = 0, 0, 1.0, 1, 0, 0.0, 2, 0, 0.0
TopounitStructureLightTrans = 0, 0, 0.5, 1, 0, 1.0, 2, 0, 1.0
```

This example gives topounit 0 a 50% effective reduction in ground shortwave
absorption and leaves topounits 1 and 2 unchanged. The companion SPRUCE OLMT
configuration uses this layout for its boardwalk/fen, hollow, and hummock
topounits, respectively.

This mechanism is more general than a hard-coded SPRUCE boardwalk fraction and
can support controlled experiments with other structures, including solar
panels. It is only a ground-shading approximation, however. It does not yet
alter canopy radiation, partition intercepted energy among reflection,
structure heating, and electrical export, calculate a structure temperature,
or redistribute precipitation. Those processes require a separate structure
energy and water balance and should not be inferred from this option.
