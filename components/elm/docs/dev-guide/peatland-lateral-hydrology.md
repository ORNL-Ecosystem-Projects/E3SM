# Peatland lateral hydrology

## Aquifer connectivity

`TopounitRegionalTarget` in the surface dataset defines the lateral aquifer
graph. A target creates one conservative, bidirectional exchange pair: the
instantaneous water flux can run in either direction according to the
hydraulic-head difference. Each pair should therefore appear only once in the
target array.

ELM honors these explicit pairs across bog/non-bog classifications. This
replaces the former hard-coded exception for the three-topounit SPRUCE layout
and makes site connectivity a surface-data decision.

When `use_humhol = .true.` and a surface dataset does not contain
`TopounitRegionalTarget`, ELM uses the following peatland ordered-layout
defaults:

| Layout | Topounit order | Target array | Aquifer pairs |
| --- | --- | --- | --- |
| Three units | fen, hollow, hummock | `0, 1, 2` | fen ↔ hollow ↔ hummock |
| Four units | fen, hollow, hummock, upland | `0, 0, 2, 1` | fen ↔ upland; hollow ↔ hummock |

Here, zero means no target and positive values are one-based local topounit
indices. Four-unit surface-data generators and existing input files should
write `0, 0, 2, 1` explicitly rather than rely on the fallback.

## Surface routing

Aquifer targets do not control surface or perched-water run-on.
`TopounitSurfaceTarget` defines a separate, one-way downhill graph. ELM checks
that each nonzero target is another active, lower-elevation topounit. When
`use_humhol = .true.` and the variable is absent, the peatland three- and
four-unit defaults are respectively `0, 1, 2` and `0, 1, 2, 1`.

For the standard four-unit layout, routing is

```text
hummock → hollow → fen
upland → fen
```

Thus the four-unit case has two disconnected groundwater exchange pairs but a
surface pathway from the hollow into the fen. Upland runoff bypasses the bog
and drains directly to the fen. IM2 hillslope hydrology continues to derive
its nearest-lower-neighbor graph from elevation rather than using the
peatland surface target.

Outside `use_humhol`, topounit count does not imply peatland classes. Missing
target variables retain the historical index-ordered fallback, and IM2 still
constructs the surface-routing graph from elevation.
