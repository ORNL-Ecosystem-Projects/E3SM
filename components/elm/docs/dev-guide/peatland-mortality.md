# Peatland shrub mortality

The peatland broadleaf deciduous boreal shrub uses its PFT-specific annual
mortality parameter, `r_mort`. This behavior is selected by the PFT name
`peatlnd_broadleaf_deciduous_boreal_shrub`; it does not require a separate
namelist option or topounit check because the PFT itself is peatland-specific.

All other PFTs retain standard ELM mortality behavior. In RD mode their annual
mortality continues to come from the soil-order calculation. In other modes
they continue to use the scalar `r_mort` value.

The parameter reader accepts both standard ELM files, where
`r_mort(allpfts)` contains one value, and peatland files, where `r_mort(pft)`
contains one value per PFT. With the scalar layout, the peatland shrub inherits
the standard scalar value. With the PFT layout, only the peatland shrub's entry
overrides the normal mortality path. This keeps peatlands-off configurations
bit-for-bit compatible while allowing the shrub mortality rate to be calibrated
independently in the standard ELM parameter file.
