# Profile Assimilation Development Notes
Eliot Quon, Shashank Yellapantula

TODOs:
- [x] Implement general input for weighting matrix
- [x] Update z normaliation in regression matrix
- [ ] Combine inputs for partial profile assimilation (fix issue with netcdf https://github.com/ewquon/amr-wind/commit/d8db4842ad8b8c6c2d123b234ec1b9d979ef8718, get rid of workaround)
- [ ] Add error profile output (to netcdf)
- [ ] Fix netCDF input format so that mesoscale data have separate ntime,nheight dimensions
      rather than an ambiguous single arraySize dimension
- [ ] Generalize polynomial fit (i.e., implement `m_ind_polyOrder`)
- [x] Generalize "Wrf" naming in code to "Meso"
- [ ] Deal with assumption that the ABL statistics class computes statistics at the cell-centers
      only on level 0
- [x] Separate ABLWrfForcing class into a separate source file
- [ ] Add r-test


## Differences from SOWFA

- Regression matrix (Z^T W Z) is formed with z/zmax in SOWFA
  Regression matrix (Z^T W Z) is formed with z*scaleFact in AMR-Wind, which was previously
    hard-coded with scaleFact=1e-3; SOWFA behavior is achieved with the `normalize_by_zmax` flag
- Blending/forcing transition in SOWFA is over:
    (assimMaxHeight, assimMaxHeight+blendThickness)
  Blending/forcing transition in AMR-Wind is over:
    (transition_height, transition_height+transition_thickness)
- Blending (from indirect) to constant/direct is enabled in SOWFA with the `blendToConst`
    and `blendToDirect` flags, respectively
  Blending in AMR-Wind is enabled by setting `forcing_transition=indirectToConstant` or
    `forcing_transition=indirectToDirect`

