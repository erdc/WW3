cp ww3_grid_a.nml ww3_grid.nml
ww3_grid
cp ww3_prnc_level.nml ww3_prnc.nml
ww3_prnc
cp ww3_prnc_wind.nml ww3_prnc.nml
ww3_prnc
cp ww3_prnc_current.nml ww3_prnc.nml
ww3_prnc
mpiexec -n 4 ww3_shel

