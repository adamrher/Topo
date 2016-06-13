&topoparams
  grid_descriptor_fname           = 'inputdata/grid-descriptor-file/fv_0.9x1.25.nc'
  intermediate_cubed_sphere_fname = '../bin_to_cube/gmted2010_modis-ncube3000.nc'
  output_fname                    = 'output/$argv[1]-$argv[2]-$argv[3]-intermediate_ncube$argv[4]-$argv[5].nc'
  externally_smoothed_topo_file   = '../cam_fv_topo-smoothing/$argv[2]-$argv[1]-$argv[3].nc'
  lsmooth_terr                    = .true.
  lexternal_smooth_terr           = .true.
  lzero_out_ocean_point_phis      = .false.
  lsmooth_on_cubed_sphere         = .false.
  ncube_sph_smooth_coarse         = 20  
  ncube_sph_smooth_fine           = 1
  lfind_ridges                    = .false.
  nwindow_halfwidth               = 14
  nridge_subsample                = 14
/
