PRO make2dslab, IRRAD=irrad

  ;; Write a 2-D slab for rhsc2d

  FIXED = 0L  &  PERIODIC = 1L
  IRRADIATED =  0L  &  ZERO = 1L  &  THERMALIZED = 2L

  lBoundCond = FIXED
  IF (keyword_set(IRRAD)) THEN $
   lBoundVal  = [IRRADIATED, IRRADIATED] $
  ELSE $
   lBoundVal  = [ZERO, ZERO]

  T_0  = 8.0E+3
  nH_0 = 1.0E+16

  n_elec_0 = 0.2 * nH_0
  vturb_0  = 2.0
  vx_0     = 0.0
  vz_0     = 0.0

  NHydr = 1L
  Nx = 50L  &  dx = 100.0 + dblarr(Nx)
  Nz = 50L  &  z  = reverse(100.0 * dindgen(Nz))
  IF (lBoundCond EQ FIXED) THEN dx[Nx-1] = 0.0

  T      = T_0 + dblarr(Nx, Nz)
  N_elec = n_elec_0 + dblarr(Nx, Nz)
  vturb  = vturb_0 + dblarr(Nx, Nz)
  vx     = dblarr(Nx, Nz)
  vz     = dblarr(Nx, Nz)
  nH     = nH_0 + dblarr(Nx, Nz, NHydr)

  openw, lun, /GET_LUN, 'slab_' + $
   string(FORMAT='(I3.3, "x", I3.3)', Nx, Nz) + $
   (keyword_set(IRRAD) ? '_irr' : '') + '.atmos', /XDR

;;-------------- Now write to file

  writeu, lun, Nx, Nz, NHydr
  writeu, lun, lBoundCond, lBoundVal
  writeu, lun, dx, z
  writeu, lun, T, n_elec, vturb, vx, vz, nH

  IF (keyword_set(IRRAD)) THEN BEGIN 

;; --- Read the radiation field from a 1-D model

    spawn, 'pwd', /NOSHELL, current_dir
    data_dir = '~/src/rh/rhf1d/run/'
    cd, data_dir

@files.common
@geometry.common
@spectrum.common

    metalFile = 'metals.out'
    moleculeFile = 'molecules.out'

    result = readInput('input.out')
    result = readgeometry('geometry.out')
    result = readatmos('atmos.out')
    atom   = readatom("atom.out")

    result = readspectrum('spectrum.out')
    ray = readray('spectrum_1.00')

;;    tabinv, spectrum.lambda, vacuumtoair(500.0), leff
;;    l500 = long(leff)
;;    [ray.i[0:l500-1], ray.i[l500+1:*]]
    irrad = ray.I
    
    irrad_left   = transpose(rebin(irrad, spectrum.Nspect, Nz)) * 0.1
    irrad_right  = transpose(rebin(irrad, spectrum.Nspect, Nz)) * 0.1
    irrad_top    = transpose(dblarr(spectrum.Nspect, Nx))
    irrad_bottom = transpose(rebin(irrad, spectrum.Nspect, Nx)) * 0.1

    writeu, lun, irrad_top, irrad_bottom, irrad_left, irrad_right

    cd, current_dir[0]
  ENDIF

  free_lun, lun
END
