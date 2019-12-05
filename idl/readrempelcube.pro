FUNCTION readvarcube, N0, N1, N2, filename

  dummy = fltarr(N0, N1, N2)

  openr, lun, /GET_LUN, filename
  readu, lun, dummy
  free_lun, lun

  return, dummy
END

;; Convert Mathias' input files to RH input atmosphere

  outputFile   = "Rempel_80G_dynamo.atmos"
  B_outputFile = "Rempel_80G_dynamo.B"

  IRRADIATED = 0L  &  ZERO = 1L  &  PLANCK = 2L

  ;; Physical constants

  WGHT_PER_H   = 1.4271D0
  AVG_MOL_WGHT = 1.2981D0
  AMU          = 1.6605402D-27
  KBOLTZMANN   = 1.380658D-23

  ;; Unit conversions

  CM_TO_M   = 1.0D-2
  KM_TO_M   = 1.0D+3
  G_TO_KG   = 1.0D-3
  DYNE_TO_N = 1.0D-5
  GAUSS_TO_TESLA = 1.0D-4

  ;; Horizontal grid spacing [km]

  DX = 8.0D0
  DY = 8.0D0
  DZ = 8.0D0

  ;; Dimensions of the original cube

  N0 = 768L   ;; x
  N1 = 180L   ;; z
  N2 = 768L   ;; y

  KMIN  = 25L
  KMAX  = 154L
  NZNEW = KMAX - KMIN + 1

  v_files   = ['velx_051000.float', 'vely_051000.float', 'velz_051000.float']
  B_files   = ['magx_051000.float', 'magy_051000.float', 'magz_051000.float']
  T_file    = 'temp_051000.float'
  dens_file = 'dens_051000.float'

  vx = (reverse(transpose(readvarcube(N0, N1, N2, v_files[0]), [0, 2, 1]), 3, $
                /OVERWRITE))[*, *, KMIN:KMAX] * (CM_TO_M / KM_TO_M)
  vy = (reverse(transpose(readvarcube(N0, N1, N2, v_files[2]), [0, 2, 1]), 3, $
                /OVERWRITE))[*, *, KMIN:KMAX] * (CM_TO_M / KM_TO_M)
  vz = (reverse(transpose(readvarcube(N0, N1, N2, v_files[1]), [0, 2, 1]), 3, $
                /OVERWRITE))[*, *, KMIN:KMAX] * (CM_TO_M / KM_TO_M)

  Bx = (reverse(transpose(readvarcube(N0, N1, N2, B_files[0]), [0, 2, 1]), 3, $
                /OVERWRITE))[*, *, KMIN:KMAX] * GAUSS_TO_TESLA
  By = (reverse(transpose(readvarcube(N0, N1, N2, B_files[2]), [0, 2, 1]), 3, $
                /OVERWRITE))[*, *, KMIN:KMAX] * GAUSS_TO_TESLA
  Bz = (reverse(transpose(readvarcube(N0, N1, N2, B_files[1]), [0, 2, 1]), 3, $
                /OVERWRITE))[*, *, KMIN:KMAX] * GAUSS_TO_TESLA

  B     = sqrt(Bx^2 + By^2 + Bz^2)
  gamma = acos(Bz / B)
  chi   = atan(By, Bx)

  T = (reverse(transpose(readvarcube(N0, N1, N2, T_file[0]), [0, 2, 1]), 3, $
              /OVERWRITE))[*, *, KMIN:KMAX]
  rho = (reverse(transpose(readvarcube(N0, N1, N2, dens_file), $
                           [0, 2, 1]), 3, /OVERWRITE))[*, *, KMIN:KMAX]
  nHtot = (rho * G_TO_KG / CM_TO_M^3) / (WGHT_PER_H * AMU)


  openw, lun, /GET_LUN, /XDR, B_outputFile
  writeu, lun, B, gamma, chi
  free_lun, lun

  openw, lun, outputFile, /GET_LUN, /XDR

  lBoundVal = [ZERO,  PLANCK]

  NHYDR = 1L
  NX    = N0
  NY    = N2
  Nz    = NZNEW

  z  = reverse(dindgen(Nz) * DZ)

  writeu, lun, NX, NY, NZ, NHYDR
  writeu, lun, lBoundVal
  writeu, lun, DX, DY, z

  writeu, lun, double(T)
  writeu, lun, dblarr(Nx, Ny, Nz)     ;; n_elec
  writeu, lun, dblarr(Nx, Ny, Nz)     ;; vturb

  ;; Remember: second dimension in the original cube is the vertical

  writeu, lun, double(vx)
  writeu, lun, double(vy)
  writeu, lun, double(vz)

  writeu, lun, double(nHtot)

  free_lun, lun
END
