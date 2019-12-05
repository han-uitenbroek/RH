FUNCTION readfullJ, L0=l0, LN=lN, J20

  DEFAULT_J_FILE   = "J.dat"
  DEFAULT_J20_FILE = "J20.out"

@geometry.common
@spectrum.common

  IF (NOT keyword_set(L0)) THEN l0 = 0
  IF (NOT keyword_set(LN)) THEN lN = spectrum.Nspect-1
  Nspect = lN - l0 + 1

  CASE geometryType OF
    "ONE_D_PLANE": BEGIN
      Jnu = dblarr(geometry.Ndep, Nspect)
      offset = 8L*l0 * geometry.Ndep
    END
    "TWO_D_PLANE": BEGIN
      Jnu = dblarr(geometry.Nx, geometry.Nz, Nspect)
      offset = 8L*l0 * (geometry.Nx*geometry.Nz)
    END
    "THREE_D_PLANE": BEGIN
      Jnu = dblarr(geometry.Nx, geometry.Ny, geometry.Nz, Nspect)
      offset = 8L*l0 * (geometry.Nx*geometry.Ny*geometry.Nz)
    END
    "SPHERICAL_SYMMETRIC": BEGIN
      Jnu = dblarr(geometry.Nradius, Nspect)
      offset = 8L*l0 * geometry.Nradius
    END
  ENDCASE

  openr, lun, /GET_LUN, DEFAULT_J_FILE
  point_lun, lun, offset
  readu, lun, Jnu
  free_lun, lun

  IF (n_params() EQ 1) THEN BEGIN
    J20 = Jnu
    openr, lun, /GET_LUN, DEFAULT_J20_FILE
    point_lun, lun, offset
    readu, lun, J20
    free_lun, lun
  ENDIF

  return, Jnu
END

