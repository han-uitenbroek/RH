PRO resample3datmos, atmos_in, STRIDE=stride, MAGNETIC=magnetic

  IF (n_params() LT 1 OR (NOT keyword_set(STRIDE))) THEN BEGIN
    print, "Usage: resample3datmos, atmos_in, [, STRIDE=stride] "
    return
  ENDIF 

  base_name = strmid(atmos_in, 0, strlen(atmos_in) - 6)

  a3d = read3datmos(atmos_in)

  a3dres = {Nx: a3d.Nx/stride,  Ny: a3d.Ny/stride,  $
            Nz: a3d.Nz,  NHydr: a3d.NHydr, boundary: a3d.boundary,  $
            dx: a3d.dx*stride,  dy: a3d.dy*stride,  z: a3d.z, $
            T: a3d.T[0:*:stride, 0:*:stride, *],  $
            n_elec: a3d.n_elec[0:*:stride, 0:*:stride, *],  $
            vturb: a3d.vturb[0:*:stride, 0:*:stride, *],  $
            vx: a3d.vx[0:*:stride, 0:*:stride, *],  $
            vy: a3d.vy[0:*:stride, 0:*:stride, *],  $
            vz: a3d.vz[0:*:stride, 0:*:stride, *],  $
            nH: a3d.nH[0:*:stride, 0:*:stride, *, *]}

  atmos_out = base_name + $
              string(FORMAT='("_stride", I2.2)', stride) + ".atmos"
  openw, lun, /XDR, /GET_LUN, atmos_out
  writeu, lun, a3dres
  free_lun, lun

  IF (keyword_set(MAGNETIC)) THEN BEGIN
    Bfile = base_name + '.B'
    openr, /XDR, /GET_LUN, lunB, Bfile
    B = dblarr(a3d.Nx, a3d.Ny, a3d.Nz)
    gamma = B
    chi = B
    readu, lunB, B, gamma, chi
    free_lun, lunB

    Bfile_res = base_name + string(FORMAT='("_stride", I2.2)', stride) + ".B"
    openw, /XDR, /GET_LUN, lunB, Bfile_res
    writeu, lunB, B[0:*:stride, 0:*:stride, *], $
            gamma[0:*:stride, 0:*:stride, *], $
            chi[0:*:stride, 0:*:stride, *]
    free_lun, lunB
  ENDIF
END
