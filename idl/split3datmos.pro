PRO split3datmos, atmos_in, ROW=row, COLUMN=column, MAGNETIC=magnetic

  IF (n_params() LT 1 OR $
      ((NOT keyword_set(ROW)) AND (NOT keyword_set(COLUMN))) OR $
      (keyword_set(ROW) AND keyword_set(COLUMN))) THEN BEGIN
    print, "Usage: split3datmos, atmos_in, atmos_out [, ROW=row] "
    print, "        [, COLUMN=column] [, /MAGNETIC]"
    return
  ENDIF 

  a3d = read3datmos(atmos_in)
  PERIODIC = 1L

  base_name = strmid(atmos_in, 0, strlen(atmos_in) - 6)

  IF (keyword_set(MAGNETIC)) THEN BEGIN
    Bfile = base_name + '.B'
    openr, /XDR, /GET_LUN, lunB, Bfile
    B = dblarr(a3d.Nx, a3d.Ny, a3d.Nz)
    gamma = B
    chi = B
    readu, lunB, B, gamma, chi
    free_lun, lunB
  ENDIF

  IF (keyword_set(ROW)) THEN BEGIN
    IF (row LT 0 OR row GE a3d.Ny) THEN BEGIN
      print, $
       string(FORMAT='("ROW must be between 0 and Ny-1 = ", I3)', a3d.Ny-1)
      return
    ENDIF

    a2d = {Nx: a3d.Nx, Nz: a3d.Nz, NHydr: a3d.NHydr, $
           boundary: [PERIODIC, a3d.boundary], $
           dx: dblarr(a3d.Nx) + a3d.dx, z: a3d.z, $
           T: reform(a3d.T[*, row, *]), $
           n_elec: reform(a3d.n_elec[*, row, *]), $
           vturb: reform(a3d.vturb[*, row, *]), $
           vx: reform(a3d.vx[*, row, *]), vz: reform(a3d.vz[*, row, *]), $
           nH: ((a3d.NHydr GT 1) ? reform(a3d.nH[*, row, *, *]) : $
                                   reform(a3d.nH[*, row, *]))}

    atmos_out = base_name + $
     string(FORMAT='("_y", I3.3)', row) + ".atmos"

    IF (keyword_set(MAGNETIC)) THEN BEGIN
      Bfile_2D = base_name + string(FORMAT='("_y", I3.3)', row) + ".B"
      openw, /XDR, /GET_LUN, lunB, Bfile_2D
      writeu, lunB, B[*, row, *], gamma[*, row, *], chi[*, row, *]
      free_lun, lunB
    ENDIF
  ENDIF

  IF (keyword_set(COLUMN)) THEN BEGIN
    IF (column LT 0 OR column GE a3d.Nx) THEN BEGIN
      print, $
       string(FORMAT='("COLUMN must be between 0 and Nx-1 = ", I3)', ad3.Nx-1)
      return
    ENDIF

    a2d = {Nx: a3d.Ny, Nz: a3d.Nz, NHydr: a3d.NHydr, $
           boundary: [PERIODIC, a3d.boundary], $
           dx: dblarr(a3d.Ny) + a3d.dy, z: a3d.z, $
           T: reform(a3d.T[column, *, *]), $
           n_elec: reform(a3d.n_elec[column, *, *]), $
           vturb: reform(a3d.vturb[column, *, *]), $
           vx: reform(a3d.vy[column, *, *]), $
           vz: reform(a3d.vz[column, *, *]), $
           nH: ((a3d.Nhydr GT 1) ? reform(a3d.nH[column, *, *, *]) : $
                                   reform(a3d.nH[column, *, *]))}

    atmos_out = strmid(atmos_in, 0, strlen(atmos_in) - 6) + $
     string(FORMAT='("_x", I3.3)', column) + ".atmos"

    IF (keyword_set(MAGNETIC)) THEN BEGIN
      Bfile_2D = base_name + string(FORMAT='("_y", I3.3)', row) + ".B"
      openw, /XDR, /GET_LUN, lunB, Bfile_2D
      writeu, lunB, $
              B[*, *, column], gamma[*, *, column], chi[*, *, column]
      free_lun, lunB
    ENDIF
  ENDIF

  openw, lun, /XDR, /GET_LUN, atmos_out
  writeu, lun, a2d
  free_lun, lun
END
