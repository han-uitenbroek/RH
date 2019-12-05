; -------- file: -------------------------- readray.pro -------------- ;

; -------- begin -------------------------- readRay.pro -------------- ;

FUNCTION readRay, fileName, FORCE_STOKES=force_stokes

@input.common
@atmos.common
@geometry.common
@spectrum.common

  WHILE (NOT existFile(fileName, UNIT=rayUnit, /XDR)) DO BEGIN
   answer = dialog_message(/QUESTION, "Find new ray file file?")
   IF (answer EQ 'Yes') THEN BEGIN
    fileName = dialog_pickfile(FILTER='*_?.??*', $
                               TITLE='Spectroscopic data file', $
                               /MUST_EXIST, /READ, FILE=fileName)
    IF (fileName EQ '') THEN return, 0
   ENDIF ELSE $
    return, 0

  ENDWHILE

  ;; --- Structure definitions --                       -------------- ;

  CASE geometryType OF 
    "ONE_D_PLANE": BEGIN
      muz = 0.0D+0
      readu, rayUnit, muz
      I = dblarr(spectrum.Nspect)
    END
    "TWO_D_PLANE": BEGIN
      mux = 0.0D+0
      muz = 0.0D+0
      readu, rayUnit, mux, muz
      I = dblarr(geometry.Nx, spectrum.Nspect)
    END
    "THREE_D_PLANE": BEGIN
      mux = 0.0D+0
      muy = 0.0D+0
      readu, rayUnit, mux, muy
      I = dblarr(geometry.Nx, geometry.Ny, spectrum.Nspect)
    END
    "SPHERICAL_SYMMETRIC": BEGIN
      muz = 0.0D+0
      readu, rayUnit, muz
      I = dblarr(spectrum.Nspect)
    END
  ENDCASE

  ;; --- Read in one go --                              -------------- ;

  Nspect = 0L
  readu, rayUnit, I, Nspect

  CASE (geometryType) OF
    "ONE_D_PLANE": BEGIN
      IF (nspect GT 0) THEN $
        d = {nspect: 0L, chi: dblarr(geometry.Ndep), S: dblarr(geometry.Ndep)}
    END
    "TWO_D_PLANE": BEGIN
      IF (nspect GT 0) THEN $
        d = {nspect: 0L, chi: dblarr(geometry.Nx, geometry.Nz), $
             S: dblarr(geometry.Nx, geometry.Nz)}
    END
    "THREE_D_PLANE": BEGIN
      IF (nspect GT 0) THEN $
        d = {nspect: 0L, chi: dblarr(geometry.Nx, geometry.Ny, geometry.Nz), $
             S: dblarr(geometry.Nx, geometry.Ny, geometry.Nz)}
    END
    "SPHERICAL_SYMMETRIC": BEGIN
      IF (nspect GT 0) THEN $
        d = {nspect: 0L, chi: dblarr(geometry.Nradius), $
             S: dblarr(geometry.Nradius)}
    END
  ENDCASE

  IF (Nspect GT 0) THEN BEGIN
    dr = replicate(d, Nspect)
    readu, rayUnit, dr
  ENDIF

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    IF (Nspect GT 0) THEN $
     ray = {mux: mux, muz: muz, I: I, $
            nspect: dr.nspect, chi: dr.chi, S: dr.S} $
    ELSE $
     ray = {mux: mux, muz: muz, I: I}
  ENDIF
  IF (geometryType EQ "THREE_D_PLANE") THEN BEGIN
    IF (Nspect GT 0) THEN $
     ray = {mux: mux, muy: muy, I: I, $
            nspect: dr.nspect, chi: dr.chi, S: dr.S} $
    ELSE $
     ray = {mux: mux, muy: muy, I: I}
  ENDIF

  IF (geometryType EQ "ONE_D_PLANE" OR $
      geometryType EQ "SPHERICAL_SYMMETRIC") THEN BEGIN
    IF (Nspect GT 0) THEN $
     ray = {muz: muz, I: I, nspect: dr.nspect, chi: dr.chi, S: dr.S} $
    ELSE $
     ray = {muz: muz, I: I}
  ENDIF

  IF (tag_present(atmos, 'STOKES') OR inputData.backgr_pol OR $
      keyword_set(FORCE_STOKES)) THEN BEGIN
    CASE (geometryType) OF
      "ONE_D_PLANE": BEGIN
        Stokes_Q = dblarr(spectrum.Nspect)
      END
      "TWO_D_PLANE": BEGIN
        Stokes_Q = dblarr(geometry.Nx, spectrum.Nspect)
      END
      "THREE_D_PLANE": BEGIN
        Stokes_Q = dblarr(geometry.Nx, geometry.Ny, spectrum.Nspect)
      END
      ELSE:
    ENDCASE
    Stokes_U = Stokes_Q
    Stokes_V = Stokes_Q
    readu, rayUnit, stokes_Q, stokes_U, stokes_V
    ray = create_struct(ray, 'stokes_Q', Stokes_Q, $
                        'stokes_U', Stokes_U, 'stokes_V', Stokes_V)
  ENDIF

  free_lun, rayUnit
  return, ray
END
; -------- end ---------------------------- readRay.pro -------------- ;
