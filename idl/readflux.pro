; -------- file: -------------------------- readflux.pro ------------- ;

; -------- begin -------------------------- readFlux.pro ------------- ;

FUNCTION readFlux, fileName

@geometry.common
@spectrum.common

  WHILE (NOT existFile(fileName, UNIT=fluxUnit, /XDR)) DO BEGIN
   answer = dialog_message(/QUESTION, "Find new flux file?")
   IF (answer EQ 'Yes') THEN BEGIN
    fileName = dialog_pickfile(FILTER='*.out', TITLE='Flux output file', $
                               /MUST_EXIST, /READ, FILE=fileName)
    IF (fileName EQ '') THEN return, 0
   ENDIF ELSE $
    return, 0

  ENDWHILE

  ;; --- Structure definitions --                       -------------- ;

  CASE geometryType OF 
    "ONE_D_PLANE": $
     flux = dblarr(spectrum.Nspect)
    "TWO_D_PLANE": $
     flux = dblarr(geometry.Nx, spectrum.Nspect)
    "THREE_D_PLANE": $
     flux = dblarr(geometry.Nx * geometry.Ny, spectrum.Nspect)
    "SPHERICAL_SYMMETRIC": $
     flux = dblarr(spectrum.Nspect)
  ENDCASE

  readu, fluxUnit, flux

  ;; --- Transpose indices so that spectral index runs fastest, and
  ;;     spatial index runs slowest --                  -------------- ;

  CASE (geometryType) OF
    "ONE_D_PLANE":
    "TWO_D_PLANE": $
      flux = transpose(flux)
    "THREE_D_PLANE": $
     flux = reform(transpose(flux), spectrum.Nspect, $
                   geometry.Nx, geometry.Ny)
    "SPHERICAL_SYMMETRIC":
  ENDCASE

  free_lun, fluxUnit
  return, 1
END
; -------- end ---------------------------- readSpectrum.pro --------- ;

