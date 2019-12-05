PRO ATMOSatlas, lambdaMin, lambdaMax, Iatlas, latlas

  ;+ Reads selected wavelength range from the ATMOS infrared Solar Atlas
  ;
  ; Parameters:
  ;   Input:    lambdaMin  -- Minimum wavelength [nm]
  ;             lambdaMax  -- Maximum wavelength [nm]
  ;
  ;   Output:   Iatlas     -- Integer array of intensities normalized to 10^4
  ;             latlas     -- Float array of corresponding wavelengths [nm]
  ;
  ;-

  if (n_params(0) lt 3) then begin
    print, "Usage: ATMOSatlas, lambdaMin, lambdaMax, Iatlas, [latlas]"
    print, "NOTE:  Wavelengths are in nanometers!!"
    return
  ENDIF
  IF (lambdaMin GE lambdaMax) THEN BEGIN
    print, "Minimum wavelength should be larger than maximum"
    return
  ENDIF

  file = 'ATMOS.atlas'
  atlasdir = getenv('RH_ATLAS_PATH') + '/ATMOS/'
  openr, unit, atlasdir + file, /GET_LUN

  IF (NOT is_big_endian()) THEN $
   little_endian = 1 $
  ELSE $
   little_endian = 0

  Nlambda = 0L
  readu, unit, Nlambda
  IF (little_endian) THEN Nlambda = swap_endian(Nlambda)

  latlas = dblarr(Nlambda)
  Iatlas = fltarr(Nlambda)
  readu, unit, latlas, Iatlas
  IF (little_endian) THEN BEGIN
    latlas = swap_endian(latlas)
    Iatlas = swap_endian(Iatlas)
  ENDIF
  free_lun, unit

  index = where(latlas GE lambdaMin  AND  latlas LE lambdaMax, count)
  IF (count GT 0) THEN BEGIN
    latlas = latlas[index]
    Iatlas = Iatlas[index]
  ENDIF ELSE BEGIN
    print, FORMAT='(" Wavelength outside range ["' + $
     ', F9.3, " ,", F9.3, "] of atlas ")', latlas[0], latlas[Nlambda-1]
    latlas = [0.0]
    Iatlas = [0.0]
  ENDELSE

  return
end
