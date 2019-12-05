PRO KPspotatlas, lambdaMin, lambdaMax, Iatlas, latlas, COLD=COLD

  ;; Returns selected wavelength range from the Kitt Peak spot atlas.
  ;;
  ;; Parameters:
  ;;   Input:    lambdaMin  -- Minimum wavelength [nm]
  ;;             lambdaMax  -- Maximum wavelength [nm]

  ;;   Output:   Iatlas     -- Float array of Intensty values 
  ;;             latlas     -- Float array of corresponding wavelengths [nm]

  if (n_params(0) lt 1) then begin
    print, "Usage: KPspotatlas, lambdaMin, lambdaMax, Iatlas[, latlas], [/COLD]"
    print, "NOTE:  Wavelengths are in nanometers!!"
    return
  ENDIF
 
  IF (keyword_set(COLD))THEN $
   atlasfile = getenv('RH_ATLAS_PATH') + '/KPspot/KPspot_cold.atlas' $
  ELSE $
   atlasfile = getenv('RH_ATLAS_PATH') + '/KPspot/KPspot.atlas'
  openr, unit, atlasfile, /GET_LUN, /XDR

  Nwave = 0L
  readu, unit, Nwave
  Iatlas = fltarr(Nwave)
  latlas = Iatlas
  readu, unit, latlas, Iatlas
  free_lun, unit

  desired = where(latlas GE lambdaMin  AND latlas LE lambdamax, count)
  IF (count LE 0) THEN BEGIN
    print, "Requested wavelength range out of bounds of atlas"
    return
  ENDIF

  Iatlas = Iatlas[desired]
  latlas = latlas[desired]
END
