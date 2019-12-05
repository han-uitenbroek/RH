PRO sumeratlas, lambdaMin, lambdaMax, SAMPLE=sample, Iatlas, latlas, $
                NETWORK=network

  ; Reads selected wavelength range from the SUMER/SOHO UV Atlas
  ;
  ; Parameters:
  ;   Input:    lambdaMin  -- Minimum wavelength [nm]
  ;             lambdaMax  -- Maximum wavelength [nm]

  ;   Output:   Iatlas     -- Integer array of intensities (counts/s/pix)
  ;             latlas     -- Float array of corresponding wavelengths [nm]

  ; Keywords:
  ;             SAMPLE     -- Wavelength sampling interval [nm]

  if (n_params(0) lt 3) then begin
    print, "Usage: sumeratlas, lambdaMin, lambdaMax, SAMPLE=sample, " + $
     "Iatlas, [latlas] [, /NETWORK] "
    print, "NOTE:  Wavelengths are in nanometers!!"
    return
  ENDIF
  IF (lambdaMin GE lambdaMax) THEN BEGIN
    print, "Minimum wavelength should be larger than maximum"
    return
  ENDIF

  cLight  = 2.99792458E+08
  dlambda = 0.00439141

  IF (keyword_set(NETWORK)) THEN BEGIN
    atlasFile = getenv('RH_ATLAS_PATH') + '/SUMER/sumer_atlas_network.dat'
    openr, unit, atlasFile, /GET_LUN, /XDR
  ENDIF ELSE BEGIN
    atlasFile = getenv('RH_ATLAS_PATH') + '/SUMER/sumer_atlas_quiet.dat'
    openr, unit, atlasFile, /GET_LUN, /XDR
  ENDELSE

  Nspect = 0L
  readu, unit, Nspect
  Iatlas = fltarr(Nspect)  &  latlas = Iatlas
  readu, unit, latlas, Iatlas
  latlas = latlas / 5.0           ;; Wavelengths are given in first order
  free_lun, unit
  valid = where(Iatlas GT 0.0)
  Iatlas = Iatlas(valid)  &  latlas = latlas(valid)
  lambda0 = latlas(0)  &  lambdaN = latlas(N_elements(latlas) - 1)

  IF (NOT keyword_set(SAMPLE)) THEN sample = dlambda
  sample = sample > dlambda

  IF ( (lambdaMin LT lambda0) OR (lambdaMin GT lambdaN) )THEN BEGIN
    print, FORMAT='("Minimum wavelength outside range: ", F5.1, ' + $
     '" - ", F6.1, " [nm]")', lambda0, lambdaN
    lambdaMin = lambda0
  ENDIF
  IF ( (lambdaMax LT lambda0) OR (lambdaMax GT lambdaN) ) THEN BEGIN
    print, FORMAT='("Maximum wavelength outside range: ", F5.1, ' + $
     '" - ", F6.1, " [nm]")', lambda0, lambdaN
    lambdaMax = lambdaN
  ENDIF

  index = where((latlas GE lambdaMin)  AND (latlas LE lambdaMax), Nrecords)
  Iatlas = Iatlas(index)

  Nsample = fix(sample / dlambda)
  Iatlas = Iatlas(Nsample*lindgen(Nrecords/Nsample))
  IF (n_params(0) EQ 4) THEN BEGIN
    latlas = latlas(index)
    latlas = latlas(Nsample*lindgen(Nrecords/Nsample))
  ENDIF

  return
END
