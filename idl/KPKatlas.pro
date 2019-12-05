PRO KPKatlas, lambdaMin, lambdaMax, SAMPLE=sample, $
              Iatlas, latlas

  ; Reads selected wavelength range from the Kohl, Parkinson, &Kurucz UV Atlas
  ;
  ; Parameters:
  ;   Input:    lambdaMin  -- Minimum wavelength [nm]
  ;             lambdaMax  -- Maximum wavelength [nm]

  ;   Output:   Iatlas     -- Integer array of intensities normalized to 10^4
  ;             latlas     -- Float array of corresponding wavelengths [nm]

  ; Keywords:
  ;             SAMPLE     -- Wavelength sampling interval [nm]

  cLight = 2.99792458E+08
  lambda0 = 224.200  &  lambdaN = 319.670
  dlambda = 0.0012

  if (n_params(0) lt 3) then begin
    print, "Usage: KPKatlas, lambdaMin, lambdaMax, SAMPLE=sample, " + $
     "Iatlas, [latlas,] LIMB=limb"
    print, "NOTE:  Wavelengths are in nanometers!!"
    return
  ENDIF
  IF (lambdaMin GE lambdaMax) THEN BEGIN
    print, "Minimum wavelength should be larger than maximum"
    return
  ENDIF

  atlasFile = getenv('RH_ATLAS_PATH') + '/KPK/KPK_atlas.dat'
  openr, unit, atlasFile, /GET_LUN

  little_endian = is_big_endian() ? 0 : 1

  Nspect = 0L
  readu, unit, Nspect
  IF (little_endian) THEN Nspect = swap_endian(Nspect)

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

  Iatlas = fltarr(Nspect)  &  latlas = Iatlas
  readu, unit, latlas, Iatlas
  free_lun, unit
  IF (little_endian) THEN BEGIN
    IF (n_params(0)EQ 4) THEN latlas = swap_endian(latlas)
    Iatlas = swap_endian(Iatlas)
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
end
