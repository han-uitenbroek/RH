PRO Hawaiiatlas, lambdaMin, lambdaMax, SAMPLE=sample, $
                 Iatlas, latlas

  ; Reads selected wavelength range from the Hawaii UV Atlas
  ;
  ;  Allen, McAllister & Jefferies 1978, 
  ;  High Resolution Atlas of the Solar Spectrum form 2678 - 2931 \AA, 
  ;  Institute for Astronomy, University of Hawaii

  ; Parameters:
  ;   Input:    lambdaMin  -- Minimum wavelength [nm]
  ;             lambdaMax  -- Maximum wavelength [nm]

  ;   Output:   Iatlas     -- Integer array of intensities normalized to 10^4
  ;             latlas     -- Float array of corresponding wavelengths [nm]

  ; Keywords:
  ;             SAMPLE     -- Wavelength sampling interval [nm]

  cLight = 2.99792458E+08

  if (n_params(0) lt 3) then begin
    print, "Usage: atlas, lambdaMin, lambdaMax, SAMPLE=sample, " + $
     "Iatlas, [latlas,] LIMB=limb"
    print, "NOTE:  Wavelengths are in nanometers!!"
    return
  ENDIF
  IF (lambdaMin GE lambdaMax) THEN BEGIN
    print, "Minimum wavelength should be larger than maximum"
    return
  ENDIF

  atlasFile = getenv('RH_ATLAS_PATH') + '/Hawaii/hawaii_UV.dat'
  openr, unit, atlasFile, /GET_LUN

  little_endian = (is_big_endian()) ? 0 : 1

  Nspect = 0L
  readu, unit, Nspect
  IF (little_endian) THEN Nspect = swap_endian(Nspect)

  ;; --- Determine wavelength sacle from two known lines -- --------- ;;

  rec1 = 33963  &  lambda1 = 284.984           ;; CrII line
  rec2 = 35132  &  lambda2 = 285.567           ;; CrII/FeII doublet

  dlambda = (lambda2 - lambda1)/(rec2 - rec1)
  lambda0 = lambda1 - rec1*dlambda
  lambdaN = lambda2 + (Nspect - rec2)*dlambda

  IF (NOT keyword_set(SAMPLE)) THEN sample = dlambda
  sample = dlambda > sample

  IF ( (lambdaMin LT lambda0) OR (lambdaMin GT lambdaN) )THEN BEGIN
    print, FORMAT='("Minimum wavelength outside range: ", F5.1, ' + $
                  '" - ", F6.1, " [nm]")', lambda0, lambdaN
    return
  ENDIF
  IF ( (lambdaMax LT lambda0) OR (lambdaMax GT lambdaN) ) THEN BEGIN
    print, FORMAT='("Maximum wavelength outside range: ", F5.1, ' + $
                  '" - ", F6.1, " [nm]")', lambda0, lambdaN
    return
  ENDIF

  recMin = long((lambdaMin - lambda0) / dlambda)
  recMax = long((lambdaMax - lambda0) / dlambda)
  Nrecords = recMax - recMin + 1

  Iatlas = fltarr(Nrecords)
  point_lun, unit, (1 + recMin) * 4
  readu, unit, Iatlas
  free_lun, unit

  Nsample = fix(sample / dlambda)
  IF (Nsample GT 1) THEN $
   Iatlas = Iatlas( Nsample*lindgen(Nrecords/Nsample) )
  IF (little_endian) THEN Iatlas = swap_endian(Iatlas)

  IF (n_params(0) EQ 4) THEN BEGIN
    sample  = Nsample * dlambda
    lambda0 = lambda0 + recMin*dlambda
    latlas  = lambda0 + sample*findgen(n_elements(Iatlas)) - 0.017
  ENDIF
  ;; --- Convert from [erg cm^-2 s^-1 nm^-1 sr^-1] to
  ;;     [J m^-2 s^-1 Hz^-1 sr^-1] --                    ------------ ;;

  Iatlas = (1.0E-3 / cLight) * (latlas)^2 * 1.0E-9 * Iatlas

  return
end
