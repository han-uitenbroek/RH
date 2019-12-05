PRO KPNOatlas, lambdaMin, lambdaMax, SAMPLE=sample, $
               Iatlas, latlas, LIMB=limb

  ; Reads selected wavelength range from the Kitt Peak Solar Atlas
  ;
  ; Parameters:
  ;   Input:    lambdaMin  -- Minimum wavelength [nm]
  ;             lambdaMax  -- Maximum wavelength [nm]

  ;   Output:   Iatlas     -- Integer array of intensities normalized to 10^4
  ;             latlas     -- Float array of corresponding wavelengths [nm]

  ; Keywords:
  ;             LIMB       -- If set intensities will be limb values (mu=0.2)
  ;             SAMPLE     -- Wavelength sampling interval [nm]

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

  atlasdir = getenv('RH_ATLAS_PATH') + '/KPNO/'
  if (keyword_set(limb)) then file = 'KPNO.Sun.Limb' else $
    file = 'KPNO.Sun.Center'
  openr, unit, atlasdir + file, /GET_LUN

  IF (NOT is_big_endian()) THEN $
   little_endian = 1 $
  ELSE $
   little_endian = 0

  interval = 0.1                ;; Wavelength interval [nm] spanning 1 record
  Nrec     = 500L               ;; Number of values in each 1 \AA interval
  dlambda  = interval / Nrec
  IF (NOT keyword_set(SAMPLE)) THEN sample = dlambda
  sample = dlambda > sample < 0.1

  Nlambda = 0
  readu, unit, Nlambda
  IF (little_endian) THEN Nlambda = swap_endian(Nlambda)

  lambda = intarr(Nlambda)
  readu, unit, lambda
  IF (little_endian) THEN lambda = swap_endian(lambda)
  lambda = float(lambda)/10.0
  lambda0 = lambda[0]
  lambdaN = lambda[Nlambda-1] + 0.1 - dlambda

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

  tabinv, lambda, [lambdaMin, lambdaMax], leff
  rec = fix(leff)
  offsetMin = fix((lambdaMin - lambda[rec[0]])/interval * Nrec)
  offsetMax = fix((lambdaMax - lambda[rec[1]])/interval * Nrec)
  Nrecords = Nrec * (rec[1] - rec[0]) + offsetMax - offsetMin

  Iatlas = intarr(Nrecords)
  point_lun, unit, (1 + Nlambda + rec[0]*Nrec + offsetMin) * 2
  readu, unit, Iatlas
  free_lun, unit

  Nsample = fix(sample / dlambda)
  IF (Nsample GT 1) THEN $
   Iatlas = Iatlas( Nsample*lindgen(Nrecords/Nsample) )
  IF (little_endian) THEN Iatlas = swap_endian(Iatlas)

  IF (n_params(0) EQ 4) THEN BEGIN
    sample = Nsample * dlambda
    lambda0 = lambda[rec[0]] + offsetMin*dlambda 
    latlas  = lambda0 + sample*findgen(n_elements(Iatlas))
  ENDIF
  Iatlas = Iatlas * 1.0E-4

  return
end
