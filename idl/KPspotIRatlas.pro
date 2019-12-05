PRO KPspotIRatlas, lambdaMin, lambdaMax, spectrum, lambda, $
  COLD=cold, HOT=hot, PHOTOSPHERE=photosphere

  atlas_file = getenv('RH_ATLAS_PATH') + '/KPspotatl5//'
  IF (keyword_set(HOT)) THEN BEGIN
    atlas_file += 'sptr3ascii'
    Nwave       = 494217L
    dummy       = fltarr(3, Nwave)
  ENDIF ELSE IF (keyword_set(COLD)) THEN BEGIN
    atlas_file += 'sptr4ascii'
    Nwave       = 500326L
    dummy       = fltarr(3, Nwave)
  ENDIF ELSE IF (keyword_set(PHOTOSPHERE)) THEN BEGIN
    atlas_file += 'sptr2ascii'
    Nwave       = 524215L
    dummy       = fltarr(3, Nwave)
  ENDIF ELSE BEGIN
    atlas_file += 'sptr1ascii'
    Nwave       = 502250L
    dummy       = fltarr(4, Nwave)
  ENDELSE    


  openr, lun, /GET_LUN, atlas_file
  readf, lun, dummy
  free_lun, lun

  wavenumber   = reform(dummy[0, *])
  spectrum     = reform(dummy[1, *])  ;;; corrected spectrum
  transmission = reform(dummy[2, *])  ;;  telluric transmission

  lambda = vacuumtoair(1.0E7 / wavenumber)

  desired = where(lambda GE lambdaMin  AND lambda LE lambdamax, count)
  IF (count LE 0) THEN BEGIN
    print, "Requested wavelength range out of bounds of atlas"
    return
  ENDIF

  spectrum = spectrum[desired]
  lambda   = lambda[desired]
END

