PRO make_B_file, atmosfile, B, GAMMA=gamma, CHI=chi, BX=Bx, BY=By, BZ=Bz, $
                 B_FILENAME=B_filename

  GAUSS_TO_TESLA = 1.0E-4

  ;; Input:
  ;;
  ;;     B -- Magnetic field strength [Gauss]
  ;; gamma -- Inclination of the magnetic field [degrees]
  ;;   chi -- Azimuth of the magnetic field [degrees]

  atmos = read1datmos(atmosfile)

  IF (atmos.ndep NE n_elements(B) AND NOT keyword_set(BX)) THEN BEGIN
    print, "Array B not of proper length: ", $
     n_elements(B), " instead of ", atmos.ndep
    return
  ENDIF
  IF (NOT keyword_set(GAMMA)) THEN gamma = dblarr(atmos.Ndep)
  IF (keyword_set(GAMMA) AND n_elements(gamma) NE atmos.Ndep) THEN BEGIN
    print, "Array keyword GAMMA not of proper length: ", $
     n_elements(gamma), " instead of ", atmos.ndep
    return
  ENDIF
  IF (NOT keyword_set(CHI)) THEN chi = dblarr(atmos.Ndep)
  IF (keyword_set(CHI) AND n_elements(chi) NE atmos.Ndep) THEN BEGIN
    print, "Array keyword chi not of proper length: ", $
     n_elements(chi), " instead of ", atmos.ndep
    return
  ENDIF

  IF (keyword_set(BX)) THEN BEGIN
    IF (n_elements(bx) NE atmos.Ndep) THEN BEGIN
      print, "Array keyword Bx not of proper length: ", $
       n_elements(bx), " instead of ", atmos.ndep
      return
    ENDIF
    B     = sqrt(Bx^2 + By^2 + Bz^2)
    chi   = atan(Bx, By)*!RADEG
    gamma = acos(Bz / B)*!RADEG
  ENDIF

  IF (NOT keyword_set(B_FILENAME)) THEN $
   B_filename = strmid(atmosfile, 0, strlen(atmosfile)-6) + ".B"

  openw, lun, B_filename, /GET_LUN, /XDR
  writeu, lun, double(B*GAUSS_TO_TESLA), $
   double(gamma/!RADEG), double(chi/!RADEG)
  free_lun, lun

END
