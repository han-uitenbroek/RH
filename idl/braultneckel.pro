PRO braultneckel, lambda_blue, lambda_red, Ilambda, lambda, Icont

  IF (n_params() LT 4) THEN BEGIN
    print, $
     "Usage: braultneckel, lambda_blue, lambda_red, Ilambda, lambda, Icont"
    return
  ENDIF
  IF (lambda_red LE lambda_blue) THEN BEGIN
    print, "LAMBDA_RED has to be larger than LAMBDA_BLUE: "
    print, "  LAMBDA_BLUE = ", lambda_blue, ",  LAMBDA_RED = ", lambda_red
    return
  ENDIF

  ANGSTROM_TO_NM = 0.1
  CLIGHT = 2.99792458E+08
  NM_TO_M = 1.0E-09
  CM_TO_M = 1.0E-02

  files  = file_search("$RH_ATLAS_PATH/LabsNeckel/file1[1-9]")
  files  = [files, file_search("$RH_ATLAS_PATH/LabsNeckel/file20")]

  Nlines = [ 156840L, 191426L, 147521L, 131128L, $
             111895L,  99141L,  96235L,  78737L,  65615L, 29424L ]

  lstart = [ 3290.0013,  4000.0041,  5000.0044,  6000.0026, $
             7000.0037,  8000.0006,  9000.0041, 10000.0080, 11000.0013, $
             12000.0063 ] * ANGSTROM_TO_NM

  lend   = [ 3999.9995,  4999.9984,  5999.9960,  6999.9948, $
             7999.9916,  8999.9948,  9999.9964, 10999.9873, 11999.9897, $
             12509.9824 ] * ANGSTROM_TO_NM

  IF (lambda_red LT lstart[0] OR $
      lambda_blue GT lend[n_elements(files)-1]) THEN BEGIN
    print, "Wavelengths outside atlas domain: ", lambda_blue, lambda_red, $
           ", atlas = [", lstart[0], "' ", lend[n_elements(files)-1], "]"
    return
  ENDIF

  tabinv, lstart, lambda_blue, index  
  nstart = fix(index)
  nend = nstart
  WHILE (lend[nend] LT lambda_red) DO nend++

  FOR n=nstart, nend DO BEGIN
    tmp = fltarr(3, Nlines[n])
    openr, lun, /GET_LUN, files[n]
    readf, lun, tmp
    free_lun, lun

    IF (n EQ nstart) THEN BEGIN
      lambda  = reform(tmp[0, *])
      Ilambda = reform(tmp[1, *])
      Icont   = reform(tmp[2, *])
    ENDIF ELSE BEGIN
      lambda  = [lambda, reform(tmp[0, *])]
      Ilambda = [Ilambda, reform(tmp[1, *])]
      Icont   = [Icont, reform(tmp[2, *])]
    ENDELSE
  ENDFOR

  lambda *= ANGSTROM_TO_NM
  tabinv, lambda, [lambda_blue, lambda_red], leff
  lleff = long(leff)
  lambda = lambda[lleff[0]:lleff[1]]

  Ilambda = Ilambda[lleff[0]:lleff[1]] / (CM_TO_M^2 * ANGSTROM_TO_NM)
  Ilambda /= (CLIGHT/(lambda^2 * NM_TO_M))

  IF (n_params() EQ 5) THEN BEGIN
    Icont = Icont[lleff[0]:lleff[1]] / (CM_TO_M^2 * ANGSTROM_TO_NM)
    Icont /= (CLIGHT/(lambda^2 * NM_TO_M))
  ENDIF
END
