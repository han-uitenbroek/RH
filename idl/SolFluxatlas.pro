PRO SolFluxatlas, lambdamin, lambdamax, flux, lambda, NM=nm, $
                  IRRADIANCE=irradiance

  IF (n_params() LT 3) THEN BEGIN
    print, " Usage: SolFluxatlas, lambdamin, lambdamax, flux, lambda, NM=nm"
    return
  ENDIF

  atlasdir = getenv('RH_ATLAS_PATH') + '/SolarFlux/'
  files = file_search(atlasdir + 'lm????.dat', COUNT=Nfile)
  wavelength = fltarr(Nfile)
  FOR n=0, Nfile-1 DO BEGIN
    length = strlen(files[n])
    wavelength[n] = float(strmid(files[n], length-8, length-5))
  ENDFOR

  IF (lambdamin LT wavelength[0]) THEN BEGIN
    print, FORMAT='(" Warning: minimum wavelength ", F9.3,' + $
     '" below that of atlas", F9.3)', lambdamin, wavelength[0]
  ENDIF
  IF (lambdamax GT  wavelength[Nfile-1] + 4.0) THEN BEGIN
    print, FORMAT='(" Warning: maximum wavelength ", F9.3,' + $
     '" above that of atlas", F9.3)', lambdamax, wavelength[Nfile-1]
  ENDIF
  IF (lambdamin GE wavelength[Nfile-1] OR $
      lambdamax LE wavelength[0]) THEN BEGIN
    print, " Error wavelength range outside domain of atlas:"
    print, FORMAT='(3X, F9.3, " < lambda < ", F9.3)', $
     wavelength[0], wavelength[Nfile-1] + 4.0
    flux = 0.0
    return
  ENDIF

  N0 = 0
  WHILE (lambdamin GT wavelength[N0+1] AND N0 LT Nfile-1) DO N0 = N0 + 1
  N1 = Nfile-1
  WHILE (lambdamax LE wavelength[N1]) DO N1 = N1 - 1
  needed = N0 + indgen(N1-N0+1)

  Nlambda = 0L
  FOR n=0, n_elements(needed)-1 DO BEGIN
    openr, lun, files[needed[n]], /GET_LUN, /XDR
    readu, lun, Nlambda
    l0 = dblarr(Nlambda)
    r0 = fltarr(Nlambda)
    f0 = fltarr(Nlambda)
    readu, lun, l0, r0, f0
    free_lun, lun

    IF (n EQ 0) THEN BEGIN
      lambda = l0
      flux = f0
    ENDIF ELSE BEGIN
      Nlast = n_elements(lambda) - 1
      non_overlap = where(l0 GT lambda[Nlast])
      lambda = [lambda, l0[non_overlap]]
      flux   = [flux, f0[non_overlap]]
    ENDELSE
  ENDFOR
  index  = where(lambda GE lambdamin AND lambda LE lambdamax)
  lambda = lambda[index]
  flux   = flux[index]

  MICROWATT = 1.0E-6
  CM_TO_M   = 1.0E-2
  NM_TO_M   = 1.0E-9
  CLIGHT    = 2.99792458E+08

  IF (keyword_set(IRRADIANCE)) THEN BEGIN

    ;; Irradiance at Earth's orbit in W m^-2 nm^-1

    flux = MIRCOWATT * flux / CM_TO_M^2
    return
  ENDIF

  ;; The atlas gives relative flux and irradiance at the earth orbit
  ;; in units of \mu W cm^-2 nm^-1

  ;; Ratio of solar orbit radius to solar radius (both in km)

  R_RATIO   = 1.495985E+8 / 6.9598E+5
  conversion = MICROWATT * R_RATIO^2 / CM_TO_M^2

  ;; Flux at the solar surface in W m^-2 nm^-1

  flux = flux * conversion

  ;; Else flux at the surface in W m^-2 Hz^-1

  IF (NOT keyword_set(NM)) THEN $
   flux = (lambda * NM_TO_M)^2 / CLIGHT * (flux / NM_TO_M)
END
