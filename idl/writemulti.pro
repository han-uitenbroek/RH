PRO writemulti, model, geometry, HLTE=HLTE, FILENAME=filename, $
                NEW_DSCALE=new_dscale

  ;; Write a MULTI-format atmosphere from structure atmos
  ;; See READ1DATMOS for the format of the atmos structure used here.
  ;;
  ;; When NEW_DSCALE is set the atmosphere is interpolated to the new scale

  ;; Example:

  ;    a = read1datmos('VAL3C_52.atmos')
  ;    ndsc1 = alog(a.cmass[0]) + $
  ;	(alog(a.cmass[25]) - alog(a.cmass[0]))*dindgen(36)/35.0

  ;    ndsc2 = alog(a.cmass[25]) + $
  ;	(alog(a.cmass[51]) - alog(a.cmass[25]))*dindgen(46)/45.0

  ;    writemulti, a, NEW_DSCALE=exp([ndsc1, ndsc2[1:*]]), /HLTE


  CM_TO_M = 1.0D-2
  KM_TO_M = 1.0D+3

  IF (n_params(0) LT 1) THEN BEGIN
    print, "Usage: writemulti,  atmos, geometry, HLTE=HLTE, "
    print, "  FILENAME=filename, NEW_DSCALE=new_dscale"
    return
  ENDIF

  IF (keyword_set(NEW_DSCALE)) THEN BEGIN
    Ndep = n_elements(new_dscale)
    IF (Ndep GT 1) THEN BEGIN
      CASE (model.scale) OF
        'MASS_SCALE': tabinv, alog(model.cmass), alog(new_dscale), d_eff
        'TAU500_SCALE': tabinv, alog(model.tau500), alog(new_dscale), d_eff
        'GEOMETRIC_HEIGHT': tabinv, model.height, new_dscale, d_eff
      ENDCASE

      T = interpolate(model.T, d_eff)
      n_elec = exp(interpolate(alog(model.n_elec), d_eff))
      v = interpolate(model.v, d_eff)
      vturb = interpolate(model.vturb, d_eff)

      IF (n_elements(model.nH) GT 1 AND NOT keyword_set(HLTE)) THEN BEGIN
        nH = dblarr(Ndep, 6)
        FOR n=0, 5 DO $
         nH[*, n] = exp(interpolate(alog(model.nH[*, n]), d_eff, /CUBIC))
      ENDIF ELSE $
       nH = 0.0D+0

      atmos = {atmosID: model.atmosID,  gravitation: model.gravitation, $
               Ndep: long(Ndep),  scale: model.scale,  cmass: new_dscale, $
               T: T,  n_elec: n_elec,  v: v,  vturb: vturb,  nH: nH}
    ENDIF
  ENDIF ELSE $ 
   atmos = model

  IF (NOT keyword_set(FILENAME)) THEN $
   openw, unit, '/dev/tty', /GET_LUN, /MORE $
  ELSE $
   openw, unit, filename, /GET_LUN

  printf, unit, FORMAT='("*  Model written with writemulti", /"*")'
  printf, unit, FORMAT='(2X, A)', atmos.atmosID

  CASE (atmos.scale) OF
    'MASS_SCALE': BEGIN
      scaleString = "Mass scale"
      scaleLabel  = "lg Column mass"
    END
    'TAU500_SCALE': BEGIN
      scaleString = "Tau 500 scale"
      scaleLabel  = " lg Tau 500 nm"
    END
    'GEOMETRIC_SCALE': BEGIN
      scaleString = "Height scale"
      scaleLabel  = " Height [km]  "
    END
  ENDCASE
  printf, unit, FORMAT='(2X, A)', scaleString

  printf, unit, alog10(atmos.gravitation), atmos.Ndep, $
   FORMAT='("*", /"* lg g [cm s^-2]", /1X, F7.4, /"*", /"* Ndep", /I4, /"*")'

  printf, unit, scaleLabel, $
   FORMAT='("*", A, 5X, "Temperature", 8X, "Ne", 9X, "V", 14X, "Vturb")'

  FOR k=0, atmos.Ndep-1 DO BEGIN
    CASE (atmos.scale) OF
      'MASS_SCALE':      depthScale = alog10(atmos.cmass[k] / 10.0)
      'TAU500_SCALE':    depthScale = alog10(atmos.tau500[k])
      'GEOMETRIC_SCALE': depthScale = atmos.height[k]
    ENDCASE

    printf, unit, FORMAT='(E17.8, 4E15.6)', $
     depthScale, atmos.T[k], atmos.n_elec[k] * CM_TO_M^3, $
     atmos.v[k], atmos.vturb[k] / 1.0E3
  ENDFOR

  IF (n_elements(atmos.nH) LE 1 OR keyword_set(HLTE)) THEN BEGIN
    printf, unit, FORMAT='("*",/"* Hydrogen populations (LTE)")'
  ENDIF ELSE BEGIN
    printf, unit, FORMAT='("*",/"* Hydrogen populations")' 
    printf, unit, FORMAT=$
     '("*",4X,"nH[1]",8X,"nH[2]",8X,"nH[3]",8X,"nH[4]",8X,"nH[5]",8X,"np")'

    FOR k=0, atmos.Ndep-1 DO $
     printf, unit, FORMAT='(6E13.5)', atmos.nH[k, *] * CM_TO_M^3
  ENDELSE

  free_lun, unit
END

