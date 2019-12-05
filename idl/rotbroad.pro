; -------- file: -------------------------- rotbroad.pro ------------- ;

; -------- begin -------------------------- rotbroad.pro ------------- ;

FUNCTION rotbroad, lambda, intensity, xmu, v_rot, NPHI=Nphi

  ;; --- Computes rotationally broadened line profile for specific mu value. 

  ;;     Inputs:     lambda:   Wavelength
  ;;              intensity:   Intensity values
  ;;                    xmu:   value of mu = cos(theta)
  ;;                  v_rot:   rotational velocity [km/s]
  
  ;;     Keyword:      NPHI:   Number of azimuthal grid points

  ;;     Han Uitenbroek
  ;; --- Last modified: Tue May 16 12:58:46 2000 --

  CLIGHT  = 2.99792458E+08
  KM_TO_M = 1.0E+03

  IF (n_params(0) LT 4) THEN BEGIN
    print, 'Usage: ROTBROAD, lambda, intensity, xmu, v_rot, NPHI=Nphi'
    return, 0
  ENDIF

  IF (xmu GE 1.0) THEN return, intensity

  IF (NOT keyword_set(NPHI)) THEN Nphi = 21

  phi  = !PI * findgen(Nphi)/float(Nphi-1)
  wphi = fltarr(Nphi) + 1.0/float(Nphi-1)
  wphi[0] = 0.5*wphi[0]  &  wphi[Nphi-1] = 0.5*wphi[Nphi-1]

  I_broad = fltarr(n_elements(lambda))

  FOR n=0, Nphi-1 DO BEGIN
    vphi = v_rot * cos(phi[n]) * sqrt(1.0 - xmu^2)
    tabinv, lambda, lambda * (1.0 + vphi * KM_TO_M / CLIGHT), lambda_eff
    I_v = interpolate(intensity, lambda_eff, /CUBIC)
    I_broad = I_broad +  wphi[n] * I_v
  ENDFOR

  return, I_broad
END
; -------- end ---------------------------- rotbroad.pro ------------- ;

; -------- begin -------------------------- rotflux.pro -------------- ;

FUNCTION rotflux, lambda, spectrum, xmu, wmu, v_rot, NPHI=Nphi

  ;; --- Computes rotationally broadened flux profile.

  ;;     Inputs:     lambda:   Wavelength
  ;;               spectrum:   Intensity values as function of lambda and mu
  ;;                    xmu:   values of mu = cos(theta)
  ;;                    wmu:   integration weights
  ;;                  v_rot:   rotational velocity [km/s]
  
  ;;     Keyword:      NPHI:   Number of azimuthal grid points

  IF (n_params(0) LT 5) THEN BEGIN
    print, 'Usage: ROTFLUX, lambda, spectrum, xmu, wmu, v_rot, NPHI=Nphi'
    return, 0
  ENDIF

  Nrays = n_elements(xmu)
  flux  = fltarr(n_elements(lambda))

  FOR n=0, Nrays-1 DO BEGIN
    IF (NOT keyword_set(NPHI)) THEN $
     I_broad = rotbroad(lambda, spectrum[*, n], xmu[n], v_rot) $
    ELSE $
     I_broad = rotbroad(lambda, spectrum[*, n], xmu[n], v_rot, NPHI=Nphi)

   flux = flux + I_broad*xmu[n]*wmu[n]
  ENDFOR

  return, flux * 2*!PI
END
; -------- end ---------------------------- rotflux.pro -------------- ;
