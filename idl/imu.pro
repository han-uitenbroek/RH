; -------- file: -------------------------- imu.pro ------------------ ;

; -------- begin -------------------------- imu.pro ------------------ ;

FUNCTION imu, spectrum, thexmu, xmu

  IF (n_params(0) LT 3) THEN BEGIN
    print, 'Usage: IMU, spectrum, thexmu, xmu'
    return, 0
  ENDIF

  ;; --- Calculate intensities at arbitrary viewing angle. -- -------- ;

  IF ((xmu LT thexmu(0))  OR  (xmu GT 1.0)) THEN BEGIN
    print, thexmu(0), 1.0, FORMAT='("Improper value of viewing angle:", ' + $
       '/"  xmu must be between ", F5.3, " and ", F5.3)'
    return, 0
  ENDIF

  ;; --- Expects spectrum(Nrays, Nspect), xmu(Nspect) -- ------------- ;

  spsize = size(spectrum)
  Nspect = spsize(1)  &  Nrays = spsize(2)
  themu  = [thexmu, 1.0 + reverse(1.0 - thexmu)]

  i_mu = fltarr(Nspect)
  FOR ns=0,Nspect-1 DO BEGIN
    int = cmprss(spectrum(ns, *))
    i_mu(ns) = spline(themu, [int, reverse(int)], xmu, 10.0)
  ENDFOR

  return, i_mu
END
; -------- end ---------------------------- imu.pro ------------------ ;
