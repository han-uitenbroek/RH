FUNCTION returnS, lambdaNo, rayNo, ACTIVE=active, BACKGROUND=background, Bp

@atmos.common
@spectrum.common
@opacity.common
@files.common

  ;; --- Returns the source function at wavelength index lambdaNo,
  ;;     and ray index rayNo --                        --------------- ;;

  IF (n_params() LT 2) THEN BEGIN
    print, "Usage: S = returnS(lambdaNo, rayNo [, Bp, /ACTIVE, /BACKGROUND])"
    return, 0
  ENDIF 
  IF (n_params() EQ 3) THEN $
   Bp = Planck(atmos.T, spectrum.lambda[lambdaNo], /HZ)

  result = openJ('J.dat')
  readJ, lambdaNo
  backgroundFile = 'background.dat'
  result = openOpacity('opacity.out')
  readOpacity, lambdaNo, rayNo
  free_lun, Junit, opacUnit, backgroundUnit

  IF (keyword_set(ACTIVE)) THEN BEGIN
    as_non_zero = where(chi_as GT 0.0, count)
    S_as = make_array(SIZE=size(chi_as), VALUE=0.0)
    IF (count GT 0) THEN $
     S_as[as_non_zero] = eta_as[as_non_zero] / chi_as[as_non_zero] $
    ELSE $
     S_as = eta_as / chi_as
    return, S_as
  ENDIF

  IF (keyword_set(BACKGROUND)) THEN $
    return, (eta_c + scatt*J) / chi_c

  return, (eta_as + eta_c + J*scatt) / (chi_as + chi_c)
END

