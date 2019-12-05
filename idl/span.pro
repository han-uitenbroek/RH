FUNCTION span, x0, x1, Npoint, DOUBLE=double

  ;; Creates an array of Npoints equidistantly spaced between x0 and x1

  IF (n_params() LT 3) THEN BEGIN
    print, " Usage: x = span(x0, x1, Npoint[, /DOUBLE])"
    return, 0
  ENDIF

  IF (keyword_set(DOUBLE)) THEN $
   xi = dindgen(Npoint) $
  ELSE $
   xi = findgen(Npoint)

   return, x0 + (x1 - x0) * xi / (Npoint - 1) 
 END
