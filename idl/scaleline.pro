FUNCTION cc_kernel, s

  s2   = s * s
  s3   = s2 * s
  u = [-0.5*(s3 + s) + s2,       1.5*s3 - 2.5*s2 + 1.0, $
       -1.5*s3 + 2.0*s2 + 0.5*s, 0.5*(s3 - s2)]

  return, u
END

FUNCTION cubeconvol, f, x

  N = n_elements(f)
  IF (N LE 3) THEN return, interpolate(f, x)

  j = long(x)
  u = cc_kernel(x - j)
  g = 0.0

  ;; --- One-dimensional interpolation with cubic convolution -- ---- ;

  IF (j EQ 0) THEN BEGIN
    g = total(f[j:j+2] * u[1:3])
    g = g + (3.0*(f[j] - f[j+1]) + f[j+2]) * u[0]
  ENDIF ELSE IF (j EQ N-2) THEN BEGIN 
    g = total(f[j-1:j+1] * u[0:2])
    g = g + (3.0*(f[j] - f[j-1]) + f[j-2]) * u[3]
  ENDIF ELSE IF ((j GT 0)  AND  (j LT N-2)) THEN $
    g = total(f[j-1:j+2] * u) $
  ELSE IF (j EQ N-1) THEN $
    g = f[N-1]

  return, g
END

FUNCTION scaleline, f, Ncx

  Nx = n_elements(f)
  x  = float(Nx-1) / (Ncx-1) * findgen(Ncx)

  g = fltarr(Ncx)
  FOR i=0, Ncx-1 DO g[i] = cubeconvol(f, x[i])

  return, g
END
